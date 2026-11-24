/*--------------------------------------------------------------------------*/
/*------------------------------ File SPT.cpp ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the SPT class.
 *
 * \version 1.00
 *
 * \date April - 2013
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Laura Galli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Luca Mencarelli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Enrico Sorbera \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni, Laura Galli, Luca Mencarelli,
 *                      Enrico Sorbera
 */

/*--------------------------------------------------------------------------*/
/*--------------------------- IMPLEMENTATION -------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "SPT.h"

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <queue>


/*--------------------------------------------------------------------------*/
/*--------------------- IMPLEMENTATION OF SPT-------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------------- PUBLIC METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*----------------------------- CONSTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

// gives default ("empty") values to all data members; LoadProblem() must
// be called before Solve() can do anything useful

SPT::SPT()
{
 Links.resize( 1 );
 nhops = 0;
 numNodes = numLinks = 0;
 s = 0;
 t = 0;
 stat = OK;
 Sol.resize( 1 ); //serve?
 DualSol.resize( 1 );
 }


/*--------------------------------------------------------------------------*/
/*----------------------------------LOAD DATA ------------------------------*/
/*--------------------------------------------------------------------------*/

// loads a brand new graph, discarding any previous instance, and sets the
// source/sink nodes for the Shortest Path computation

void SPT::LoadProblem( int nnodes , int nlinks , std::vector< SPTLink > links ,
                       int sourcenode , int sinknode )
{

 clean_up();

 stat = OK;
 s = sourcenode;
 t = sinknode;
 numNodes = nnodes;
 numLinks = nlinks;
 Sol.clear(); DualSol.clear(); nhops = 0;
 if(nnodes <= 0 || nlinks < 0 || links.size() != std::size_t(nlinks) ||
    s < 0 || s >= nnodes || t < 0 || t >= nnodes) { stat = Error; return; }
 for(const auto & e : links)
  if(e.startnode < 0 || e.startnode >= nnodes || e.endnode < 0 || e.endnode >= nnodes) {
   stat = Error; return;
  }
 copyDataArrays( links );
 }

/*--------------------------------------------------------------------------*/

// updates only the arc costs (topology, source and sink are unchanged)
// and clears any previously computed solution

void SPT::updCosts( std::vector< SPTLink > links )
{
 int i;

 Sol.clear();
 DualSol.clear();

 //delete[] Sol;
 //delete[] DualSol;
 //Links = new SPT::SPTLink[numLinks];

 nhops = 0;
 if(links.size() != Links.size()) { Links.clear(); stat = Error; return; }
 stat = OK;

 for( i = 0 ; i < numLinks ; i++ )
  Links[ i ].cost = links[ i ].cost;
 }

/*--------------------------------------------------------------------------*/
/*--------------------------- GET METHODS-----------------------------------*/
/*--------------------------------------------------------------------------*/

/// returns the optimal s-t path, as arc indices in source-to-sink order

std::vector< int > SPT::getPath( void )
{
 return( Sol );
 }

/*--------------------------------------------------------------------------*/

/// returns the status of the last Solve() call

SPT::Status SPT::getStatus()
{
 return( stat );
 }

/*--------------------------------------------------------------------------*/

/// returns the dual solution (shortest-distance label of every node)

std::vector< double > SPT::getLabel()
{
 return( DualSol );
 }

/*--------------------------------------------------------------------------*/

/// returns the number of arcs (hops) of the optimal s-t path

int SPT::getNHops()
{
 return( nhops );
 }

/*--------------------------------------------------------------------------*/
/*----------------------------------SOLVE-----------------------------------*/
/*--------------------------------------------------------------------------*/

/// FIFO label correction, with a bounded Bellman-Ford fallback if a
/// predecessor walk reaches n arcs. In the fallback, after k rounds labels describe
/// walks of at most k arcs. An improvement in round n signals a reachable
/// negative cycle (or uncertifiable floating-point arithmetic). Unreachable
/// components are never relaxed. No tolerance masks small negative costs.
void SPT::Solve()
{
 Sol.clear();
 DualSol.clear();
 nhops = 0;
 stat = Error;
 try {
  checkDeadline();
  if(numNodes <= 0 || s < 0 || s >= numNodes || t < 0 || t >= numNodes ||
     numLinks < 0 || Links.size() != std::size_t(numLinks)) return;
  for(const auto & e : Links)
   if(e.startnode < 0 || e.startnode >= numNodes ||
      e.endnode < 0 || e.endnode >= numNodes || !std::isfinite(e.cost)) return;

  const double inf = Inf<double>();
  std::vector<double> distance(numNodes,inf);
  std::vector<int> previous(numNodes,-1), previousArc(numNodes,-1);
  distance[s] = 0;
  std::queue<int> pending;
  std::vector<bool> queued(numNodes,false);
  std::vector<int> depth(numNodes,0);
  pending.push(s); queued[s] = true;
  bool fallback = false;
  while(!pending.empty() && !fallback) {
   checkDeadline();
   const int u = pending.front(); pending.pop(); queued[u] = false;
   unsigned scanned = 0;
   for(int j : adjOut[u]) {
    // Keep deadline checks bounded even for a very high-degree vertex,
    // without reading the clock on every individual arc of small graphs.
    if((scanned++ & 255) == 0) checkDeadline();
    const auto & e = Links[j];
    const double candidate = distance[u] + e.cost;
    if(!std::isfinite(candidate)) return;
    if(candidate < distance[e.endnode]) {
     distance[e.endnode] = candidate;
     previous[e.endnode] = u;
     previousArc[e.endnode] = j;
     depth[e.endnode] = depth[u] + 1;
     if(depth[e.endnode] >= numNodes) { fallback = true; break; }
     if(!queued[e.endnode]) {pending.push(e.endnode); queued[e.endnode] = true;}
    }
   }
  }
  // Depth is only a trigger, never a negative-cycle certificate. Confirm
  // with n synchronous Bellman-Ford rounds to avoid false queue-count alarms.
  if(fallback) {
   distance.assign(numNodes,inf); distance[s] = 0;
   previous.assign(numNodes,-1); previousArc.assign(numNodes,-1);
  }
  for(int round = 0; fallback && round < numNodes; ++round) {
   checkDeadline();
   auto nextDistance = distance;
   bool changed = false;
   for(int j = 0; j < numLinks; ++j) {
    if((j & 255) == 0) checkDeadline();
    const auto & e = Links[j];
    if(distance[e.startnode] == inf) continue;
    const double candidate = distance[e.startnode] + e.cost;
    if(!std::isfinite(candidate)) return;
    if(candidate < nextDistance[e.endnode]) {
     nextDistance[e.endnode] = candidate;
     previous[e.endnode] = e.startnode;
     previousArc[e.endnode] = j;
     changed = true;
    }
   }
   if(!changed) break;
   if(round == numNodes - 1) { stat = Unbounded; return; }
   distance.swap(nextDistance);
  }
  if(distance[t] == inf) { stat = Infeasible; return; }

  std::vector<bool> seen(numNodes,false);
  std::vector<int> reverse;
  for(int v = t; v != s;) {
   checkDeadline();
   if(v < 0 || v >= numNodes || seen[v] || reverse.size() >= std::size_t(numNodes))
    return;
   seen[v] = true;
   const int arc = previousArc[v], pred = previous[v];
   if(arc < 0 || arc >= numLinks || pred < 0 || pred >= numNodes ||
      Links[arc].startnode != pred || Links[arc].endnode != v) return;
   reverse.push_back(arc);
   v = pred;
  }
  Sol.assign(reverse.rbegin(),reverse.rend());
  nhops = int(Sol.size());
  DualSol = std::move(distance);
  stat = OK;
 } catch(const DCRDeadline::Expired &) {
  stat = Interrupted;
  throw;
 }
}


/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

// releases all dynamically allocated (container) memory

SPT::~SPT()
{
 clean_up();
 }

/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

// deep-copies links (startnode/endnode/cost) into the Links data member
// and rebuilds the out-arc adjacency list adjOut from the topology alone;
// numLinks/numNodes must already be set

void SPT::copyDataArrays( std::vector< SPTLink > links )
{
 int i , j;

 Links.resize( numLinks );

 for( j = 0 ; j < numLinks ; j++ ) {
  Links[ j ].startnode = links[ j ].startnode;
  Links[ j ].endnode = links[ j ].endnode;
  Links[ j ].cost = links[ j ].cost;
  //no need to copy Link.rstar at this level.
  }

 adjOut.clear();
 adjOut.resize( numNodes );
 for( j = 0 ; j < numLinks ; j++ )
  adjOut[ Links[ j ].startnode ].push_back( j );
 }

/*--------------------------------------------------------------------------*/

// clears all container data members (Links, Sol, DualSol, adjOut)

void SPT::clean_up()
{
 //cout<<"Links"<<endl;
 Links.clear();
 //cout<<"Sol"<<endl;
 Sol.clear();
 //cout<<"DualSol"<<endl;
 DualSol.clear();
 adjOut.clear();
 }


/*--------------------------------------------------------------------------*/
/*------------------------------ End SPT.cpp -------------------------------*/
/*--------------------------------------------------------------------------*/
