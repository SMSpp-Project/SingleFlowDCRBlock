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


/*--------------------------------------------------------------------------*/
/*--------------------- IMPLEMENTATION OF SPT-----------------------*/
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

/// computes the Shortest Path Tree from the source node s
/** FIFO-queue label-correcting algorithm (a variant of Bellman-Ford,
 * sometimes attributed to D'Esopo-Pape): distance[] is initialized to
 * +Infinity everywhere but the source (distance[ s ] = 0), and nodes
 * are processed out of a FIFO queue Q[] (initially containing only s):
 * when node h is extracted, every arc leaving h is relaxed, and
 * whenever the relaxation improves the distance label of the head node
 * "next", that node is appended to the tail of the queue -- unless it
 * is already the head, the tail, or somewhere in the middle of it
 * (Q[ next ] == -1 tests exactly this), in which case it is left where
 * it is and will naturally be reprocessed with the updated label. Since
 * a node can be enqueued several times (which is what allows the
 * algorithm, unlike Dijkstra, to correctly handle negative arc costs),
 * the loop terminates only when the queue becomes empty, i.e., no
 * further improvement is possible; this implicitly assumes the graph
 * has no negative-cost cycle reachable from s (which is not checked
 * for directly: an infinite loop would ensue in that case).
 *
 * Once distances have stabilized, the shortest s-t path is
 * reconstructed backwards from t by repeatedly following previous[]
 * (the predecessor labels set while relaxing) until s is reached; if at
 * any point the arc from previous[ i ] to i cannot be found (getLink()
 * returns -1 and sets stat = Error), the reconstruction stops early,
 * signalling that t is not reachable from s. The resulting arc indices
 * are stored, in source-to-sink order, in Sol (see getPath()), and the
 * final node distances in DualSol (see getLabel()). */

void SPT::Solve()
{
 int i , j , h;
 //int maxIter = pow(numNodes,2);
 // int iterat=0;  //if it exceeds n^2, then there will be at least one
 // negative-cost cycle (?)
 double * distance = new double[ numNodes ];
 int * previous = new int[ numNodes ];
 int * Q = new int[ numNodes ]; //FIFO queue
 double dist , inf = Inf< double >();
 int next;
 int HEAD; //FIFO head
 int TAIL; //FIFO tail

 //initialization of the structures
 for( i = 0 ; i < numNodes ; i++ ) {
  distance[ i ] = inf;
  previous[ i ] = -1;
  Q[ i ] = -1;
  }
 distance[ s ] = 0;
 HEAD = s;
 TAIL = s;

 //solve SHORTEST PATH
 //////////////////////////////////////
 while( HEAD != -1 ) {
  //extract node from queue Q
  h = HEAD;
  HEAD = Q[ HEAD ];
  Q[ h ] = -1;
  if( HEAD == -1 )
   TAIL = -1;

  //consider all neighbours of h
  for( j = 0 ; j < numLinks ; j++ ) {
   if( Links[ j ].startnode == h ) {
    next = Links[ j ].endnode;
    dist = Links[ j ].cost + distance[ h ];

    if( dist < distance[ next ] ) {
     distance[ next ] = dist;
     previous[ next ] = h;

     //insert node in queue Q
     if( HEAD == -1 ) {
      HEAD = next;
      TAIL = next;
      }
     else if( ( HEAD != next ) && ( TAIL != next ) &&
                ( Q[ next ] == -1 ) ) {
      Q[ TAIL ] = next;
      TAIL = next;
      }

     } //if (distance label update)

    } //if (neighbour arc)

   } //for (all arcs)
  } //while (head != -1)
 //Set up the variables for the get methods.

 ///  /////is this ok? or new[nnodes] and ignore what was done in load pr, or
 /// delete and new??
 DualSol.resize( numNodes );
 //cout<<"total distance in SP = "<<distance[t]<<endl;

 for( j = 0 ; j < numNodes ; j++ )
  DualSol[ j ] = distance[ j ]; //is the for loop needed?

 i = t;
 int * indices = new int[ numNodes ];
 int index;
 nhops = 0; //needed for subsequent calls

 while( i != s && stat == OK ) {
  index = getLink(
   previous[ i ] ,
   i ); //if no link is found, it concludes that the graph is disconnected
  indices[ nhops ] =
   index; //copies the indices of the optimal path in reverse
  i = previous[ i ];
  nhops++; //saves the number of these indices
  }

 if( stat == OK ) {
  Sol.resize( nhops );

  for( j = 0 ; j < nhops ; j++ ) {
   Sol[ nhops - j - 1 ] =
    indices[ j ]; //takes the optimal path, putting it back in s-t order.
   }
  }
 else {
  // the sink is unreachable from the source in the current graph (or the
  // backtrack could not complete): leave no stale path behind for
  // getPath()/getNHops() to hand out. Without this, a caller that does
  // not check getStatus() (or one who does check it but only nhops/path
  // afterwards) would silently reuse whatever a previous, unrelated, and
  // possibly differently-sized successful Solve() left in Sol, mixing
  // arcs from two different solves into one bogus "path"
  Sol.clear();
  nhops = 0;
  }

 delete[] indices;
 delete[] previous;
 delete[] distance;
 delete[] Q;
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

// deep-copies links (startnode/endnode/cost) into the Links data member;
// numLinks must already be set

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
 }

/*--------------------------------------------------------------------------*/

// returns the index (in Links) of the arc with tail == from, head == to;
// if no such arc exists (the graph is disconnected there, which for the
// caller of Solve() typically means it was over-restricted, e.g. to a
// reduced graph, when trying to reconstruct the s-t path) sets
// stat = Error and returns -1

int SPT::getLink( int from , int to )
{
 int i;

 for( i = 0 ; i < numLinks ; i++ ) {
  if( Links[ i ].startnode == from && Links[ i ].endnode == to )
   return i;
  }
 //cout<<"rmin too large: the graph has become disconnected\n";
 stat = Error;
 return -1;
 }


/*--------------------------------------------------------------------------*/

// clears all container data members (Links, Sol, DualSol)

void SPT::clean_up()
{
 //cout<<"Links"<<endl;
 Links.clear();
 //cout<<"Sol"<<endl;
 Sol.clear();
 //cout<<"DualSol"<<endl;
 DualSol.clear();
 }


/*--------------------------------------------------------------------------*/
/*------------------------------ End SPT.cpp -------------------------------*/
/*--------------------------------------------------------------------------*/
