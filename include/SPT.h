/*--------------------------------------------------------------------------*/
/*---------------------------- File SPT.h ----------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class SPT, which computes a Shortest Path Tree from a
 * single source on a directed graph.
 *
 * SPT is the elementary building block used by DCRLagrangianSolver to
 * evaluate the Lagrangian relaxation of the (single-flow) Delay-Constrained
 * Routing problem: for a given value of the Lagrangian multiplier, the
 * subproblem reduces to finding a shortest path (in terms of "reduced"
 * arc costs) from the flow's source to every other node, which is exactly
 * what this class does.
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

/*<SPT for directed graphs with a single source*/

/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef SPT_H
#define SPT_H

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/
#include <limits>

#include <vector>

/*--------------------------------------------------------------------------*/
/*------------------------------- CLASS SPT --------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// Shortest Path Tree from a single source, directed graph, general costs
/** SPT computes a Shortest Path Tree from a given source node s to all
 * other nodes of a directed graph with arbitrary (possibly negative) real
 * arc costs, and in particular the shortest s-t path to a given sink node
 * t. The algorithm implemented (see Solve()) is the classical FIFO-queue
 * label-correcting method (a variant of the Bellman-Ford algorithm often
 * attributed to D'Esopo and Pape): distance labels are relaxed arc by arc,
 * a node is (re-)enqueued in a FIFO queue whenever its label improves, and
 * it is not re-inserted while it is already in the queue; this correctly
 * handles negative arc costs (as needed when SPT is used, by
 * DCRLagrangianSolver, to solve a Lagrangian-relaxed shortest path
 * subproblem where "reduced costs" can be negative), unlike Dijkstra's
 * label-*setting* algorithm which would require nonnegative costs. The
 * method does not explicitly detect negative-cost cycles: if the graph
 * (implicitly) becomes disconnected between the source and some node
 * needed to reconstruct the path -- which is how an unbounded instance
 * manifests itself here -- getStatus() reports Error. */

class SPT
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

public:
/*--------------------------------------------------------------------------*/
/*--------------------------- Inf() and Eps() ------------------------------*/
/*--------------------------------------------------------------------------*/
 /** Very small class to simplify extracting the "+ infinity" value for a
    basic type; just use Inf<type>(). */

 template < typename T > class Inf
 {
 public:
  Inf() {}
  operator T() { return( std::numeric_limits< T >::max() ); }
 };

/*--------------------------------------------------------------------------*/
 /** Very small class to simplify extracting the "machine epsilon" for a
    basic type; just use Eps<type>(). */

 template < typename T > class Eps
 {
 public:
  Eps() {}
  operator T() { return( std::numeric_limits< T >::epsilon() ); }
  };

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
 /** @name Public types
  *  @{ */

 /// status of the Shortest Path computation
 enum Status {
  OK ,         ///<  solver found a feasible solution
  Infeasible , ///<  problem infeasible
  Unbounded ,  ///<  problem unbounded
  Error       ///<  the graph is disconnected between s and some node
              /// needed to reconstruct the s-t path
  };

/*--------------------------------------------------------------------------*/

 /// a simplified arc for SPT: only cost, start and end node
 struct SPTLink
 {
  int startnode; ///< index of the tail node of the arc
  int endnode;   ///< index of the head node of the arc
  double cost;   ///< (possibly negative) cost/length of the arc

  double rstar; // not very useful in SP, but avoids a double check
                // on input and output for the Lagrangian solution.
  };

 /** @} ---------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/
 /** @name Constructor
  *  @{ */

 /// constructor of SPT: initializes all data members to "empty"

 SPT( void );

 /** @} ---------------------------------------------------------------------*/
/*---------------------------------- METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/
 /** @name Loading the data of the problem
  *  @{ */

 /// loads a new graph, source and sink node
 /** Loads a directed graph with nnodes nodes and nlinks arcs (links,
  * given as a vector of SPTLink), and sets sourcenode/sinknode as the
  * source/sink for the Shortest Path computation performed by Solve().
  * Any previously loaded instance is discarded. */

 virtual void LoadProblem( int nnodes , int nlinks ,
                           std::vector< SPTLink > links , int sourcenode ,
                           int sinknode );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// updates the arc costs, keeping the same topology
 /** Replaces the cost of every arc with the corresponding entry of
  * links (the topology, i.e., startnode/endnode of each arc, must be
  * unchanged), clearing any previously computed solution; meant to be
  * used to solve a sequence of instances differing only in the arc
  * costs (e.g., across Lagrangian iterations), without the overhead of
  * reloading the whole topology. */

 virtual void updCosts( std::vector< SPTLink > links );

 /** @} ---------------------------------------------------------------------*/
/*----------------------------------SOLVE-----------------------------------*/
/*--------------------------------------------------------------------------*/
 /** @name Solving the problem
  *  @{ */

 /// computes the Shortest Path Tree from the source node
 /** Runs a FIFO-queue label-correcting (Bellman-Ford-type) algorithm to
  * compute shortest distances from the source node to every other node,
  * and reconstructs the shortest s-t path (see getPath()). See the
  * general notes of the class for the details of the algorithm; sets
  * getStatus() to Error if the sink cannot be reached backwards from
  * itself down to the source (a disconnected/unbounded instance). */

 virtual void Solve();

 /** @} ---------------------------------------------------------------------*/
/*---------------------------------GET RESULTS------------------------------*/
/*--------------------------------------------------------------------------*/
 /** @name Reading the results
  *  @{ */

 /// returns the number of arcs (hops) of the optimal s-t path

 virtual int getNHops( void );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns the status of the last Solve() call

 virtual Status getStatus( void );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns the optimal s-t path as a list of arc indices
 /** Returns the sequence of arc indices (into the Links loaded by
  * LoadProblem()/updCosts()) forming the shortest path from the source
  * to the sink, in source-to-sink order. Note that if SPT is used on a
  * "reduced" subgraph (a subset of arcs of some larger graph), these
  * indices refer to the reduced graph and do *not* in general coincide
  * with the corresponding arc indices in the original, full graph. */

 virtual std::vector< int > getPath( void );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns the dual solution of SP, i.e., the shortest-distance labels
 /** Returns, for every node, its shortest distance from the source (the
  * dual variables/node potentials of the Shortest Path LP). */

 virtual std::vector< double > getLabel( void );

 /** @} ---------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/
 /** @name Destructor
  *  @{ */

 /// destructor of SPT: releases all dynamically allocated memory

 ~SPT();

 /** @} ---------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

private:
 std::vector< SPTLink > Links;  ///< the arcs of the graph
 int numNodes;                  ///< number of nodes of the graph
 int numLinks;                  ///< number of arcs of the graph
 int nhops;                     ///< number of arcs of the optimal s-t path
 int s;                         //sourceindex.
 int t;                         //sinkindex.
 std::vector< int > Sol;        //indices of the s-t path that solves SP.
 std::vector< double > DualSol; //dual solution of SP (labels).
 Status stat; //to check for the presence of negative-cost cycles.


/*--------------------------------------------------------------------------*/
 /// deep-copies links into the Links data member
 void copyDataArrays( std::vector< SPTLink > links );

/*--------------------------------------------------------------------------*/
 /// returns the index (in Links) of the arc with tail = from, head = to,
 /// or sets stat = Error and returns -1 if no such arc exists
 int getLink( int from , int to );

 /// clears all container data members (Links, Sol, DualSol)
 void clean_up( void );

 }; // end( class( SPT ) )

/*--------------------------------------------------------------------------*/

#endif

/*--------------------------------------------------------------------------*/
/*----------------------------- End File SPT.h -----------------------------*/
/*--------------------------------------------------------------------------*/
