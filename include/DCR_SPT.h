/*--------------------------------------------------------------------------*/
/*--------------------------- File DCR_SPT.h--------------------------------*/
/*--------------------------------------------------------------------------*/

/** @file
 * Header file for the class DCR_SPT, which implements
 * two heuristics for Single-Flow Single-Path (SFSP) DCR problems.
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
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef DCR_SPT_H
#define DCR_SPT_H

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "DCR.h"

/*--------------------------------------------------------------------------*/
/*--------------------------- CLASS DCR_SPT --------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// two Shortest-Path-Tree based heuristics for Single-Flow, Single-Path DCR
/** This class implements two heuristics for Single-Flow Single-Path (SFSP)
 *  DCR problems, both due to Orda and collaborators (the "Extended Routing
 *  Algorithm", ERA): for a Single-Flow DCR instance, only a *single* path
 *  from source to sink needs to be found, whose delay (comprised of a
 *  transmission term MTU / r_min, depending on the rate r_min reserved on
 *  every arc of the path, plus fixed per-arc/per-node delays, plus a term
 *  MTU / r_min * burst accounting for the flow's burstiness) must not
 *  exceed the flow's deadline, while minimizing the total routing cost.
 *
 *  Both heuristics restrict attention to a finite set of candidate values
 *  for r_min, namely the (distinct) capacities of the flow's links: for
 *  each such candidate rmin, only links whose capacity is >= rmin can be
 *  used (the "reduced graph"), and a Shortest Path Tree from the source is
 *  computed on it (see DCRheurERAI() and DCRheurERAH()) to check whether a
 *  feasible s-t path exists and, if so, at what cost; the best (lowest
 *  cost) feasible path found over all candidate values of rmin is
 *  returned. ERA-I (DCRheurERAI(), heur == 1) prices each arc of a
 *  candidate path using its own individual capacity as reserved rate,
 *  while ERA-H (DCRheurERAH(), heur == 2) additionally optimizes, for each
 *  candidate path, a single common reserved rate r0 shared by all its arcs
 *  (chosen as small as the deadline allows), typically yielding cheaper
 *  (and only single-path, mixed-integer) solutions. Since only single-path
 *  solutions are considered, the class cannot represent or solve the
 *  continuous (multi-path) relaxation of the DCR problem: DCRgetPSolNPaths()
 *  and DCRgetPSolPath() are therefore not supported. */

class DCR_SPT : public DCR
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
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/
 /** @name Constructors
    @{ */

 DCR_SPT( void );

 /**< Constructor: gives some default values to all the data
   members of the class.*/

 /** @} */
/*--------------------------------------------------------------------------*/
/*----------------------------------READ DATA ------------------------------*/
/*--------------------------------------------------------------------------*/

 /** @name Loading the data of the problem
    @{ */

 virtual void DCRloadProblem( int nnodes , int nlinks , int nflows ,
                              DCRFlow * flows , DCRLink * links ,
                              DCRNode * nodes , double MTU ,
                              DCRDelay deltype );

 /** @} */

/*--------------------------------------------------------------------------*/
/*----------------------------------SOLVE-----------------------------------*/
/*--------------------------------------------------------------------------*/
 /** @name Solving the problem
    @{ */

 virtual DCRStatus DCRsolve( void );

 /**< Solves the problem using two heuristics. */

 /** @} */

/*--------------------------------------------------------------------------*/
/*----------------------------------GET RESULTS-----------------------------*/
/*--------------------------------------------------------------------------*/

 /** @name Reading solver output
  @{ */

 virtual double DCRgetObj( void );

 virtual int DCRgetUBSolNPaths( int k );

 /**< Returns 1 if a (mixed-)integer solution exists, since we only
  *  solve *single-path* instances (see UBsol).
  * \param k flow index
  */

 virtual int DCRgetUBSolPath( int k , int p , int * X , double * R );

 /**< Returns the number of hops in path 0 of flow k in the current
  * (mixed-)integer solution,
  * if any (see UBsol). It only works for p=0, since we only solve
  * *single-path* instances.
  * X is an array of int of size (n-1) allocated by the user to store
  * the indices of the links in path p. R is an array of double
  * of size (n-1) allocated by the used to store the corresponding rates.
  * \param k flow index
  * \param p path index (only 0)
  * \param X pointer to an array of int of size (numNodes-1)
  * \param R pointer to an array of double of size (numNodes-1) */

 virtual int DCRgetPSolNPaths( int k );

 /**< This function is not supported because heuristics
  * cannot solve the continuous relaxation of the problem.*/

 virtual void DCRgetPSolPath( int k , int p , double * X , double * R );

 /**< This function is not supported because heuristics
  * cannot solve the continuous relaxation of the problem.*/

 /** @} */

/*--------------------------------------------------------------------------*/
/*----------------------------------SETTERS---------------------------------*/
/*--------------------------------------------------------------------------*/
 /** @name Other initializations
    @{ */

 virtual void DCRsetHeur( char );
 /**< Sets the type of heuristic to use.
  * This method does not appear in the base-class interface.*/

 /** @} */

/*--------------------------------------------------------------------------*/
/*----------------------------------MODIFIERS-------------------------------*/
/*--------------------------------------------------------------------------*/

 /** @name Changing the data of the problem
    @{ */

 virtual void DCRcloseArcs( int * whch , int na );
 /**< Closes arcs of the network, for all flows.
  * \param whch indices of the arcs to be closed*/

 virtual void DCRcloseArcs( int k , int * whch , int na );
 /**< Closes arcs of the network for a specific flow.
  * \param which indices of the arcs to be closed
  * \param k index of the flow to which the process applies*/

 /** @} */

/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/
 /** @name Destructor
    @{ */

 ~DCR_SPT();

 /**< Frees up dinamically allocated memory */
 /** @} */

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 /** @name Private part of the class
    @{ */

private:
 //network data
 DCRNode * Nodes; ///< array of numNodes node data (node delays)
 DCRLink * Links; ///< array of numLinks link data (capacity, speed,
                  /// delay, cost)
 DCRFlow * Flows; ///< the (single) flow to be routed, Flows[ 0 ]

 //can be either 1=ERA-I or 2=ERA-H
 int heur; ///< which heuristic to use, set by DCRsetHeur()

 //solution info
 int * Xsol;    ///< indices (in Links) of the arcs of the best path
                /// found so far, in reverse (sink-to-source) order
 double * Rsol; ///< reserved rate on each arc of Xsol, same order
 double objval; ///< cost of the best feasible path found so far
 int nhops;     ///< number of arcs (hops) of the best path found

 //private methods

 /// index (in Links) of the arc from "from" to "to", -1 if there is none
 int DCRgetLink( int from , int to );

 /// runs the ERA-I heuristic (see DCR_SPT.cpp for the details)
 void DCRheurERAI();

 /// runs the ERA-H heuristic (see DCR_SPT.cpp for the details)
 void DCRheurERAH();

 //memory management methods

 /// deep-copies flows/links/nodes into the Flows/Links/Nodes members
 void copyDataArrays( DCRFlow * flows , DCRLink * links , DCRNode * nodes );

 /// releases all dynamically allocated memory
 void clean_up( void );
 //solution management methods

 /** @} */
 };

#endif

/*--------------------------------------------------------------------------*/
/*-------------------------- End File DCR_SPT.h ----------------------------*/
/*--------------------------------------------------------------------------*/
