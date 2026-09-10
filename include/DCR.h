/*--------------------------------------------------------------------------*/
/*--------------------------- File DCR.h -----------------------------------*/
/*--------------------------------------------------------------------------*/

/** @file
 * Header file for the abstract (pure virtual) base class DCR, which
 * defines a standard interface for solvers of Delay Constrained Routing
 * (DCR) problems.
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

#ifndef DCR_H
#define DCR_H


/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <chrono>
#include <exception>
#include <ostream>

/*--------------------------------------------------------------------------*/
/*------------------------- CLASS DCRtimer ---------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// minimal wall-clock timer for the DCR solvers
/** Very small class measuring the time spent by a DCR solver: Start() and
 * Stop() bracket one run of the algorithm, Read() returns the total time
 * of all the runs since the timer was constructed or last ReSet(), the
 * current one included if the clock is ticking. The time is the wall-clock
 * one, since the split between user and system time is not available in a
 * portable way. */

class DCRtimer
{

/*--------------------------------------------------------------------------*/

public:

/*--------------------------------------------------------------------------*/

 /// constructor: the clock is not ticking and no time is accumulated

 DCRtimer( void ) : ticking( false ) , elapsed( 0 ) { }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// starts the clock, if it is not ticking already

 void Start( void ) {
  if( ! ticking ) {
   ticking = true;
   last = clock_type::now();
   }
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// stops the clock, accumulating the time of the run that just ended

 void Stop( void ) {
  if( ticking ) {
   elapsed += since_last();
   ticking = false;
   }
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// stops the clock and forgets all the time accumulated so far

 void ReSet( void ) {
  ticking = false;
  elapsed = 0;
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// total time, in seconds, accumulated so far
 /** Returns the total time, in seconds, of all the runs of the clock since
  * the last ReSet(); if the clock is ticking, the current run is counted
  * in without stopping it. */

 double Read( void ) const {
  return( elapsed + ( ticking ? since_last() : 0 ) );
  }

/*--------------------------------------------------------------------------*/

private:

/*--------------------------------------------------------------------------*/

 using clock_type = std::chrono::steady_clock;

 /// time, in seconds, elapsed since the clock was last started

 double since_last( void ) const {
  return( std::chrono::duration< double >( clock_type::now() - last
                                           ).count() );
  }

 bool ticking;                 ///< true if the clock is ticking
 double elapsed;               ///< time accumulated by the past runs
 clock_type::time_point last;  ///< when the current run started

 };  // end( class( DCRtimer ) )

/*--------------------------------------------------------------------------*/
/*--------------------------- CLASS DCR ------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// standard abstract interface for solvers of Delay Constrained Routing
/// problems
/** This class defines a standard abstract interface for solvers of
 *  Delay Constrained Routing (DCR) problems. A DCR problem consists of
 *  an (Integer) Multi-Commodity Flow part plus a max-delay constraint part.
 * The Multi-Commodity Flow part can be single/multi flow and single/multi
 * path.
 * The max-delay part is defined via network calculus and different delay
 * formula
 *  can be used according to different network traffic shapers (e.g.,
 *  Strictly-Rate-Proportional, Weakly-Rate-Proportional, Frame-Based).
 *
 * More precisely, the data of the problem consist of a (directed) network
 * G = ( N , A ), |N| = numNodes and |A| = numLinks, where each node i has a
 * processing delay d_i (DCRNode::delay) and each link (i, j) has a
 * transmission speed s_{ij} (DCRLink::speed), a propagation/queueing delay
 * del_{ij} (DCRLink::delay), a mutual capacity C_{ij} (DCRLink::capacity) and
 * a unit routing cost c_{ij} (DCRLink::cost). A set of numFlows flows has to
 * be routed on G; each flow k is described by a source node o_k
 * (DCRFlow::sourcenode), a sink node t_k (DCRFlow::sinknode), a token-bucket
 * traffic profile with burst b_k (DCRFlow::burst) and rate rho_k
 * (DCRFlow::rate), a deadline D_k (DCRFlow::deadline) and, possibly, a
 * per-arc cost/capacity override (DCRFlow::costs, DCRFlow::caps). Solving
 * the DCR problem means finding, for each flow k, a routing (a single path
 * or, in the multi-path variants, a set of paths) together with the local
 * transmission rate r_{ij}^k > 0 reserved for k on each arc (i, j) it uses,
 * so as to
 * \f[
 *  \min \sum_{ k } \sum_{ (i,j) \in A } c_{ij} \, r_{ij}^k \, x_{ij}^k
 * \f]
 * subject to the customary multi-commodity flow conservation constraints,
 * to the arc capacity constraints
 * \f[
 *  \sum_{ k } r_{ij}^k \, x_{ij}^k \leq C_{ij} \quad (i,j) \in A
 * \f]
 * and, for every flow k, to a worst-case end-to-end delay constraint
 * \f[
 *  \Delta_k( x^k , r^k ) \leq D_k
 * \f]
 * where the exact form of \f$ \Delta_k() \f$ depends on the network
 * calculus formula selected via DCRDelay (dtype); for instance, under a
 * Strictly-Rate-Proportional (SRP) traffic shaper the worst-case delay
 * experienced by flow k when routed on a (single) path P_k is
 * \f[
 *  \Delta_k( x^k , r^k ) = \frac{ b_k }{ \min_{ (i,j) \in P_k } r_{ij}^k } +
 *   \sum_{ (i,j) \in P_k } \left( \frac{ MTU }{ r_{ij}^k } +
 *   \frac{ MTU }{ s_{ij} } + del_{ij} + d_i \right)
 * \f]
 * i.e., the sum of the "burst delay" (governed by the smallest rate
 * reserved along the path) and, for each hop, of the packetization,
 * transmission, propagation and node processing delays. Concrete
 * derived classes (e.g., DCR_SPT for the Single-Flow Single-Path case, or
 * the Lagrangian-relaxation-based solvers built around
 * DCRLagrangianSolver) implement DCRsolve() according to the particular
 * combination of (single/multi flow, single/multi path, delay formula)
 * that they support. */

class DCR
{

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/
public:
/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
 /** @name Public types
  *  @{ */

 /// the possible outcomes of a call to DCRsolve()
 /** DCR solver status*/
 enum DCRStatus {
  OK = 0 ,     ///<  solver found a feasible solution
  Stopped ,    ///<  solver stopped
  Infeasible , ///<  problem infeasible
  Unbounded ,  ///<  problem unbounded
  Error       ///<  solver error
 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the network calculus formula used to compute the worst-case delay
 /** DCR delay formula: selects which traffic-shaping model is used to
  * turn the (burst, rate) profile of a flow and the per-hop network data
  * (speed, delay, node delay) into the worst-case end-to-end delay
  * \f$ \Delta_k() \f$ that has to be compared against the flow deadline
  * (see the class general notes above). */
 enum DCRDelay {
  SRP = 0 , ///< Strictly Rate Proportional
  WRT ,     ///< Weakly Rate Proportional
  FB       ///< Frame Based
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the data describing one flow to be routed
 /** DCR network flow: source and sink node, token-bucket (burst, rate)
  * traffic profile, delay deadline and, optionally, per-arc cost/capacity
  * overrides (costs[] and caps[], each of size numLinks) used in place of
  * the "global" DCRLink::cost / DCRLink::capacity when this flow is
  * considered; a nullptr means that the global values are used instead. */
 struct DCRFlow
 {
  int sourcenode;  ///< flow source
  int sinknode;    ///< flow sink
  double burst;    ///< flow burst
  double rate;     ///< flow rate
  double deadline; ///< flow deadline
  double * costs;  ///< arc-flow costs
  double * caps;   ///< arc-flow individual capacity
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the data describing one (directed) link of the network
 /** DCR network link: start/end node, transmission speed, mutual
  * capacity (shared among all the flows routed on the link), and the
  * propagation/queueing delay and unit routing cost of the link. */
 struct DCRLink
 {
  int startnode;   ///< link start-node
  int endnode;     ///< link end-node
  double speed;    ///< link speed
  double capacity; ///< link mutual capacity
  double delay;    ///< link delay
  double cost;     ///< link cost
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the data describing one node of the network
 /** DCR network node: currently only carries the node processing delay
  * that contributes to the worst-case end-to-end delay of any flow
  * traversing the node. */
 struct DCRNode
 {
  double delay; ///< node delay
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// exception class thrown by (some of) the methods of DCR
 /** Very small class for DCR exceptions
  */
 class DCRException : public std::exception
 {
 public:
  DCRException( const char * const msg = 0 ) { errmsg = msg; }

  const char * what( void ) const throw() { return( errmsg ); }

 private:
  const char * errmsg;
  };

 /** @} ---------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/
 /** @name Constructors
    @{ */

 /// constructor of DCR: initializes all fields to default values
 DCR( void )

 {
  numNodes = numLinks = numFlows = 0;
  MTU = 0;
  dtype = SRP; //i.e., SRP=strictly rate proportional delay

  optEps = fsbEps = 0;

  log = 0;
  verbosity = 0;

  timer = 0;

  tlimit = 0;

  Psol = false;
  UBsol = false;
  }
 /**< Constructor of the base class: gives some default values to the data
   structure of the base class, that can be changed later in the constructors
   of the derived classes.*/

 /** @} ---------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------------- SET METHODS-----------------------------------*/
/*--------------------------------------------------------------------------*/

 /** @name Other initializations
    @{ */

 /// sets the verbosity of the log messages produced by the solver
 virtual void DCRsetVerbosity( int lvl = 0 ) { verbosity = lvl; }

 /**< Sets the level of verbosity (the higher the more verbose).
  * \param lvl int expressing the level of verbosity */

/*--------------------------------------------------------------------------*/

 /// sets the output stream used for the solver's log messages
 //FIXME: would be better to use an *ostream* type to make it more general
 virtual void DCRsetLog( std::ostream * outs ) { log = outs; }

 /**< Creates a log outstream.
  * \param outs ostream pointer*/

/*--------------------------------------------------------------------------*/

 /// creates, resets or destroys the timer used to measure the solution time
 virtual void DCRsetTime( bool timeON = true )
 {
  if( timeON )
   if( timer )
    timer->ReSet();
   else
    timer = new DCRtimer();
  else
   delete timer;
  }

 /**< If timeON is true sets or resets the timer, if false deletes the timer.
  * \param timeON bool value */

/*--------------------------------------------------------------------------*/

 /// starts (or re-starts) the timer, if any was set with DCRsetTime()
 virtual void DCRstartTime( void )
 {
  if( timer == 0 )
   throw(
    DCR::DCRException( "DCR::DCRstartTimer(): failed to start timer" ) );

  timer->Start();
  }

 /**< If timer was set or re-set, starts or re-starts timer ticking. */


/*--------------------------------------------------------------------------*/

 /// stops the timer, if any was set with DCRsetTime()
 virtual void DCRstopTime( void )
 {
  if( timer == 0 )
   throw( DCR::DCRException( "DCR::DCRstopTimer(): failed to stop timer" ) );

  timer->Stop();
  }

 /**< If timer was set, stops timer ticking. */

/*--------------------------------------------------------------------------*/

 /// sets the maximum wall-clock time (in seconds) allowed to the solver
 virtual void DCRsetTimeLimit( long secs ) { tlimit = secs; }

 /**< Sets a time limit for the solver.
  * \param secs long expressing timelimit in seconds */

/*--------------------------------------------------------------------------*/

 /// sets the optimality tolerance used to decide if a solution is optimal
 virtual void DCRsetOptEps( double OE = 0 ) { optEps = OE; }

 /**< In many cases, only an "approximate" solution of the problem is possible;
   alternatively, only an "approximate" solution may be required for the
   purposes of the caller (in order to save time).
   The exact meaning of "approximate" is solver-dependent, but the more
   common ways in which this happens are

   - either the value of the solution is not exactly optimal;

   - or the constraints are not exactly satisfied.

   SetOptEps() tells that any solution that is OE-optimal w.r.t. the value of
   the objective function can be considered optimal.
  * \param OE double expressing optimality tolerance
  */

/*--------------------------------------------------------------------------*/

 /// sets the feasibility tolerance used to decide if a solution is feasible
 virtual void DCRsetFsbEps( double FE = 0 ) { fsbEps = FE; }

 /**< In many cases, only an "approximate" solution of the problem is possible;
   alternatively, only an "approximate" solution may be required for the
   purposes of the caller (in order to save time).
   The exact meaning of "approximate" is solver-dependent, but the more
   common ways in which this happens are

   - either the value of the solution is not exactly optimal;

   - or the constraints are not exactly satisfied.

   SetFsbEps() tells that any solution where the violation of the
   constraints is not larger than FE can be considered feasible.
  * \param FE double expressing feasibility tolerance
  */

 /** @} ---------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR SOLVING THE PROBLEM ---------------------*/
/*--------------------------------------------------------------------------*/
 /** @name Solving the problem
    @{ */


 /// solves the DCR problem loaded with DCRloadProblem()
 virtual DCRStatus DCRsolve( void ) = 0;

 /**< DCR solver, returns a DCRStatus value */

 /** @} ---------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/
 /** @name Reading solver output
  @{ */

 /// returns the objective function value of the current solution
 virtual double DCRgetObj( void ) = 0;

 /**< Returns objective function value in the current solution, if any
  * (see Psol and UBsol).*/

/*--------------------------------------------------------------------------*/

 /// returns the number of paths used by flow k in the (mixed-)integer solution
 virtual int DCRgetUBSolNPaths( int k ) = 0;

 /**< Returns the number of paths of flow k in the current
  * (mixed-)integer solution,
  * if any (see UBsol).
  * \param k flow index
  */

 /// returns path p of flow k in the current (mixed-)integer solution
 virtual int DCRgetUBSolPath( int k , int p , int * X , double * R ) = 0;

 /**< Returns the number of hops in path p of flow k in the current
  * (mixed-)integer solution,
  * if any (see UBsol).
  * X is an array of int of size (n-1) allocated by the user to store
  * the indices of the links in path p. R is an array of double
  * of size (n-1) allocated by the used to store the corresponding rates.
  * \param k flow index
  * \param p path index
  * \param X pointer to an array of int of size (numNodes-1)
  * \param R pointer to an array of double of size (numNodes-1) */

 /// returns the number of paths used by flow k in the current continuous
 /// solution
 virtual int DCRgetPSolNPaths( int k ) = 0;

 /**< Returns the number of (possibly "splitted") "paths" of flow k in
  * the current *continuous* solution,
  * if any (see Psol).
  * Note that in the case of a continuous solution the concept of path no
  * longer exists,
  * because the flow can be splitted. Yet, since we also consider multi-path
  * versions of the problem
  * it makes sense to have different paths also for the continuous case.
  * \param k flow index*/

 /// returns path p of flow k in the current continuous solution
 virtual void DCRgetPSolPath( int k , int p , double * X , double * R ) = 0;

 /**< X is an array of double of size numLinks allocated by the user to store
  * the values of the path (x_ij) variables in the current *continuous*
  * solution,
  * if any (see Psol).
  * R is an array of double of size numLinks allocated by the user to store
  * the values of the rates (r_ij) variables in the current *continuous*
  * solution.
  * \param k flow index
  * \param p path index
  * \param X pointer to an array of double of size numLinks
  * \param R pointer to an array of double of size numLinks */

 /// advances to the next available solution, if any
 virtual void getNewSol()
 {
  Psol = false;
  UBsol = false;
  }
 /**< Moves to the next solution if multiple solutions are given.
  *   If no more solutions are available, sets all solution flags to zero. */

/*--------------------------------------------------------------------------*/

 /// returns the elapsed solution time (see the timer set by DCRsetTime())
 virtual double DCRgetTime( void ) { return( timer ? timer->Read() : 0 ); }

 /**< Returns the elapsed time. If the clock is ticking, returns the *total*
    time since the last Start() without stopping the clock; otherwise,
    returns the total elapsed time of all the past runs of the clock since
    the last reset. */

 /** @} ---------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*----------------------------------LOAD METHOD ----------------------------*/
/*--------------------------------------------------------------------------*/

 /** @name Loading the data of the problem
    @{ */

 /// loads a DCR instance (network + flows) from memory
 virtual void DCRloadProblem( int nnodes , int nlinks , int nflows ,
                              DCRFlow * flows , DCRLink * links ,
                              DCRNode * nodes , double MTU ,
                              DCRDelay deltype ) = 0;

 /**< Reads data from memory and creates a DCR problem.
  * \param nnodes int expressing the number of nodes in the network
  * \param nlinks int expressing the number of links in the network
  * \param nflows int expressing the number of flows in the network
  * \param flows pointer to an array of DCRFlow
  * \param links pointer to an array of DCRLink
  * \param nodes pointer to an array of DCRNode
  * \param MTU double expressing maximum transmit unit
  * \param deltype DCRdelay enum expressing type of delay formula to be used*/

 /** @} ---------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*----------------------------------GET METHODS-----------------------------*/
/*--------------------------------------------------------------------------*/

 /** @name Reading the data of the problem
    @{ */

 /// returns the number of nodes in the network
 virtual int DCRgetNumNodes( void ) const { return( numNodes ); }

 /**< Returns number of nodes in the network. */
/*--------------------------------------------------------------------------*/

 /// returns the number of links in the network
 virtual int DCRgetNumLinks( void ) const { return( numLinks ); }

 /**< Returns number of links in the network. */

/*--------------------------------------------------------------------------*/

 /// returns the number of flows to be routed
 virtual int DCRgetNumFlows( void ) const { return( numFlows ); }

 /**< Returns number of flows in the network. */

/*--------------------------------------------------------------------------*/

 /// returns the Maximum Transmit Unit (MTU) of the network
 virtual double DCRgetMTU( void ) const { return( MTU ); }

 /**< Returns Maximum Transmit Unit (MTU) of the network. */

 /** @} ---------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*----------------------------------MODIFIERS-------------------------------*/
/*--------------------------------------------------------------------------*/

 //TODO: change costs, demand, capacities, open and close arcs.

 /** @name Changing the data of the problem
    @{ */

 /// closes a set of arcs of the network for all the flows
 /** \param whch pointer to an array of na arc indices to be closed
  * \param na number of arcs to be closed */
 virtual void DCRcloseArcs( int * whch , int na ) = 0;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// closes a set of arcs of the network for a specific flow only
 /** \param k index of the flow to which the closing applies
  * \param whch pointer to an array of na arc indices to be closed
  * \param na number of arcs to be closed */
 virtual void DCRcloseArcs( int k , int * whch , int na ) = 0;

 /** @} ---------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

 /** @name Destructor
    @{ */

 /// destructor of DCR: deletes the timer, if any
 virtual ~DCR( void ) { delete timer; }

 /**< Frees up dinamically allocated memory */

 /** @} ---------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 /** @name Protected fields of the class
    @{ */

protected:
 int numNodes; ///< number of nodes in the network
 int numLinks; ///< number of links in the network
 int numFlows; ///< number of flows in the network

 double MTU; ///< maximum transmit unit

 double optEps; ///< optimality tolerance
 double fsbEps; ///< feasibility tolerance

 std::ostream * log;  ///< log output
 int verbosity; ///< verbosity level

 DCRtimer * timer;  ///< timer
 long tlimit;       ///< time limit

 DCRDelay dtype; ///< network delay formula used

 bool Psol;  ///< true if continuous solution was found, false otherwise
 bool UBsol; ///< true if (mixed)-integer solution was found, false otherwise
 };

/** @} ---------------------------------------------------------------------*/

#endif

/*--------------------------------------------------------------------------*/
/*-------------------------- End File DCR.h --------------------------------*/
/*--------------------------------------------------------------------------*/
