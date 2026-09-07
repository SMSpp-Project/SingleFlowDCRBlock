/*--------------------------------------------------------------------------*/
/*--------------------------- File DCRLagrangianSolver.h ----------------------------*/
/*--------------------------------------------------------------------------*/

/** @file
 * Header file for the class DCRLagrangianSolver, which solves the
 * Lagrangian relaxation (w.r.t. the delay constraint) of the Single-Flow
 * Single-Path (SFSP) Delay Constrained Routing problem for a Strictly
 * Rate Proportional (SRP) traffic shaper [see DCR.h], for a *given* value
 * of the minimum per-hop rate r_min, by means of a two-cut line search
 * over the Lagrangian multiplier that repeatedly invokes a Shortest Path
 * Tree solver [see SPT.h].
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
 * \copyright &copy; by Antonio Frangioni
 */ 

/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef DCR_LAGRANGIAN_SOLVER_H
#define DCR_LAGRANGIAN_SOLVER_H

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "DCR.h"
#include "SPT.h"

#include <vector>

using namespace std;

/*--------------------------------------------------------------------------*/
/*--------------------------- CLASS DCRLagrangianSolver --------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// Lagrangian-relaxation solver for the Single-Flow Single-Path SRP DCR problem
/** DCRLagrangianSolver computes the Lagrangian dual, with respect to the
 * delay constraint, of the Single-Flow Single-Path (SFSP) Delay
 * Constrained Routing problem under a Strictly Rate Proportional (SRP)
 * traffic shaper [see DCR.h for the general DCR model], for a single flow
 * with source/sink nodes, burst and deadline given by a DCR::DCRFlow, and
 * for a *given* value r_min of the minimum rate that has to be reserved
 * on every arc of the path (r_min is typically driven from the outside,
 * e.g. by DCR_Benders / BenBound, which explores a discrete set of
 * candidate values). Only arcs whose (possibly temporarily reduced, see
 * closeArcs() / openArcs()) capacity is >= r_min can appear on the path;
 * this "reduced graph" is recomputed by setReducedGraph() whenever
 * r_min or the arc capacities change.
 *
 * Denoting by P a source-sink path in the reduced graph and by r_{ij} the
 * rate reserved on each arc (i, j) of P, the problem solved is
 * \f[
 *  \min \sum_{ (i,j) \in P } cost_{ij} \, r_{ij}
 * \f]
 * \f[
 *  \frac{ burst }{ r_{min} } - deadline +
 *  \sum_{ (i,j) \in P } \left( \frac{ MTU }{ r_{ij} } +
 *  \frac{ MTU }{ speed_{ij} } + delay_{ij} + delay_i \right) \leq 0
 * \f]
 * \f[
 *  r_{min} \leq r_{ij} \leq capacity_{ij} \quad (i,j) \in P
 * \f]
 * i.e., minimize the routing cost of a single source-sink path subject to
 * the SRP worst-case delay bound being respected. Dualizing the delay
 * constraint with a multiplier \f$ \lambda \geq 0 \f$ yields, for a fixed
 * \f$ \lambda \f$, the Lagrangian function
 * \f[
 *  L( \lambda ) = \lambda \left( \frac{ burst }{ r_{min} } - deadline
 *  \right) + \min_{ P , r } \sum_{ (i,j) \in P } \left[ cost_{ij} \, r_{ij}
 *  + \lambda \left( \frac{ MTU }{ r_{ij} } + \frac{ MTU }{ speed_{ij} } +
 *  delay_{ij} + delay_i \right) \right]
 * \f]
 * which, since the inner minimization is separable over the arcs of P,
 * decomposes into: (i) for every arc (i, j) of the reduced graph, a
 * closed-form choice of the locally optimal rate r*_{ij}(\lambda) in
 * [r_min, capacity_{ij}] (see getCost() / setSPTcosts()), giving an
 * arc "Lagrangian cost"; (ii) a Shortest Path Tree computation (via SPT,
 * which can handle the possibly negative arc costs arising for \f$
 * \lambda = 0 \f$ or small \f$ \lambda \f$) from the source to the sink
 * on these costs, whose value, plus the constant term above, is
 * \f$ L( \lambda ) \f$.
 *
 * \f$ L( \lambda ) \f$ is concave in \f$ \lambda \f$: every path/rate
 * combination (\f$ \alpha , \beta \f$) produced by a Shortest Path Tree
 * solve at some trial \f$ \lambda \f$, with \f$ \alpha \f$ the routing
 * cost and \f$ \beta \f$ the (possibly negative) residual delay slack of
 * that path, defines an affine "cut" \f$ \alpha + \lambda \beta \f$ that
 * majorizes \f$ L( \lambda ) \f$ and touches it at the \f$ \lambda \f$
 * where it was generated. The Lagrangian dual
 * \f$ \max_{ \lambda \geq 0 } L( \lambda ) \f$ is then found by a simple
 * *two-cut line search* (Solve()): the class keeps the best cut with
 * non-negative slope (pCut) and the best one with negative slope (mCut);
 * at each iteration the multiplier \f$ \lambda \f$ is moved to the
 * intersection of the two lines, a new cut is generated there by solving
 * the Shortest Path Tree subproblem, and it replaces whichever of pCut /
 * mCut has the same sign of slope; the process stops when the value of
 * the two-cut intersection (InterVal) and the newly computed \f$ L(
 * \lambda ) \f$ (ObjVal) coincide up to the relative tolerance eps. Every
 * time a path/rate combination with \f$ \beta \leq 0 \f$ (i.e., already
 * delay-feasible) is found, its cost is recorded as a primal heuristic
 * value (HeurVal), a valid upper bound on the true SFSP DCR optimum; when
 * \f$ \lambda = 0 \f$ is itself dual-optimal (checked cheaply by is0opt()
 * with at most two Shortest Path Tree solves) there is no duality gap and
 * HeurVal is the true optimum.
 *
 * To support the repeated re-solves with different r_min values that are
 * typical of a Benders/branch-and-bound exploration (see updrmin()),
 * DCRLagrangianSolver keeps a bounded cache (Cuts, up to maxCutSize
 * entries) of the cuts generated by previous solves, together with the
 * associated rate solutions; Reopt() tries to "warm-start" the line
 * search by discarding the cuts that are no longer feasible for the new
 * r_min, correcting the slope of the survivors (whose constant delay term
 * burst / r_min has changed) and resuming the two-cut search on the
 * cached cuts alone, only falling back to a full re-initialization
 * (Inizial()) if this is not possible. */

class DCRLagrangianSolver
{
 
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/	
	public:

/*--------------------------------------------------------------------------*/
/*--------------------------- Inf() and Eps() ------------------------------*/
/*--------------------------------------------------------------------------*/
/// helper to extract the "+ infinity" value for a basic type
/** Very small class to simplify extracting the "+ infinity" value for a
    basic type; just use Inf<type>(). */

 template <typename T>
  class Inf {
   public:
  Inf() {}
  operator T() { return( std::numeric_limits<T>::max() ); }
  };

/*--------------------------------------------------------------------------*/
/// helper to extract the "machine epsilon" for a basic type
/** Very small class to simplify extracting the "machine epsilon" for a
    basic type; just use Eps<type>(). */

 template <typename T>
  class Eps {
   public:
  Eps() {}
  operator T() { return( std::numeric_limits<T>::epsilon() ); }
  };

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public types
 *  @{ */

/// the possible outcomes of a call to Solve()
  enum LAGStatus
   {
    OK,   ///<  solver found a feasible solution
    Infeasible,   ///<  problem infeasible
    Unbounded,    ///<  problem unbounded
    Error
   };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/// a single affine "cut" (supporting line) of the Lagrangian function
/** Represents the affine function \f$ q + \lambda m \f$ generated by a
 * Shortest Path Tree solve at some trial multiplier: m is the slope
 * (the residual delay slack \f$ \beta \f$ of the corresponding path/rate
 * solution) and q is the intercept (its routing cost \f$ \alpha \f$).
 * pCut and mCut, the two currently best cuts with non-negative and
 * negative slope respectively, define the current two-cut approximation
 * of the Lagrangian dual function used by the line search in Solve(). */
  struct LinearCut
  {
     	double m; ///< slope of the cut (the delay slack beta)
	    double q; ///< intercept of the cut (the routing cost alpha)
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/// a saved cut together with the rate solution that generated it
/** Structure used to support the "warm-start" reoptimization implemented
 * in Reopt(): besides the cut itself (c), it stores the rate solution
 * RSol (of size RSolsize) that produced it, and the value of r_min that
 * was in force when the cut was generated (rmin), so that the cut can
 * later be checked for continued feasibility and its slope corrected if
 * r_min has since changed. */
  struct Cut_Val
  {
      LinearCut c; ///< the cut (slope/intercept)
      vector<double> RSol;  ///< the rate solution that generated the cut
      int RSolsize; ///< the size of RSol (number of hops of the path)

      double rmin; ///< the r_min in force when the cut was generated
  };

/** @} ---------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor
 *  @{ */

/// constructor of DCRLagrangianSolver: initializes all fields to default values
   DCRLagrangianSolver(void);

/** @} ---------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*---------------------------------- METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Solving the problem
 *  @{ */

/// solves the Lagrangian dual by the two-cut line search described above
/** Runs the two-cut line search over the Lagrangian multiplier lambda
 * (alternating Shortest Path Tree solves with the update of the pCut /
 * mCut supporting lines) until the intersection value and the newly
 * computed Lagrangian function value agree within the relative tolerance
 * eps [see the class general notes]. Requires LoadProblem() to have been
 * called beforehand; the outcome can be read with getStatus(), and the
 * various components of the optimal cut/path/rate solution with the
 * get...() methods below. If reoptflag is set (see updrmin()), Reopt()
 * is first attempted to warm-start the search from the cached cuts. */
   void Solve();

/** @} ---------------------------------------------------------------------*/
/*-------------------- Methods for reading the solution ---------------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the solution
 *  @{ */

/// returns the status of the last call to Solve() / isFeasible()
   LAGStatus getStatus();

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/// returns the optimal Lagrangian multiplier (the pCut/mCut intersection)
   double getLambda();

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/// returns the value of the two-cut intersection at the optimal lambda
/** Returns the value of the two-cut approximation (InterVal) at the
 * optimal lambda found by the line search; at convergence this
 * coincides, up to the tolerance eps, with getObjVal(). */
   double getInterVal();

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/// returns the Lagrangian function value at the optimal lambda
   double getObjVal();

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/// returns the optimal path (as arc indices in the reduced graph)
/** Returns the optimal path found by the last Shortest Path Tree solve of
 * Inizial()/Solve(), as a vector of arc indices in the reduced graph [see
 * getRedGraphPos() to map them back to the original graph]. This is
 * meaningful only when the optimal lambda is 0, i.e. when the delay
 * constraint is not active at the optimum. */
   vector<int> getPath();

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/// returns the number of arcs (hops) of the optimal path
   int getnumHops();

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/// returns the best delay-feasible primal solution value found so far
/** Returns the cost of the best path/rate combination encountered during
 * the line search that was already delay-feasible (\f$ \beta \leq 0 \f$),
 * i.e., a valid upper bound (primal heuristic value) on the true SFSP DCR
 * optimum [see the class general notes]. */
   double getHeurVal();

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/// returns the slope (delay slack beta) of the cut at the optimum
   double getOptBeta();

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/// returns the vector of the r_ij rates of the optimal solution
   vector<double> getRSol();

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/// returns the arc costs f_ij corresponding to the r_ij of getRSol()
   vector<double> getRSolCosts();

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/// returns the positions of the reduced graph arcs in the original graph
/** Returns, for each arc of the current reduced graph [see
 * setReducedGraph()], its index in the original (full) graph; used e.g.
 * to check the Bellman optimality conditions on the full graph. */
   vector<int> getRedGraphPos();

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/// returns the number of arcs of the reduced graph
   int getCardRedGraph();

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/// returns the positions, in the original graph, of the arcs *not* in the reduced graph
/** Returns the complement of getRedGraphPos(): the indices, in the
 * original graph, of the arcs that have been excluded from the reduced
 * graph because their capacity is below r_min; used for checks related
 * to the non-convexity of the per-arc Lagrangian cost. */
   vector<int> getComplGraph();

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/// returns a convex combination of the positive- and negative-slope solutions
/** Returns a rate vector obtained as a convex combination of the rate
 * solutions associated with the currently best positive-slope and
 * negative-slope cuts (solpos / solneg), weighted so as to (approximately)
 * satisfy the delay constraint; used as an auxiliary check/heuristic
 * solution (e.g. for comparison against a MILP solver such as CPLEX). */
   vector<double> getCheckSol();

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/// returns the positions (in the graph) of the positive-slope check solution
   vector<int> getCheckSolPos();

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/// returns the positions (in the graph) of the negative-slope check solution
   vector<int> getCheckSolNeg();

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/// returns the Shortest Path Tree node potentials of the optimal solution
   vector<double> getSPLabels();

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/// returns the number of line-search iterations performed by Solve()
   int getNumIte();

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
   //vector<Cut_Val> getCuts();
   /*<returns the set of cuts currently under consideration (is this needed?)*/

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/// returns the elapsed solution time (see the timer set by DCRsetTime())
   double DCRgetTime( void );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/// returns the vector of the r_ij rates of the optimal solution, indexed
/// over *all* the arcs of the (unreduced) network (0 if the arc is unused)
   vector<double> getRSOLS() { return(RSOLS); }

/** @} ---------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------- Methods for checking feasibility ---------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for checking feasibility
 *  @{ */

/// checks whether the Lagrangian dual has a feasible solution, without solving it
/** Runs only the initialization phase (Inizial()) without performing the
 * full line search, and reports whether the Lagrangian relaxation turned
 * out to be infeasible (i.e., whether even the best achievable routing
 * cost exceeds the trivial upper bound LimitVal). Cheaper than a full
 * Solve() when only feasibility is of interest; since Inizial() is run,
 * a subsequent Solve() call will skip it (see inizialflag).
 * \return true if the problem is infeasible, false if a feasible solution
 *         (possibly requiring further line-search iterations) exists */
   bool isFeasible ();

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/// checks whether lambda = 0 is already the optimal Lagrangian multiplier
/** Solves the Shortest Path Tree subproblem at lambda = 0 and, if that
 * path already violates the delay constraint, also at a small lambda =
 * 1e-6, comparing the resulting Lagrangian function values to decide
 * whether increasing lambda beyond 0 would improve the dual value. This
 * costs at most two Shortest Path Tree solves and is meant to be called
 * right after LoadProblem(), typically to quickly detect (e.g. from
 * BenBound, during Benders reoptimization) that no relaxation gap can
 * arise for the current data.
 * \return true if lambda = 0 is NOT optimal (i.e., a positive lambda
 *         should be searched for), false if lambda = 0 is already optimal */
   bool is0opt();

/** @} ---------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*----------------------------------LOAD METHODS ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Loading the data of the problem
 *  @{ */

/// loads a full DCR instance (network, single flow and r_min) from memory
/** Copies the network and flow data into the internal data structures and
 * forwards them to the internal SPT object, resetting all solution flags.
 * Must be called before Solve() and before any of the get...() methods,
 * otherwise they will either throw or return default/stale values.
 * \param nnodes number of nodes of the network
 * \param nlinks number of links of the network
 * \param flow the (single) DCR::DCRFlow to be routed
 * \param links pointer to an array of nlinks DCR::DCRLink
 * \param nodes pointer to an array of nnodes DCR::DCRNode
 * \param mtu the Maximum Transmit Unit of the network
 * \param r_min the minimum rate to be reserved on every arc of the path */
   //FIXME: should this also require the desired precision? (implemented with machine precision)
   void LoadProblem (int nnodes, int nlinks, DCR::DCRFlow flow, DCR::DCRLink *links, DCR::DCRNode *nodes, double mtu, double r_min);

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/// loads a new flow, keeping the previously loaded network data
/** Alternative, lighter-weight loading method to be used only *after* the
 * network has already been loaded with LoadProblem(memory): it changes
 * the flow being routed (and the corresponding arc costs) while reusing
 * the existing network data structures, and reopens all the arcs that had
 * been closed with closeArcs() [see the FIXME in the .cpp for the
 * possibility of not doing so].
 * \param flow the new DCR::DCRFlow to be routed */
   void LoadProblem(DCR::DCRFlow flow);

/** @} ---------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*----------------------------------MODIFIERS-------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Changing the data of the problem
 *  @{ */

/// changes r_min without reloading the rest of the problem data
/** Sets a new value for r_min; if it actually differs from the current
 * one (or the problem had not been solved yet), all solution flags are
 * reset and reoptflag is set so that the next Solve() call will attempt
 * to warm-start via Reopt() instead of a full re-initialization.
 * \param r_min the new minimum rate to be reserved on every arc */
   void updrmin(double r_min);

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/// closes a set of arcs, i.e., sets their (working) capacity to 0
/** Closes the arcs whose indices are passed in arcs by setting their
 * entry in ModCaps (the working copy of the arc capacities) to 0; the
 * original capacity stored in Links is left untouched, so that
 * openArcs() can later restore it.
 * \param arcs vector of the indices of the na arcs to be closed
 * \param na number of arcs to be closed */
   void closeArcs(vector<int> arcs, int na);

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/// reopens a set of previously closed arcs
/** Reopens the arcs whose indices are passed in arcs by restoring their
 * entry in ModCaps (the working copy of the arc capacities) to the
 * original capacity stored in Links.
 * \param arcs pointer to an array of the indices of the na arcs to reopen
 * \param na number of arcs to be reopened */
   void openArcs(int * arcs, int na);

/** @} ---------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- TIMING METHODS --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Handling the timer
 *  @{ */

/// sets a time limit (in seconds) for the solver
   void DCRsetTimeLimit(long secs);

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/// creates, resets or destroys the timer used to measure the solution time
   void DCRsetTime(bool timeON);

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/// starts (or re-starts) the timer, if any was set with DCRsetTime()
   void DCRstartTime();

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/// stops the timer, if any was set with DCRsetTime()
   void DCRstopTime();

/** @} ---------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Destructor
 *  @{ */

/// destructor of DCRLagrangianSolver: releases all dynamically allocated memory
	~DCRLagrangianSolver();

/** @} ---------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Private methods of the class
 *  @{ */

/// deep-copies the network and flow data into the internal structures
/** Allocates Links and Nodes and copies into them the data pointed to by
 * links/nodes, and stores flow into Flow; also initializes ModCaps (the
 * working copy of the arc capacities used by closeArcs()/openArcs()) to
 * the original capacities. Called by LoadProblem(memory) after
 * clean_up(). */
    void copyDataArray(DCR::DCRFlow flow, DCR::DCRLink *links, DCR::DCRNode *nodes);

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/// computes a trivial worst-case upper bound (LimitVal) on the routing cost
/** Sums, over all the arcs of the current reduced graph with positive
 * cost, the product cost_ij * capacity_ij, i.e., the cost of saturating
 * every positive-cost arc to its capacity: a (very loose but valid)
 * upper bound on the cost of any feasible routing, used by Inizial() to
 * detect that the Lagrangian relaxation (and hence the original problem)
 * is infeasible. */
    void getLimitVal();

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/// returns the Lagrangian cost of reserving rate r on arc lindex
/** Computes \f$ \lambda \, \bar{l}_{ij} + cost_{ij} \, r + \lambda \,
 * MTU / r \f$, where \f$ \bar{l}_{ij} = MTU / speed_{ij} + delay_{ij} +
 * delay_i \f$ is the rate-independent part of the per-hop delay of arc
 * lindex of the reduced graph (RedGraLinks) at the current value of
 * lambda [see the class general notes].
 * \param r the rate reserved on the arc
 * \param lindex index of the arc in the reduced graph
 * \return the corresponding Lagrangian arc cost, to be fed to the
 *         Shortest Path Tree solver */
    double getCost(double r, int lindex);

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/// initializes the line search by finding the first negative-slope cut
/** Builds the reduced graph and solves the Shortest Path Tree subproblem
 * at increasing values of lambda (starting at lambda = 1, doubling at
 * each step) until either a path with non-positive delay slack (beta) is
 * found, in which case lagstat is set to OK and the search can start
 * from lambda = 0, or the running Lagrangian value exceeds LimitVal, in
 * which case lagstat is set to Infeasible, or the Shortest Path Tree
 * itself fails (negative-cost cycle / disconnected sink), in which case
 * lagstat is also set to Infeasible. In every case both pCut and mCut are
 * updated with the cuts found along the way. */
    void Inizial();

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/// updates the best positive- or negative-slope cut with a new one
/** Given the (alpha, beta) = (intercept, slope) pair of a newly computed
 * cut, replaces pCut (if beta >= 0) or mCut (if beta < 0) with it: the
 * two-cut approximation of the Lagrangian function used by the line
 * search in Solve() and Reopt() is thereby tightened.
 * \param alpha the intercept (routing cost) of the new cut
 * \param beta the slope (delay slack) of the new cut */
    void UpdCut(double alpha, double beta);

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/// (re)builds the reduced graph, keeping only the arcs with capacity >= r_min
/** Scans ModCaps (the working arc capacities, reflecting closeArcs() /
 * openArcs()) and copies into RedGraLinks/Linksp the arcs whose capacity
 * is not smaller than r_min, i.e. those that can legally carry the flow
 * at the required minimum rate; RedGraPos records, for each arc of the
 * reduced graph, its index in the original graph. Must be called again
 * whenever r_min or the arc capacities change. */
    void setReducedGraph();

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/// computes, for the current lambda, the locally optimal rate and cost of every reduced-graph arc
/** For each arc of the reduced graph, chooses the rate r* in [r_min,
 * capacity] minimizing the per-arc Lagrangian cost at the current lambda
 * [see getCost()]: if the arc cost is negative the minimum is attained at
 * r* = capacity, otherwise at r* = clip( sqrt( lambda * MTU / cost ) ,
 * r_min , capacity ); the resulting (rstar, cost) pairs are stored in
 * Linksp, ready to be fed to SPT::updCosts() / SPT::Solve(). */
    void setSPTcosts();

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/// attempts to warm-start the line search from the cached cuts after r_min changed
/** Discards, from the cut cache Cuts, every cut whose associated rate
 * solution is no longer feasible for the new r_min, corrects the slope of
 * the surviving cuts (whose constant delay term burst / r_min has
 * changed) and, if at least one cut of each sign survives, resumes the
 * two-cut line search restricted to the cached cuts (verifying, at each
 * iteration, that no other cached cut gives a lower value at the current
 * lambda). Falls back to signalling that a full Inizial() is needed
 * whenever fewer than two cuts survive, no cut of one of the two signs
 * survives, or the restricted search does not converge within a fixed
 * number of iterations.
 * \return 1 if the caller must fall back to a full re-initialization
 *         (Inizial()), 0 if the warm-start succeeded and pCut/mCut/lambda
 *         are ready for Solve() to resume the line search */
    int Reopt();

/** @} ---------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/	
	
	protected:

/** @name Protected fields of the class
 *  @{ */

	LinearCut pCut; ///< best cut with non-negative slope defining the optimal solution
	LinearCut mCut; ///< best cut with negative slope defining the optimal solution
  vector<Cut_Val> Cuts; ///< cache of the cuts currently in play (for Reopt())
  double InterVal; ///< value of the pCut/mCut intersection
  double lambda;    ///< current/optimal Lagrangian multiplier (intersection point)
  double ObjVal;   ///< value of the Lagrangian function at lambda
  double LimitVal; ///< trivial worst-case upper bound on the routing cost
  double HeurVal; ///< best delay-feasible primal solution value found so far
  double optBeta; ///< slope of the Lagrangian function at the optimum

  vector<double> RSOLS; ///< optimal r_ij, indexed over all arcs (0 if unused)

  int num_ite; ///< number of line-search iterations performed by Solve()

  // data describing the input instance, kept around to be passed to SPT
  double rmin; ///< the minimum rate to be reserved on every arc of the path
  double MTU; ///< the Maximum Transmit Unit of the network
  DCR::DCRNode * Nodes; ///< array of numNodes DCR::DCRNode (network node data)
  DCR::DCRLink * Links; ///< array of numLinks DCR::DCRLink (network link data)
  DCR::DCRFlow  Flow; ///< the (single) flow being routed
  int numNodes,numLinks; ///< number of nodes / links of the network
  int cardRedGraph; ///< number of arcs of the reduced graph
  vector<double> ModCaps; ///< working (possibly closed/reopened) arc capacities
  //int * closedArcs; //vector of the arcs that have been closed
  //int nclar; //commented out because it might be useful to reopen them, but it isn't really necessary

  // data produced by the Shortest Path Tree solves
  int nhops; ///< number of arcs (hops) of the optimal path
  vector<DCR::DCRLink> RedGraLinks; ///< the DCR arcs of the reduced graph
  vector<int> RedGraPos; ///< positions of the reduced graph arcs in the full graph
  vector<int> ComplGraph; ///< positions, in the full graph, of the arcs NOT in the reduced graph
  vector<int> OptPath; ///< positions (in the reduced graph) of the optimal path
  int solvedflag; ///< true if the current data have already been solved
  int inizialflag; ///< true if Inizial() has already been run on the current data
  int reoptflag; ///< true if Solve() should attempt Reopt() (r_min just changed)
  //double * XSol;
  int numHops; ///< number of arcs (hops) of the optimal path (see nhops)
  vector<double> RSol; ///< the r_ij of the optimal solution
  vector<double> RSolCosts; ///< the arc costs f_ij associated with the r*_ij of RSol
  vector<double> SPLabels; ///< Shortest Path Tree node potentials of the optimal solution

  int maxCutSize; ///< maximum number of cuts kept in the Cuts cache for Reopt()

  // auxiliary data used to build a CPLEX-comparable heuristic check solution
  vector<double> solpos; ///< rate solution of the best positive-slope cut
  vector<double> solneg; ///< rate solution of the best negative-slope cut
  double betaneg; ///< slope of the best negative-slope cut
  double betapos; ///< slope of the best positive-slope cut
  int nonposflag; ///< true if no positive-slope cut is available
  vector<int> checkSolPos; ///< positions (in the graph) of solpos
  vector<int> checkSolNeg; ///< positions (in the graph) of solneg

  double eps; ///< relative precision used as convergence tolerance (machine epsilon)
  SPT spt;     ///< the Shortest Path Tree solver used by the line search
  vector<SPT::SPTLink> Linksp; ///< simplified arcs (cost only) fed to spt
  SPT::Status spstat; ///< outcome of the last spt.Solve() call
  LAGStatus lagstat;  ///< outcome of the Lagrangian relaxation (see getStatus())

  OPTtimers *timer;  ///< timer
  long tlimit; ///< time limit

/// releases all dynamically allocated memory and resets the solution data
/** Deletes Links and Nodes and clears all the vectors holding solution and
 * cut-cache data, so that the object is left in the same state as right
 * after construction (net of the timer). Called by LoadProblem(memory)
 * before loading a new instance, and by the destructor. */
  void clean_up();

/** @} ---------------------------------------------------------------------*/

};

#endif

/*--------------------------------------------------------------------------*/
/*-------------------------- End File DCRLagrangianSolver.h ---------------------------*/
/*--------------------------------------------------------------------------*/