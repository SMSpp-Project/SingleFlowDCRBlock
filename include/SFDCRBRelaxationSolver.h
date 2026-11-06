/*--------------------------------------------------------------------------*/
/*-------------------- File SFDCRBRelaxationSolver.h -----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the *concrete* class SFDCRBRelaxationSolver, which
 * implements the RelaxationSolver concept [see ChangeSolver.h] for a
 * SingleFlowDCRBlock, by means of the "Benders with nested Lagrange" scheme
 * of SingleFlowDCRBendersSolver, so that it can solve the nodes of an
 * enumeration tree (as in BranchAndXSolver) branching on the arcs.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Filippo Magi \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni, Filippo Magi
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __SFDCRBRelaxationSolver
 #define __SFDCRBRelaxationSolver
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "SingleFlowDCRBendersSolver.h"

#include "ChangeSolver.h"

#include <atomic>

/*--------------------------------------------------------------------------*/
/*-------------------------- NAMESPACE & USING -----------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*------------------------------- CLASSES ----------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*-------------------- CLASS SFDCRBRelaxationSolver ------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// RelaxationSolver of a SingleFlowDCRBlock that branches on the arcs
/** The relaxation of a node is solved by SingleFlowDCRBendersSolver, which
 * gives a lower bound on the cost of the flow and, when it finds one, a
 * routing that meets the deadline, i.e., a feasible solution of the true
 * problem [see has_true_var_solution()].
 *
 * On acyclic graphs a distinct master pair closes one outgoing arc per
 * child [see branch()]. Otherwise a binary x_a=0 / x_a=1 split covers the
 * parent, including circulations. The closing and requirement Changes (SingleFlowDCRBlockChange) are applied to the Solver and not to
 * the Block [see apply()], which then keeps the arcs out of every routing
 * at the next compute() [see arc_closed()]; the arcs closed or deleted in
 * the Block are kept out as well. A node in which the sink cannot be
 * reached from the source on the arcs left open is infeasible, and
 * compute() says so without running the Benders scheme at all.
 *
 * Children first screen against a finite incumbent using interval bounds.
 * Warm multiplier hints are suspended during compute(), with the configured
 * option restored on exit.
 *
 * An open gap first uses the defining paths at BenBound's final master
 * candidate. If no distinct pair is available, refine_lb() bisects the
 * weakest interval and checks its lower endpoint with the ORIGINAL flow.
 * It stops as soon as that evaluation provides two distinct valid paths,
 * or a bound closes the node. Equal paths alone do not certify global
 * optimality. If neither closure nor a pair is obtained within the budget,
 * the complete binary fallback takes over. With intRoutingFallback=0,
 * compute() returns kStopIter instead. Fully fixed routings use a convex
 * rate solver, returning kLowPrecision if its primal/dual gap stays open.
 *
 * Node profiles report refine_enabled=1. [refine_profile] records
 * reason=branch_pair_found for early exits, or the bound/budget stop reason;
 * original-flow probes are timed separately from interval solves.
 * Records are identified by solver address and a local compute sequence
 * (run); timers are disabled when no log stream is set. */

class SFDCRBRelaxationSolver : public SingleFlowDCRBendersSolver ,
                               public RelaxationSolver
{

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

public:

/*--------------------------------------------------------------------------*/
/*--------------------- PUBLIC METHODS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 *  @{ */

/*--------------------------------------------------------------------------*/
 /// constructor: no arc is closed

 SFDCRBRelaxationSolver( void ) : SingleFlowDCRBendersSolver() {
  collect_branch_paths( true );
  }

/*--------------------------------------------------------------------------*/
 /// destructor: does nothing

 ~SFDCRBRelaxationSolver() override = default;

/** @} ---------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */

 /// set the (pointer to the) Block that the Solver has to solve
 /** As SingleFlowDCRBendersSolver::set_Block(); besides, all the arcs are
  * open for the Solver [see apply()]. */

 void set_Block( Block * block ) override;

 /*--------------------------------------------------------------------------*/

 using SingleFlowDCRBendersSolver::set_par;  // keep the other set_par()

 /// Refinement budgets are per tranche; zero intRefineSplits disables
 /// splitting. intRefineStall=16 stops after consecutive nonimproving splits
 /// (zero disables this test). intRoutingFallback=1 enables binary branching.
 /// dblRefineTime=1 bounds a tranche between oracle calls; zero is unlimited.
 /// The heap survives unchanged-node computes and is invalidated by model
 /// notifications and applied/undone node Changes.
 enum { intRefineSplits = SingleFlowDCRBendersSolver::intMaxCuts + 1,
        intRefineStall, intRoutingFallback, intLastRABaX };
 enum { dblRefineTime = Solver::dblLastAlgPar, dblLastRABaX };
 void set_par( idx_type par , int value ) override;
 idx_type get_num_int_par() const override { return intLastRABaX; }
 int get_int_par( idx_type par ) const override;
 int get_dflt_int_par( idx_type par ) const override;
 idx_type int_par_str2idx( const std::string & name ) const override;
 const std::string & int_par_idx2str( idx_type par ) const override;
 idx_type get_num_dbl_par() const override { return dblLastRABaX; }
 double get_dflt_dbl_par( idx_type par ) const override {
  return par == dblRefineTime ? 1.0 :
   SingleFlowDCRBendersSolver::get_dflt_dbl_par( par );
 }
 idx_type dbl_par_str2idx( const std::string & name ) const override {
  return name == "dblRefineTime" ? dblRefineTime :
   SingleFlowDCRBendersSolver::dbl_par_str2idx( name );
 }
 const std::string & dbl_par_idx2str( idx_type par ) const override {
  static const std::string name = "dblRefineTime";
  return par == dblRefineTime ? name :
   SingleFlowDCRBendersSolver::dbl_par_idx2str( par );
 }

 /// set the double parameters of the Solver
 /** dblRelAcc [1e-6] controls closure of the relative gap in compute().
  * dblUpCutOff is the global incumbent supplied by the driver: a lower
  * bound reaching it
  * stops computation with kCutOff, without preparing a branching pair.
  * dblAbsAcc is stored; the other parameters are passed to the base. */

 void set_par( idx_type par , double value ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] double get_dbl_par( idx_type par ) const override;

 /// Invalidate verified solution data when model data change.
 void add_Modification( sp_Mod & mod ) override;

 /// Cost of the current verified candidate, reusing the per-node selection.
 OFValue get_var_value() override;


/** @} ---------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE MODEL ----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Solving the relaxation of the current node
 *  @{ */

 /// solve the relaxation of the current node
 /** Loads the instance with the arcs closed by the Solver, then runs the
  * "Benders with nested Lagrange" scheme of SingleFlowDCRBendersSolver [see
  * SingleFlowDCRBendersSolver::compute()] on it. Returns kInfeasible if the
  * sink cannot be reached from the source on the open arcs, or if the
  * scheme finds the node infeasible. Children with a finite incumbent
  * first try an independent bound and may skip BenBound. Returns kError
  * if the Block has no source or
  * no sink or if the scheme ends in error, kCutOff if the lower bound
  * reaches a finite dblUpCutOff, and kOK otherwise. */

 int compute( bool changedvars = true ) override;

/** @} ---------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Accessing the found solutions (if any)
 *  @{ */

 /// return a valid lower bound on the cost of the flow in the node
 /** The lower bound of SingleFlowDCRBendersSolver, or the independent
  * certificate if it is larger. The certificate is also available when
  * early screening skips BenBound. It is obtained by bisecting the domain
  * of r_min and taking the minimum of
  * the lower bounds on the resulting intervals. Numerical gap closure
  * never replaces the calculated LB with the candidate's UB.
  * +Inf if the node is infeasible, -Inf if no bound is available. */

 OFValue get_lb( void ) override;

/*--------------------------------------------------------------------------*/
 /// return a valid lower bound on the cost of the flow (as get_lb())

 OFValue get_true_lb( void ) override { return( get_lb() ); }

/*--------------------------------------------------------------------------*/
 /// return a valid upper bound on the cost of the flow (as get_ub())
 /** The cost of the feasible routing found [see
  * SingleFlowDCRBendersSolver::get_ub()], +Inf if there is none. */

 OFValue get_true_ub( void ) override { return( get_ub() ); }

/*--------------------------------------------------------------------------*/
 /// tells whether a feasible routing of the true problem is available
 /** The routing found by the relaxation is a solution of the true problem
  * whenever it is there at all [see
  * SingleFlowDCRBendersSolver::has_var_solution()], since the arcs closed
  * by the Solver are open in the Block. */

 bool has_true_var_solution( void ) override {
  return( has_var_solution() );
  }

/*--------------------------------------------------------------------------*/
 /// write the feasible routing in the Variable of the Block

 void get_true_var_solution( Configuration * solc = nullptr ) override {
  get_var_solution( solc );
  }

/*--------------------------------------------------------------------------*/
 /// physically construct the Solution holding the feasible routing

 Solution * get_true_solution( Configuration * solc = nullptr ) override {
  return( get_Solution( solc ) );
  }

/** @} ---------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/
/** @name Changing the data of the model
 *  @{ */

 /// branch on the first arc where two paths of the relaxation split
 /** Uses the defining paths at BenBound's final master candidate. Equal
  * paths are retained as a result, and compute() closes the node only
  * when its global bound meets the verified primal cost. If no distinct
  * pair is available and the gap remains open, a complete binary split is
  * used. The outgoing-pair split is restricted to acyclic graphs, where
  * binary unit flows cannot contain circulations. */

 std::vector< Change * > branch( void ) override;

/*--------------------------------------------------------------------------*/
 /// apply the Change; SingleFlowDCRBlockChange are applied to the Solver
 /** A SingleFlowDCRBlockChange, as those produced by branch(), closes or
  * opens the arc for the Solver alone, the Block being left as it is, and
  * the next compute() solves the relaxation with it; with doUndo, the
  * returned Change brings the arc back to the state it had for the Solver.
  * Any other Change is applied to the Block. */

 Change * apply( Change * chg , bool doUndo = false ) override;

/** @} ---------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

protected:

/*--------------------------------------------------------------------------*/
/*--------------------------- PROTECTED METHODS ----------------------------*/
/*--------------------------------------------------------------------------*/

 /// true if arc i is closed by the Solver, or closed or deleted in the Block
 /** Tells SingleFlowDCRBendersSolver::load_BenBound() which arcs to keep
  * out of every routing [see SingleFlowDCRBendersSolver::arc_closed()]. */

 bool arc_closed( Index i ) override;

 /// A reload invalidates the verified candidate even outside compute().
 void load_BenBound() override;

 /// Cache the base class's feasibility checks and selected candidate cost.
 SolSource pick_solution_source( double feps = 1e-5 ) override;


/*--------------------------------------------------------------------------*/
/*---------------------------- PROTECTED FIELDS ----------------------------*/
/*--------------------------------------------------------------------------*/

 std::vector< bool > v_closed;  ///< the arcs closed by the Solver

 bool f_infeasible = false;     ///< the node has no path on the open arcs

 /// the certificate of the node [see get_lb()], -Inf if there is none
 OFValue f_cert = - SMSpp_di_unipi_it::Inf< OFValue >();

 double f_RelAcc = 1e-6;        ///< the parameter dblRelAcc

 double f_AbsAcc = 1e-6;        ///< the parameter dblAbsAcc

 OFValue f_UpCutOff = SMSpp_di_unipi_it::Inf< OFValue >(); ///< dblUpCutOff

 /// the middle point of the r_min interval that gives the lower bound
 /** Set by the fallback refine_lb(), 0 if it has not been run;
  * retained for refinement diagnostics. */

 double f_zstar = 0;

 double f_za = 0; ///< the left end of the r_min interval that gives the bound

 int f_refine_splits = 200; ///< maximum bisections per tranche
 int f_refine_stall = 16;   ///< consecutive nonimproving splits before branching
 int f_routing_fallback = 1; ///< complete binary branching when no path pair
 double f_refine_time = 1; ///< seconds per tranche, zero means no time limit
 std::vector< bool > v_required; ///< x_a=1 is a real restriction
 bool f_terminal = false;
 SingleFlowDCRBlock::Vec_double f_terminal_x, f_terminal_r;
 OFValue f_terminal_cost = SMSpp_di_unipi_it::Inf< OFValue >();
 void get_solution_vectors( SingleFlowDCRBlock::Vec_double & X,
                            SingleFlowDCRBlock::Vec_double & R,
                            bool heuristic = false ) override;
 bool solution_matches_node( bool heuristic );
 bool acyclic();
 bool routing_fixed();
 int solve_fixed_routing();
 void polish_solution();

 /// distinct paths cached by compute() for branch(), empty if none found
 std::vector< int > v_branch_plus , v_branch_minus;

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

private:

 struct BoundInterval { double a , b , lb; };
 /// Current node only: retain the interval heap between screening and fallback.
 std::vector< BoundInterval > v_bound_intervals;
 std::size_t f_partition_revision = std::numeric_limits<std::size_t>::max();
 bool f_bound_only = true; ///< BenBound has not been solved for this node

 /// A notification can arrive from another thread; only compute/query
 /// methods access the actual cache. Model notifications advance its key.
 std::atomic< bool > f_bound_dirty{false};
 std::atomic< std::size_t > f_solution_revision{ 0 };
 std::size_t f_cached_revision = std::numeric_limits< std::size_t >::max();
 SolSource f_solution_source = SolSource::None;
 double f_solution_feps = 0;
 OFValue f_solution_cost = 0;
 void invalidate_solution_cache() {
  f_solution_revision.fetch_add( 1 , std::memory_order_relaxed );
  }

 /// Per-compute diagnostics; wall-clock timers are enabled only with f_log.
 struct Profile {
  double probe_seconds = 0 , refine_seconds = 0;
  std::size_t probes = 0 , refine_calls = 0 , interval_solves = 0;
  std::size_t refine_splits = 0;
  std::size_t cache_hits = 0;
  std::size_t solution_checks = 0 , solution_cache_hits = 0;
 } f_profile;
 std::size_t f_compute_id = 0; ///< local sequence; identify with solver pointer

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE METHODS ------------------------------*/
/*--------------------------------------------------------------------------*/

 /// use a distinct master pair or refine until the first pair is found
 /** Called by compute() before its final gap check. Refinement also stops
  * when a global bound closes the node, or its budget is exhausted. */

 void prepare_branch();

/*--------------------------------------------------------------------------*/

 /// lower bound on the cost of the flows whose r_min is in [ a , b ]
 /** Solves the Lagrangian relaxation in which every rate is at least a and
  * the burst term is that of b, which contains all the routings with r_min
  * in [ a , b ]; +Inf if it is empty, -Inf if nothing can be said. */

 double interval_lb( double a , double b );

/*--------------------------------------------------------------------------*/
 /// refine the lower bound of the node by bisecting the domain of r_min
 /** Bisects [ max{ rho , burst / deadline } , largest capacity ] keeping the
  * interval with the smallest lower bound, until it reaches dblUpCutOff
  * or the local upper bound within tolerance, or a distinct original-flow
  * branching pair is found. Also stops on a narrow interval or intRefineSplits
  * bisections. With cutoff_only, uses the incumbent alone and performs no
  * branching probes: BenBound has not been solved on this node yet. Uses
  * 16 bisections before BenBound, extended by at most eight when the bound
  * is improving and close to the incumbent. max_splits limits an explicit
  * pre-branch attempt; -1 uses the default budget. A positive max_seconds
  * also limits elapsed time between interval solves. The heap is retained for
  * the later fallback, avoiding repeated interval solves. The least
  * lower bound found goes in f_cert, and the middle point of the interval
  * that attains it in f_zstar. */

 void refine_lb( bool cutoff_only = false , int max_splits = -1 ,
                 double max_seconds = -1 );

 /*--------------------------------------------------------------------------*/

 /// true if the sink can be reached from the source on the open arcs

 bool sink_reachable( void );

/*--------------------------------------------------------------------------*/
 /// compute the certificate of the node [see get_lb()] in f_cert

 //void certify( void );

/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;
 // insert SFDCRBRelaxationSolver in the factory

/*--------------------------------------------------------------------------*/

 }; // end( class( SFDCRBRelaxationSolver ) )

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* SFDCRBRelaxationSolver.h included */

/*--------------------------------------------------------------------------*/
/*------------------- End File SFDCRBRelaxationSolver.h --------------------*/
/*--------------------------------------------------------------------------*/
