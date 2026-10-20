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
 * The children of a node close one arc each [see branch()]. The closing
 * Changes (SingleFlowDCRBlockChange) are applied to the Solver and not to
 * the Block [see apply()], which then keeps the arcs out of every routing
 * at the next compute() [see arc_closed()]; the arcs closed or deleted in
 * the Block are kept out as well. A node in which the sink cannot be
 * reached from the source on the arcs left open is infeasible, and
 * compute() says so without running the Benders scheme at all.
 *
 * When a gap remains, compute() reuses the intervals computed by BenBound
 * to locate two distinct paths for branch(). Only if these probes fail
 * does it refine an independent lower bound [see get_lb()]. This happens
 * before returning to the driver, since refinement can close the node.
 *
 * When a log stream is set [see Solver::set_log()], compute(), branch()
 * write diagnostic information on it, including the number of probes and
 * whether refinement was needed; nothing is written otherwise. */

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

 SFDCRBRelaxationSolver( void ) : SingleFlowDCRBendersSolver() {}

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

 /// set the double parameters of the Solver
 /** dblRelAcc [1e-6] controls closure of the relative gap in compute() and
  * stopping of the fallback refinement in compute(). dblAbsAcc is stored;
  * the other parameters are passed to SingleFlowDCRBendersSolver. */

 void set_par( idx_type par , double value ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] double get_dbl_par( idx_type par ) const override;

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
  * scheme finds the node infeasible, kError if the Block has no source or
  * no sink or if the scheme ends in error, and kOK otherwise. */

 int compute( bool changedvars = true ) override;

/** @} ---------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Accessing the found solutions (if any)
 *  @{ */

 /// return a valid lower bound on the cost of the flow in the node
 /** The lower bound of SingleFlowDCRBendersSolver, or the independent
  * certificate if it is larger. The latter is computed only as a fallback
  * in compute(), by bisecting the domain of r_min and taking the minimum of
  * the lower bounds on the resulting intervals. compute() reports LB == UB
  * when a finite upper bound is reached within dblRelAcc.
  * +Inf if the node is infeasible, -Inf if the last compute() did not
  * solve it. */

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
 /** Uses the pair cached by compute(), which probes the envelope minima
  * of up to eight BenBound intervals with the smallest lower bounds,
  * then BenBound's reported r_min, then the endpoints and midpoints of
  * those intervals. Duplicate candidates are tried only once. If no pair
  * is found, refine_lb() can close the gap or provide further candidates.
  * A node closed by refinement is reported to the driver by compute(),
  * so branch() need not be called for it. The two children close,
  * respectively, the arcs of the two paths
  * that leave the last node they have in common. A simple path uses at most
  * one of the two arcs, hence every routing of the node is in one of the
  * children. The child that closes the arc of the path not meeting the
  * deadline comes first, since it keeps the other path. An exception is
  * thrown if the node has no relaxation, or if the two paths cannot be
  * found. */

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

/*--------------------------------------------------------------------------*/
/*---------------------------- PROTECTED FIELDS ----------------------------*/
/*--------------------------------------------------------------------------*/

 std::vector< bool > v_closed;  ///< the arcs closed by the Solver

 bool f_infeasible = false;     ///< the node has no path on the open arcs

 /// the certificate of the node [see get_lb()], -Inf if there is none
 OFValue f_cert = - SMSpp_di_unipi_it::Inf< OFValue >();

 double f_RelAcc = 1e-6;        ///< the parameter dblRelAcc

 double f_AbsAcc = 1e-6;        ///< the parameter dblAbsAcc

 /// the middle point of the r_min interval that gives the lower bound
 /** Set by the fallback refine_lb(), 0 if it has not been run;
  * prepare_branch() probes it after the BenBound candidates have failed,
  * unless refinement has already closed the node. */

 double f_zstar = 0;

 double f_za = 0; ///< the left end of the r_min interval that gives the bound

 int f_MaxIv = 200;   ///< the maximum number of bisections of refine_lb()

 /// distinct paths cached by compute() for branch(), empty if none found
 std::vector< int > v_branch_plus , v_branch_minus;

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

private:

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE METHODS ------------------------------*/
/*--------------------------------------------------------------------------*/

 /// cache a branching pair, refining the bound only if the probes fail
 /** Called by compute() before its final gap check, so a node closed by
  * refinement is pruned by the driver before branch() is requested. */

 void prepare_branch( void );

/*--------------------------------------------------------------------------*/

 /// lower bound on the cost of the flows whose r_min is in [ a , b ]
 /** Solves the Lagrangian relaxation in which every rate is at least a and
  * the burst term is that of b, which contains all the routings with r_min
  * in [ a , b ]; +Inf if it is empty, -Inf if nothing can be said. */

 double interval_lb( double a , double b );

/*--------------------------------------------------------------------------*/
 /// refine the lower bound of the node by bisecting the domain of r_min
 /** Bisects [ max{ rho , burst / deadline } , largest capacity ] keeping the
  * interval with the smallest lower bound, until it is not below the upper
  * bound, it is narrow, or f_MaxIv bisections have been done; the least
  * lower bound found goes in f_cert, and the middle point of the interval
  * that attains it in f_zstar. */

 void refine_lb( void );

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
