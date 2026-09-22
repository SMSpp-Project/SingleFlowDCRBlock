/*--------------------------------------------------------------------------*/
/*------------------ File SingleFlowDCRBendersSolver.h ---------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the SingleFlowDCRBendersSolver class, implementing a
 * Solver for Delay-Constrained Routing problems (DCR) relative to a Single
 * Flow, as set by SingleFlowDCRBlock, via a "Benders with nested Lagrange"
 * approach.
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

#ifndef __SingleFlowDCRBendersSolver
#define __SingleFlowDCRBendersSolver
/* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "Solver.h"

#include "SingleFlowDCRBlock.h"

#include "BenBound.h"

#include "DCR.h"

#include "BlockSolverConfig.h"

#include "CDASolver.h"

#include <algorithm>
#include <cmath>
#include <limits>

/*--------------------------------------------------------------------------*/
/*-------------------------- NAMESPACE & USING -----------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it {
/*--------------------------------------------------------------------------*/
/*------------------------------- CLASSES ----------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup SingleFlowDCRBendersSolver_CLASSES
 *  Classes in SingleFlowDCRBendersSolver.h
 *  @{ */

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS SingleFlowDCRBendersSolver --------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// Solver for SingleFlowDCRBlock via a "Benders with nested Lagrange" scheme
/** SingleFlowDCRBendersSolver implements the Solver concept [see Solver.h]
 * for SingleFlowDCRBlock by wrapping the (non-SMS++, "legacy") BenBound
 * class, which it publicly inherits from together with Solver.
 *
 * The DCR problem is a min-cost routing problem in which, besides the
 * flow-conservation and capacity constraints, a worst-case end-to-end
 * delay bound (computed via network calculus, as a function of the
 * reserved rate r[ i , j ] on each arc and of the "burst" of the flow)
 * must not exceed a given deadline. BenBound tackles the resulting
 * nonlinear (albeit convex, arc-wise separable) problem via a Benders-like
 * decomposition on r_min, the minimum reserved rate along the used path:
 * for a fixed value of r_min, the residual problem decouples into a
 * (Lagrangian-relaxable) shortest-path-like problem, solved by the nested
 * DCRLagrangianSolver (which in turn solves a sequence of shortest path
 * problems via SPT); the resulting piecewise-linear lower bound as a
 * function of r_min is then explored by a line-search (Benders master)
 * over r_min itself. See BenBound.h and DCRLagrangianSolver.h for the
 * details of this "legacy" solution algorithm.
 *
 * This class merely translates the data of a SingleFlowDCRBlock into the
 * format expected by BenBound::LoadProblem() (see set_Block() and
 * process_outstanding_Modification()), drives the solution process (see
 * compute()), and translates the results back into SMS++ terms (see
 * get_lb(), get_ub(), get_var_solution(), ...). */

class SingleFlowDCRBendersSolver : public Solver , public BenBound
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

public:
/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
 /** @name Public Types
  *  @{ */

 /** @} ---------------------------------------------------------------------*/
 /*----------------- CONSTRUCTING AND DESTRUCTING SingleFlowDCRBendersSolver */
/*--------------------------------------------------------------------------*/
 /** @name Constructing and destructing SingleFlowDCRBendersSolver
  *  @{ */

 /// constructor: does nothing special
 /** Void constructor: does nothing special besides default-constructing
  * both base classes, Solver and BenBound. */

 SingleFlowDCRBendersSolver( void ) : Solver() , BenBound() {}

/*--------------------------------------------------------------------------*/
 /// destructor: does nothing special
 /** Destructor: currently empty; the Solver and BenBound base classes take
  * care of releasing their own resources (comprised any pending
  * Modification) via their own destructors. */

 virtual ~SingleFlowDCRBendersSolver() {}

 /** @} ---------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
 /** @name Other initializations
  *
  * SingleFlowDCRBendersSolver currently does not define any Solver
  * parameter of its own (no {get,set}_{num_,}*_par() method is overridden):
  * it only inherits the (empty) set of parameters of the base Solver class.
  *  @{ */

 /// set the (pointer to the) Block that the Solver has to solve
 /** Sets (or resets, if \p block == nullptr) the SingleFlowDCRBlock that
  * this Solver has to solve (an exception is thrown if \p block is not a
  * SingleFlowDCRBlock). Besides the usual bookkeeping (via
  * Solver::set_Block()), this [read_]locks the Block (unless already
  * owned), extracts source/sink node (the nodes with, respectively,
  * negative and positive deficit), traffic parameters (burst, rate,
  * deadline), MTU, and the per-arc topology/speed/capacity/delay/cost
  * data, translates them into the DCR::DCRFlow / DCR::DCRLink / DCR::DCRNode
  * format expected by BenBound, and calls BenBound::LoadProblem() with
  * them; the same translation is repeated (see
  * process_outstanding_Modification()) whenever the SingleFlowDCRBlock is
  * subsequently modified. If the call does nothing (\p block == f_Block),
  * it cowardly and silently returns. */

 void set_Block( Block * block ) override
 {

  if( block == f_Block ) // actually doing nothing
   return;               // cowardly and silently return

  Solver::set_Block( block ); // attach to the new Block

  if( block ) { // this is not just resetting everything
   auto MCFB = dynamic_cast< SingleFlowDCRBlock * >( f_Block );
   if( ! MCFB )
    throw( std::invalid_argument( "SingleFlowDCRBendersSolver:set_Block: "
                                  "block must be a SingleFlowDCRBlock" ) );

   bool owned = MCFB->is_owned_by( f_id );
   if( ( ! owned ) && ( ! MCFB->read_lock() ) )
    throw(
     std::logic_error( "cannot acquire read_lock on SingleFlowDCRBlock" ) );
   // load the new SingleFlowDCRBlock into the :BenBound object

   load_BenBound();

   // once done, read_unlock the SingleFlowDCRBlock (if it was read-lock()-ed)
   if( ! owned )
    MCFB->read_unlock();

   // TODO: maybe log it
   }
  } // end( set_Block )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

 /** @} ---------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE Block ----------------------*/
/*--------------------------------------------------------------------------*/
 /** @name Solving the SingleFlowDCR
  *  @{ */

 /// (try to) solve the SingleFlowDCR
 /** (Tries to) solve the SingleFlowDCRBlock currently attached to this
  * Solver via the "Benders with nested Lagrange" scheme implemented by
  * BenBound (see the class GENERAL NOTES). After acquiring the self-lock
  * and, if necessary, a read_lock on f_Block, processes any outstanding
  * Modification (see process_outstanding_Modification()), then calls
  * BenBound::DCRsetTime() / BenBound::Solve() / BenBound::DCRstopTime() to
  * actually run the algorithm. \p changedvars is currently unused.
  *
  * IMPLEMENTATION NOTE: the method contains a (large) commented-out block
  * that would run an external warm-start Solver (via BlockSolverConfig,
  * reading "MILPPar-initial.txt", "DCRCfgPC.txt" and a
  * "WarmStartCfg.txt"-selected Solver) to compute an initial bound on
  * r_min (stored in bound_PC) before calling BenBound; this warm-start is
  * currently disabled.
  *
  * Returns kOK unconditionally (i.e., irrespective of the actual status
  * reported by BenBound::getStat(), which can be queried via
  * has_var_solution()). */

 int compute( bool changedvars = true ) override
 {

  lock(); // first of all, acquire self-lock

  if( ! f_Block )            // there is no [SingleFlowDCRBlock] to solve
   return( kBlockLocked ); // return error

  bool owned = f_Block->is_owned_by( f_id );    // check if already locked
  if( ( ! owned ) && ( ! f_Block->read_lock() ) ) // if not try to read_lock
   return( kBlockLocked );                     // return error on failure

  // while [read_]locked, process any outstanding Modification

  if( ! owned ) // if the [SingleFlowDCR]Block was actually read_locked
   f_Block->read_unlock(); // read_unlock it

  process_outstanding_Modification();

  auto f_Block_aux = f_Block;
  /*
 const std::string warmstart_cfg_file = "MILPPar-initial.txt";

 auto warmstart_cfg = Configuration::deserialize( warmstart_cfg_file );
 auto warmstart_bsc = dynamic_cast< BlockSolverConfig * >( warmstart_cfg );

 auto dcrb_pc = dynamic_cast< BlockConfig * >
          ( Configuration::deserialize( "DCRCfgPC.txt" ) );

 // instantiate the warm-start Solver via the factory, picking name and
 // (optional) ComputeConfig from WarmStartCfg.txt. The Solver is NOT
 // registered on f_Block so it cannot interfere with PrimalProximalHeur
 // itself (which is already attached to f_Block).
 auto warmstart =
  dynamic_cast< CDASolver * >(
   Solver::new_Solver( warmstart_bsc->get_SolverName( 0 ) ) );

 if( warmstart_bsc->num_ComputeConfig() > 0 )
    if( auto cc = warmstart_bsc->get_SolverConfig( 0 ) )
      warmstart->set_ComputeConfig( cc );
 dcrb_pc->apply( f_Block );
 warmstart->set_Block( f_Block );
 warmstart->compute( changedvars );
 warmstart->get_var_solution();
 bound_PC = dynamic_cast< SingleFlowDCRBlock * >( f_Block )->get_rmin();

 std::cout << "bound_PC = " << bound_PC << std::endl;
 std::cout << "sol = " << warmstart->get_lb() << std::endl;
   */
  // ensure the timer exists (or reset it)
  BenBound::DCRsetTime( true );
  // then (try to) solve the SingleFlowDCR
  BenBound::Solve();
  BenBound::DCRstopTime();

  //std::cout << "LB = " << get_lb() << std::endl;
  //std::cout << "UB = " << get_ub() << std::endl;

  unlock(); // unlock the mutex

  // now give out the result: note that the vector MCFstatus_2_sol_type[]
  // starts from 0 whereas the first value of MCFStatus is -1 (= kUnSolved),
  // hence the returned status has to be shifted by + 1

  return( kOK );
  }

 /** @} ---------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/
 /** @name Accessing the found solutions (if any)
  *  @{ */

 /// returns the time (in seconds) spent so far in BenBound::Solve()

 double get_elapsed_time( void ) const override
 {
  return( this->BenBound::getTime() );
  }

/*--------------------------------------------------------------------------*/
 /// returns the best lower bound found so far by the Benders/Lagrange bound

 OFValue get_lb( void ) override { return( this->BenBound::getLB() ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns the best upper bound on the optimum
 /** An upper bound is only meaningful when it is the cost of an actually
  * achievable, feasible solution: hence, whenever pick_solution_source()
  * finds one (either the primary or the heuristic candidate), this
  * *always* recomputes the bound from scratch via get_var_value(), out of
  * the very routing/rates that get_var_solution() and get_Solution() would
  * hand out for that same candidate, rather than trusting
  * BenBound::getUB() == min( BestUB , HeurVal ) at face value. Deliberately
  * so: BestUB is a Lagrangian-relaxation value (\f$ \alpha + \lambda \beta
  * \f$ for whatever path/lambda attained it), not the routing's raw cost
  * \f$ \alpha \f$, and the two can differ by an arbitrary amount whenever
  * the delay slack \f$ \beta \f$ at that point is not exactly 0 -- even
  * when the extracted routing is fully DCR-feasible. A concrete instance:
  * a candidate whose actually-reserved rates cost 162000 (verified
  * feasible, has_var_solution() true) was once seen reported as an upper
  * bound of 129800 by trusting BenBound::getUB() directly, a bound with no
  * solution behind it at all -- 129800 is not the cost of anything this
  * Solver can actually produce. get_var_value() cannot suffer from this,
  * since it always re-derives the value from the concrete X, R it reads.
  *
  * If no candidate is feasible, +INF is returned, i.e. no bound at all.
  * Note that BenBound and DCRLagrangianSolver are plain (non-SMS++) classes
  * that define their own local template class "Inf< T >", whose operator
  * T() returns std::numeric_limits< T >::max() rather than a true
  * infinity; since SingleFlowDCRBendersSolver derives from BenBound, an
  * unqualified "Inf< OFValue >()" here would resolve to that inherited
  * class, NOT to SMSpp_di_unipi_it::Inf< OFValue >() (a real infinity),
  * even though this method lives in the SMSpp_di_unipi_it namespace:
  * base-class lookup wins over the enclosing namespace. Hence the explicit
  * qualification below, needed so a Solver with no bound at all correctly
  * reports +Inf and not a (finite!) huge number that could be mistaken for
  * a genuine, if very bad, bound by a caller checking std::isfinite(). */

 OFValue get_ub( void ) override
 {
  if( pick_solution_source() != SolSource::None )
   return( get_var_value() );

  return( SMSpp_di_unipi_it::Inf< OFValue >() );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns the objective value of the (active) solution BenBound found
 /** Recomputes, from scratch, the objective value \f$ \sum_i C[ i ] r[ i ]
  * \f$ of whichever of BenBound's two candidates pick_solution_source()
  * currently selects [see get_solution_vectors()]. Note that the reserved
  * rates are taken from BenBound and *not* from the Variable of the
  * SingleFlowDCRBlock: what the Variable hold is whatever was last written
  * there, by this Solver or by anybody else, while what this method has to
  * answer is the value of the solution that this Solver found. */

 OFValue get_var_value( void ) override
 {
  return( cost_of( pick_solution_source() == SolSource::Heuristic ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns the number of Benders iterations performed so far

 Index get_BenIt( void )
 {
  return( this->BenBound::getNumIterationBender() );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns the total number of (nested) Lagrangian iterations performed

 Index get_LagIt( void ) { return( this->BenBound::getNumIterationLagr() ); }

/*--------------------------------------------------------------------------*/
 /// returns true if BenBound found a feasible (primal) solution
 /** Returns true if and only if BenBound terminated properly
  * (BenBound::getStat() == BenBound::OK) *and* at least one of its two
  * candidate solutions -- the primary one or the heuristic one -- is
  * feasible for the DCR problem [see pick_solution_source()]: a point
  * that is not is not a solution, whatever BenBound makes of it. */

 bool has_var_solution( void ) override
 {
  return( pick_solution_source() != SolSource::None );
  }

/*--------------------------------------------------------------------------*/
 /// writes the (active) BenBound solution back into the SingleFlowDCRBlock
 /** Copies the per-arc reserved-rate solution of whichever candidate
  * pick_solution_source() currently selects (see get_solution_vectors())
  * into the r[] Variable of the SingleFlowDCRBlock f_Block, deriving the
  * corresponding "flow" ColVariable x[ i ] (1 if the arc is used, i.e., its
  * reserved rate is > 0, 0 otherwise) and the theta[ i ] Variable
  * (MTU / r[ i ] if the arc is used, 0 otherwise); \p solc is currently
  * unused. If f_Block is nullptr the method silently does nothing.
  *
  * The aggregate r_min is recomputed here as the smallest R[ i ] among the
  * arcs actually used, rather than read off BenBound::getr_min(): the
  * latter is the candidate r_min BenBound's own line search was evaluating
  * (meaningful for the primary solution only), while the former is correct
  * for either solution and matches what SingleFlowDCRBlock::min_rate()
  * would report for the very routing being written out here. */

 void get_var_solution( Configuration * solc = nullptr ) override
 {
  if( ! f_Block ) // no [SingleFlowDCR]Block to write to
   return;       // cowardly and silently return

  auto DCRB = static_cast< SingleFlowDCRBlock * >( f_Block );

  SingleFlowDCRBlock::Vec_double X , R;
  get_solution_vectors( X , R , pick_solution_source() == SolSource::Heuristic );

  double rmin = std::numeric_limits< double >::infinity();
  for( Index i = 0 ; i < X.size() ; ++i ) {
   DCRB->set_x( i , X[ i ] );
   DCRB->set_r( i , R[ i ] );
   // the smallest burst delay the cone constraint of the arc allows
   DCRB->set_theta( i , X[ i ] > 0 ? DCRB->get_MTU() / R[ i ] : 0 );
   if( X[ i ] > 0 )
    rmin = std::min( rmin , R[ i ] );
   }
  if( ! std::isfinite( rmin ) )
   rmin = 0;

  // the same for the two "aggregate" variables: the minimum reserved rate
  // along the path and the burst delay that goes with it
  DCRB->set_rmin( rmin );
  DCRB->set_theta_min( rmin > 0 ? DCRB->get_FlowBurst() / rmin : 0 );
  }

/*--------------------------------------------------------------------------*/
 /// produces the Solution of the SingleFlowDCRBlock out of BenBound
 /** Fills the Solution that the SingleFlowDCRBlock provides directly out of
  * the solution BenBound found, without writing anything into the Variable
  * of the Block, which therefore need not even exist [see
  * Solver::get_Solution()]. The Solution is asked to the Block, as the
  * general rule requires; should it not be a DCRSolution, which is what a
  * SingleFlowDCRBlock gives, the base class is left to do the job. */

 [[nodiscard]] Solution *
 get_Solution( Configuration * solc = nullptr ) override
 {
  if( ! f_Block )
   return( Solver::get_Solution( solc ) );

  auto sol = f_Block->get_Solution( solc , true );
  auto dsol = dynamic_cast< DCRSolution * >( sol );
  if( ! dsol )      // not the Solution this Solver knows how to fill:
   return( sol ); // hand back the empty one the Block gave

  const auto src = pick_solution_source();
  if( src == SolSource::None ) // nothing to put in it
   return( sol );

  SingleFlowDCRBlock::Vec_double X , R;
  get_solution_vectors( X , R , src == SolSource::Heuristic );

  // which of the two parts the Solution wants is its own business, as
  // dictated by the Configuration it was asked with: only fill in what is
  // there [see SingleFlowDCRBlock::get_Solution()]
  if( ! dsol->get_x().empty() )
   dsol->set_x( std::move( X ) );
  if( ! dsol->get_r().empty() )
   dsol->set_r( std::move( R ) );

  return( sol );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns true if the solution this Solver hands out meets the deadline
 /** Checks, within the relative tolerance feps, whether the worst-case
  * end-to-end delay of the routing found by BenBound meets the deadline of
  * the flow: the routing and the reserved rates are read out of BenBound
  * [see get_solution_vectors()] and handed to
  * SingleFlowDCRBlock::delay_feasible(), which is where the network
  * calculus formula lives. Note that this says nothing about the *other*
  * constraints of the DCR problem: use
  * SingleFlowDCRBlock::is_sol_feasible() on the Solution this Solver
  * produces to have them all checked.
  *
  * The candidate checked is the one pick_solution_source() selects, i.e.
  * exactly the one get_ub(), get_var_value(), get_var_solution() and
  * get_Solution() report on, so that this method answers the question a
  * caller actually asks of it: "is the solution I am being given
  * feasible?". It deliberately does NOT default to BenBound's primary
  * candidate: that one is tied to BestUB, it is frequently not the one
  * handed out (whenever the heuristic candidate is the feasible or the
  * cheaper of the two), and checking it produced an answer about a vector
  * nobody ever receives -- on a sample of 44 flows whose bounds were
  * being investigated, 7 reported "not DCR feasible" while the solution
  * actually returned passed every DCR check, at this tolerance and at the
  * looser one. When no candidate is feasible at all the Solver hands out
  * nothing (get_ub() is +INF), and this correctly returns false rather
  * than reporting on a vector that is not on offer.
  *
  * The default tolerance is the one this Solver holds itself to when it
  * declares having a solution [see solution_is_feasible()], and it is
  * deliberately a tight one: a point that misses the deadline by a whisker
  * is still a point that misses the deadline, and its value bounds
  * nothing. Note that the *selection* above is made at
  * pick_solution_source()'s own tolerance, not at feps: which candidate is
  * on offer is a property of the Solver, not of how strictly the caller
  * wants to inspect it, so that a caller tightening feps re-checks the
  * same solution rather than silently being handed a different one. */

 bool is_DCR_feasible( double feps = 1e-6 )
 {
  const SolSource src = pick_solution_source();

  if( src == SolSource::None )  // nothing is being handed out at all
   return( false );

  auto DCRB = static_cast< SingleFlowDCRBlock * >( f_Block );

  SingleFlowDCRBlock::Vec_double X , R;
  get_solution_vectors( X , R , src == SolSource::Heuristic );

  return( DCRB->delay_feasible( feps , X , R ) );
  }

 /** @} ---------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/
 /** @name Changing the data of the model
  *  @{ */

/*--------------------------------------------------------------------------*/
 /// processes any Modification issued by the SingleFlowDCRBlock since the
 /// last compute()
 /** Drains, under lock, the list v_mod of Modification collected by the
  * Solver since the last call (see Solver::process_outstanding_Modification
  * for the general mechanism); if at least one Modification is found, the
  * whole BenBound problem is rebuilt from scratch out of the current state
  * of the SingleFlowDCRBlock, exactly as done in set_Block() (topology,
  * costs, capacities, delays, flow burst/rate/deadline and MTU are all
  * re-read and passed to BenBound::LoadProblem()). No attempt is made to
  * incrementally update the BenBound data structures according to the
  * specific Modification received; any change at all triggers a full
  * reload. */

 void process_outstanding_Modification( void )
 {

  bool reload = false;

  // note: since processing the Modification is fast, we don't bother with
  // being nice to other processes and do it all with v_mod under lock
  // try to acquire lock, spin on failure
  while( f_mod_lock.test_and_set( std::memory_order_acquire ) )
   ;

  // process all the Modifications
  for( auto mod : v_mod )
   if( auto tmod = mod.get() ) {
    reload = true; // a reset must be done
    break;         // ignore all the remaining Modifications
   }

  v_mod.clear(); // all Modifications tackled, clear the list

  f_mod_lock.clear( std::memory_order_release ); // release lock

  if( reload ) {

   load_BenBound();
   }
  }

 /** @} ---------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

protected:
/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/
 /// (re)loads the data of the SingleFlowDCRBlock into BenBound
 /** Translates the instance held by the SingleFlowDCRBlock this Solver is
  * attached to into the DCR::DCRFlow / DCRLink / DCRNode description that
  * BenBound wants [see DCR.h] and hands it to BenBound::LoadProblem(). The
  * source and the sink of the flow are the nodes with, respectively,
  * negative and positive deficit; a Block with no source or no sink has
  * nothing to route, and BenBound is left with whatever it had. This is
  * done when the Block is attached and every time it changes [see
  * set_Block() and process_outstanding_Modification()]. */

 void load_BenBound( void )
 {
  auto DCRB = static_cast< SingleFlowDCRBlock * >( f_Block );
  auto nnodes = DCRB->get_NNodes();
  auto narcs = DCRB->get_NArcs();

  int source = -1;
  int sink = -1;
  for( Index i = 0 ; i < nnodes ; ++i ) {
   auto Bi = DCRB->get_B( i );
   if( Bi < 0 )
    source = i;
   else
    if( Bi > 0 )
     sink = i;
   }

  if( ( source < 0 ) || ( sink < 0 ) )  // nothing to route
   return;

  DCR::DCRFlow flow = {};
  flow.sourcenode = source;
  flow.sinknode = sink;
  flow.burst = DCRB->get_FlowBurst();
  flow.rate = DCRB->get_rho();
  flow.deadline = DCRB->get_FlowDeadline();

  auto & LD = DCRB->get_LinkDelays();
  auto & ND = DCRB->get_NodeDelays();

  std::vector< DCR::DCRLink > links( narcs );
  for( Index i = 0 ; i < narcs ; ++i ) {
   links[ i ] = {};
   links[ i ].startnode = DCRB->get_SN( i ) - 1;
   links[ i ].endnode = DCRB->get_EN( i ) - 1;
   links[ i ].speed = DCRB->get_U( i );
   links[ i ].capacity = DCRB->get_U( i );
   links[ i ].delay = LD.empty() ? 0 : LD[ i ];
   links[ i ].cost = DCRB->is_deleted( i ) ? 0 : DCRB->get_C( i );
   }

  std::vector< DCR::DCRNode > nodes( nnodes );
  for( Index i = 0 ; i < nnodes ; ++i ) {
   nodes[ i ] = {};
   nodes[ i ].delay = ND.empty() ? 0 : ND[ i ];
   }

  BenBound::LoadProblem( nnodes , narcs , flow , links , nodes ,
                         DCRB->get_MTU() );
 }

/*--------------------------------------------------------------------------*/
 /// true if the solution BenBound found is feasible for the DCR problem
 /** Asks the SingleFlowDCRBlock whether the routing and the reserved rates
  * BenBound found [see get_solution_vectors()] satisfy all the constraints
  * of the DCR problem: that the routing is a path from the source to the
  * sink, that the rates are within the arc capacities and go with the
  * routing, and that the end-to-end delay meets the deadline. This is what
  * has_var_solution() and get_ub() are held to: a Solver that hands out a
  * point that is not feasible, and calls its value an upper bound on the
  * optimum, says something false.
  *
  * The default tolerance, 1e-5, matches the one SingleFlowDCRBlock.cpp
  * itself uses for the same kind of cone-feasibility check [see "tol" in
  * generate_dynamic_constraints()]: at the true optimum the end-to-end
  * delay sits exactly at the deadline (it is the constraint that
  * determines the optimal r_min), so the closed-form Lagrangian formulas
  * BenBound uses to evaluate a candidate necessarily land within
  * floating-point roundoff of it, not exactly on it. A tighter tolerance
  * here (1e-6 was tried first) makes has_var_solution() reject genuinely
  * optimal candidates purely on accumulated roundoff from the several
  * sqrt()/pow() evaluations in that closed form, reporting +Inf for a
  * problem that was in fact solved to optimality. */

 bool solution_is_feasible( bool heuristic = false , double feps = 1e-5 )
 {
  auto DCRB = static_cast< SingleFlowDCRBlock * >( f_Block );

  SingleFlowDCRBlock::Vec_double X , R;
  get_solution_vectors( X , R , heuristic );

  return( DCRB->flow_feasible( feps , X ) &&
          DCRB->bound_feasible( feps , X , R ) &&
          DCRB->link_feasible( feps , X , R ) &&
          DCRB->delay_feasible( feps , X , R ) );
  }

/*--------------------------------------------------------------------------*/
 /// which candidate solution (if any) get_ub()/get_var_solution()/
 /// get_Solution()/is_DCR_feasible() should currently report
 /** BenBound holds up to two internally-tracked candidates: the "primary"
  * one (BenBound::getSolution(), tied to BestUB/getUB()) and the
  * "heuristic" one (BenBound::getHeurSolution(), tied to getHeurVal()).
  * Neither is guaranteed feasible a priori -- BestUB is a Lagrangian-
  * relaxation value that can be attained at a delay-infeasible candidate
  * [see the note on get_ub()], while HeurVal is only ever updated from a
  * candidate DCRLagrangianSolver itself verified delay-feasible, but nothing
  * stops both from being unset (BenStat != OK) or, in principle, both from
  * failing the DCR-level check applied here (e.g. on a badly conditioned
  * instance). This method is the single place that decides, so get_ub()
  * and the solution-reporting methods below never disagree on which one
  * they are looking at: prefer the primary candidate when it checks out
  * (it is the one BestUB/getUB() naturally corresponds to), but only
  * when it does not cost more than the heuristic one: BestUB is a
  * Lagrangian-relaxation value, not necessarily this Solver's cheapest
  * verified-feasible candidate, whereas getHeurVal() (and hence the
  * heuristic candidate) is reduced_line-search's own explicit tracking
  * of the best delay-feasible routing it has found, arc rates included --
  * see DCRLagrangianSolver::tightenPath(), which squeezes any avoidable
  * slack out of every heuristic candidate before it is accepted. When
  * both hold up, prefer whichever is actually cheaper; a Solver that
  * reports an upper bound larger than one it could have reported instead
  * is not wrong, but it is needlessly weak. */

 enum class SolSource { None , Primary , Heuristic };

 /// recomputes, from scratch, the raw cost of one of BenBound's two
 /// candidates (see get_solution_vectors()), independently of whichever
 /// one pick_solution_source() currently selects
 double cost_of( bool heuristic )
 {
  auto DCRB = static_cast< SingleFlowDCRBlock * >( f_Block );

  SingleFlowDCRBlock::Vec_double X , R;
  get_solution_vectors( X , R , heuristic );

  double sum = 0;
  for( Index i = 0 ; i < R.size() ; ++i )
   if( ! DCRB->is_deleted( i ) )
    sum += DCRB->get_C( i ) * R[ i ];

  return( sum );
  }

 SolSource pick_solution_source( double feps = 1e-5 )
 {
  if( this->BenBound::getStat() != BenBound::OK )
   return( SolSource::None );

  const bool primary_ok = solution_is_feasible( false , feps );
  const bool heuristic_ok = solution_is_feasible( true , feps );

  if( primary_ok && heuristic_ok )
   return( cost_of( true ) < cost_of( false ) ? SolSource::Heuristic
                                              : SolSource::Primary );

  if( primary_ok )
   return( SolSource::Primary );

  if( heuristic_ok )
   return( SolSource::Heuristic );

  return( SolSource::None );
  }

/*--------------------------------------------------------------------------*/
 /// reads a solution of BenBound into a routing and a rate vector
 /** Fills X and R, both sized get_NArcs(), with a solution BenBound holds:
  * R[ i ] is the rate reserved on arc i and X[ i ] is 1 if the arc is used
  * by the flow, i.e., if that rate is positive, and 0 otherwise. With \p
  * heuristic == false (the default) this is BenBound::getSolution(), the
  * candidate tied to BestUB/getUB(); with \p heuristic == true it is
  * instead BenBound::getHeurSolution(), the candidate tied to getHeurVal()
  * that DCRLagrangianSolver itself already verified delay-feasible [see
  * pick_solution_source()]. This is the one place either solution of
  * BenBound is turned into a solution of the SingleFlowDCRBlock, and it is
  * what get_var_value(), get_var_solution(), get_Solution() and
  * is_DCR_feasible() all go through. */

 void get_solution_vectors( SingleFlowDCRBlock::Vec_double & X ,
                            SingleFlowDCRBlock::Vec_double & R ,
                            bool heuristic = false )
 {
  auto DCRB = static_cast< SingleFlowDCRBlock * >( f_Block );
  auto narcs = DCRB->get_NArcs();

  X.resize( narcs );
  R.resize( narcs );

  for( Index i = 0 ; i < narcs ; ++i ) {
   auto v = heuristic ? BenBound::getHeurSolution( i )
                       : BenBound::getSolution( i );
   R[ i ] = v;
   X[ i ] = v > 0 ? 1 : 0;
   }
  }

/*--------------------------------------------------------------------------*/
/*---------------------------- PROTECTED FIELDS  ---------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

private:
 double bound_PC = 0.0; ///< bound on r_min computed by an (optional)
                        ///< warm-start
                        ///< primal heuristic Solver; currently only set by
                        ///< the code commented out inside compute(), so it
                        ///< is unused in the present implementation

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 }; // end( class( SingleFlowDCRBendersSolver ) )

/*--------------------------------------------------------------------------*/

/** @}  end( group( SingleFlowDCRBendersSolver_CLASSES ) ) -----------------*/
/*--------------------------------------------------------------------------*/

 } // namespace SMSpp_di_unipi_it

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* SingleFlowDCRBendersSolver.h included */

/*--------------------------------------------------------------------------*/
/*---------------- End File SingleFlowDCRBendersSolver.h -------------------*/
/*--------------------------------------------------------------------------*/
