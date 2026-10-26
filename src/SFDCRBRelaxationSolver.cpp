/*--------------------------------------------------------------------------*/
/*------------------- File SFDCRBRelaxationSolver.cpp ----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the SFDCRBRelaxationSolver class, the RelaxationSolver
 * of a SingleFlowDCRBlock that solves the nodes by means of
 * SingleFlowDCRBendersSolver and branches on the arcs.
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
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "SFDCRBRelaxationSolver.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <queue>

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register SFDCRBRelaxationSolver to the Solver factory

SMSpp_insert_in_factory_cpp_0( SFDCRBRelaxationSolver );

/*--------------------------------------------------------------------------*/
/*--------------------------- AUXILIARY CLASSES ----------------------------*/
/*--------------------------------------------------------------------------*/

namespace {

using ProfileClock = std::chrono::steady_clock;

double profile_seconds( ProfileClock::time_point start )
{
 return( std::chrono::duration< double >( ProfileClock::now() - start ).count() );
}

}  // end( anonymous namespace )

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS OF SFDCRBRelaxationSolver -----------------*/
/*--------------------------------------------------------------------------*/

void SFDCRBRelaxationSolver::set_Block( Block * block )
{
 if( block && ( ! dynamic_cast< SingleFlowDCRBlock * >( block ) ) )
  throw( std::invalid_argument( "SFDCRBRelaxationSolver::set_Block: block "
                                "must be a SingleFlowDCRBlock" ) );

 invalidate_solution_cache();
 f_instance_fresh = false;
 v_closed.clear();  // before the base loads the instance [see arc_closed()]
 f_infeasible = false;
 f_cert = - SMSpp_di_unipi_it::Inf< OFValue >();
 f_zstar = 0;
 f_za = 0;
 v_bound_intervals.clear();
 f_bound_only = true;
 v_branch_plus.clear();
 v_branch_minus.clear();

 SingleFlowDCRBendersSolver::set_Block( block );

 if( block )
  v_closed.assign( static_cast< SingleFlowDCRBlock * >( block )->get_NArcs() ,
                   false );

 }  // end( SFDCRBRelaxationSolver::set_Block )

/*--------------------------------------------------------------------------*/

void SFDCRBRelaxationSolver::set_par( idx_type par , double value )
{
 if( par == dblRelAcc )
  f_RelAcc = value;
 else if (par == dblAbsAcc)
  f_AbsAcc = value;
 else if( par == dblUpCutOff )
  f_UpCutOff = value;
 else
  SingleFlowDCRBendersSolver::set_par( par , value );

 }  // end( SFDCRBRelaxationSolver::set_par( double ) )

/*--------------------------------------------------------------------------*/

double SFDCRBRelaxationSolver::get_dbl_par( idx_type par ) const
{
 if( par == dblRelAcc )
  return( f_RelAcc );
 else if (par == dblAbsAcc)
  return( f_AbsAcc );
 else if( par == dblUpCutOff )
  return( f_UpCutOff );
 return( SingleFlowDCRBendersSolver::get_dbl_par( par ) );

 }  // end( SFDCRBRelaxationSolver::get_dbl_par )

/*--------------------------------------------------------------------------*/

void SFDCRBRelaxationSolver::load_BenBound()
{
 invalidate_solution_cache();
 f_bound_only = true;
 SingleFlowDCRBendersSolver::load_BenBound();
}

void SFDCRBRelaxationSolver::add_Modification( sp_Mod & mod )
{
 SingleFlowDCRBendersSolver::add_Modification( mod );
 invalidate_solution_cache();
}

SingleFlowDCRBendersSolver::SolSource
SFDCRBRelaxationSolver::pick_solution_source( double feps )
{
 const auto revision = f_solution_revision.load( std::memory_order_relaxed );
 if( revision == f_cached_revision && feps == f_solution_feps ) {
  if( f_log )
   ++f_profile.solution_cache_hits;
  return( f_solution_source );
  }
 if( f_log )
  ++f_profile.solution_checks;
 const auto source = f_bound_only ? SolSource::None :
  SingleFlowDCRBendersSolver::pick_solution_source( feps );
 f_solution_source = source;
 f_solution_cost = source == SolSource::None ? 0 :
  cost_of( source == SolSource::Heuristic );
 f_solution_feps = feps;
 f_cached_revision = revision;
 return( source );
}

Solver::OFValue SFDCRBRelaxationSolver::get_var_value()
{
 const auto source = pick_solution_source();
 // Preserve the base accessor's raw-primary-cost behavior with no candidate.
 return( source == SolSource::None ? cost_of( false ) : f_solution_cost );
}

/*--------------------------------------------------------------------------*/

int SFDCRBRelaxationSolver::compute( bool changedvars )
{
 // Hints can make the outer search much more expensive even at the root.
 // Suspend them for this compute, restoring the configured option on every
 // return path (and on exceptions), without changing the base solver.
 struct RestoreWarmLambda {
  BenBound & solver;
  bool enabled;
  ~RestoreWarmLambda() { solver.setWarmLambda( enabled ); }
 } restore_warm{ *this , get_int_par( intWarmLambda ) != 0 };
 BenBound::setWarmLambda( false );

 invalidate_solution_cache();
 const auto start = f_log ? ProfileClock::now() : ProfileClock::time_point{};
 ++f_compute_id;
 f_profile = {};
 double setup_seconds = 0 , benbound_seconds = 0 , prepare_seconds = 0;
 OFValue lb_before = - SMSpp_di_unipi_it::Inf< OFValue >();
 OFValue ub_before = SMSpp_di_unipi_it::Inf< OFValue >();
 bool results_available = false;
 auto finish = [ & ]( int status , const char * reason ) {
  if( f_log ) {
   const double total = profile_seconds( start );
   const double prepare_other = std::max( 0.0 , prepare_seconds -
    f_profile.probe_seconds - f_profile.refine_seconds );
   f_log->precision( 12 );
   *f_log << "[node_profile] solver=" << this << " run=" << f_compute_id
    << " status=" << status << " reason=" << reason
    << " warm_lambda=0"
    << " total_s=" << total << " setup_s=" << setup_seconds
    << " benbound_s=" << benbound_seconds
    << " probes_s=" << f_profile.probe_seconds
    << " refine_s=" << f_profile.refine_seconds
    << " prepare_other_s=" << prepare_other
    << " other_s=" << std::max( 0.0 , total - setup_seconds -
                                   benbound_seconds - prepare_seconds )
    << " probes=" << f_profile.probes
    << " cache_hits=" << f_profile.cache_hits
    << " solution_checks=" << f_profile.solution_checks
    << " solution_cache_hits=" << f_profile.solution_cache_hits
    << " refine_calls=" << f_profile.refine_calls
    << " interval_solves=" << f_profile.interval_solves
    << " refine_splits=" << f_profile.refine_splits
    << " benIt=" << ( results_available ? BenBound::getNumIterationBender() : 0 )
    << " lagIt=" << ( results_available ? BenBound::getNumIterationLagr() : 0 )
    << " LB_before=" << lb_before
    << " LB=" << ( ( results_available || ! v_bound_intervals.empty() ) ? get_lb() :
       ( status == kInfeasible ? SMSpp_di_unipi_it::Inf< OFValue >() : lb_before ) )
    << " UB=" << ub_before << " cutoff=" << f_UpCutOff
    << " master_rmin=" << ( results_available ? branch_result().rmin : 0 )
    << " paths_available=" << ( results_available && branch_result().available )
    << " paths_equal=" << ( results_available && branch_result().available &&
                              branch_result().plus == branch_result().minus )
    << " dual_converged=" << ( results_available && branch_result().dual_converged )
    << " pair_ready=" << ( results_available &&
       ! v_branch_plus.empty() && ! v_branch_minus.empty() )
    << '\n';
   }
  return( status );
  };

 if( ! f_Block )
  return( finish( kError , "no_block" ) );

 f_cert = - SMSpp_di_unipi_it::Inf< OFValue >();
 f_zstar = 0;
 f_za = 0;
 v_bound_intervals.clear();
 f_bound_only = true;
 v_branch_plus.clear();
 v_branch_minus.clear();

 const bool owned = f_Block->is_owned_by( f_id );
 if( ( ! owned ) && ( ! f_Block->read_lock() ) )
  return( finish( kBlockLocked , "block_locked" ) );

 // a flow needs a source and a sink [see load_BenBound()]
 auto DCRB = static_cast< SingleFlowDCRBlock * >( f_Block );
 bool source = false;
 bool sink = false;
 for( Index i = 0 ; i < DCRB->get_NNodes() ; ++i ) {
  source |= ( DCRB->get_B( i ) < 0 );
  sink |= ( DCRB->get_B( i ) > 0 );
  }
 const bool routable = source && sink;

 // Reuse a load made by set_Block() or modification processing. A solved
 // instance must be loaded again; apply() also invalidates it on arc changes.
 if( routable ) {
  process_outstanding_Modification();
  if( ! f_instance_fresh )
   load_BenBound();
  f_infeasible = ! sink_reachable();
  }

 if( ! owned )
  f_Block->read_unlock();

 if( f_log )
  setup_seconds = profile_seconds( start );

 if( ! routable )
  return( finish( kError , "no_source_or_sink" ) );

 if( f_log )
 {
   f_log->precision( 17 );
   *f_log << "[compute] solver=" << this << " run=" << f_compute_id
            << " rho=" << f_flow.rate << " burst=" << f_flow.burst
            << " deadline=" << f_flow.deadline << " infeasible="
            << f_infeasible << std::endl;
 }

 if( f_infeasible )    // no path left: nothing to solve
  return( finish( kInfeasible , "unreachable" ) );

 // Screen children against the incumbent before the full Benders solve.
 // This uses only the freshly loaded graph and independent valid bounds;
 // no primal solution or status from a previous node is consulted.
 if( std::isfinite( f_UpCutOff ) &&
     std::any_of( v_closed.begin() , v_closed.end() , []( bool x ) { return x; } ) ) {
  const auto screen_start = f_log ? ProfileClock::now() : ProfileClock::time_point{};
  refine_lb( true );
  if( f_log )
   prepare_seconds += profile_seconds( screen_start );
  if( get_lb() >= f_UpCutOff )
   return( finish( kCutOff , "cutoff_before_benbound" ) );
  }

 const auto ben_start = ProfileClock::now();
 const int rc = SingleFlowDCRBendersSolver::compute( changedvars );
 f_bound_only = false;
 const double ben_elapsed = profile_seconds( ben_start );
 if( f_log )
  benbound_seconds = ben_elapsed;
 if( rc != kOK )
  return( finish( rc , "base_status" ) );

 switch( BenBound::getStat() ) {
  case( BenBound::OK ):
   break;
  case( BenBound::Infeasible ):
   return( finish( kInfeasible , "benbound_infeasible" ) );
  default:
   return( finish( kError , "benbound_error" ) );
  }

 results_available = true;
 lb_before = get_lb();
 if( f_log )
  ub_before = get_ub();

 // Once the node cannot beat the driver's incumbent, neither a branching
 // pair nor closure against the local UB is needed. Retain the actual LB.
 if( std::isfinite( f_UpCutOff ) && ( get_lb() >= f_UpCutOff ) )
  return( finish( kCutOff , "cutoff_after_benbound" ) );

 OFValue ub = f_log ? ub_before : get_ub();
 ub_before = ub;
 OFValue tol = std::max( f_AbsAcc ,
  f_RelAcc * std::max( OFValue( 1 ) , std::abs( ub ) ) );

 // Prepare the branching pair while still in compute(): the fallback can
 // prove the node closed, and the driver must see that bound before it
 // decides whether to call branch(). A finite UB is needed for closure.
 if( ! std::isfinite( ub ) || ( ub - get_lb() > tol ) ) {
  const auto prepare_start = f_log ? ProfileClock::now() : ProfileClock::time_point{};
  prepare_branch( std::min( 1.0 , 0.5 * ben_elapsed ) );
  if( f_log )
   prepare_seconds += profile_seconds( prepare_start );
  }

 ub = get_ub();
 ub_before = ub;
 tol = std::max( f_AbsAcc ,
  f_RelAcc * std::max( OFValue( 1 ) , std::abs( ub ) ) );

 // closed within the tolerance: report LB == UB exactly, so that a driver
 // with zero tolerance does not branch on roundoff
 if( std::isfinite( ub ) && ( ub - get_lb() <= tol ) )
  f_cert = std::max( f_cert , ub );

 if( std::isfinite( f_UpCutOff ) && ( get_lb() >= f_UpCutOff ) )
  return( finish( kCutOff , "cutoff_after_gap_check" ) );

 if( ( ! std::isfinite( ub ) || ( ub - get_lb() > tol ) ) &&
     ( v_branch_plus.empty() || v_branch_minus.empty() ) )
  return( finish( kStopIter , "master_unresolved_no_distinct_paths" ) );

 return( finish( kOK , std::isfinite( ub ) && ( ub - get_lb() <= tol ) ?
                 "gap_closed" : "gap_open" ) );

 }  // end( SFDCRBRelaxationSolver::compute )

/*--------------------------------------------------------------------------*/

Solver::OFValue SFDCRBRelaxationSolver::get_lb( void )
{
 if( f_infeasible || ( ! f_bound_only &&
                       BenBound::getStat() == BenBound::Infeasible ) )
  return( SMSpp_di_unipi_it::Inf< OFValue >() );
 if( f_bound_only || BenBound::getStat() != BenBound::OK )
  return( f_cert );
 return( std::max( SingleFlowDCRBendersSolver::get_lb() , f_cert ) );

 }  // end( SFDCRBRelaxationSolver::get_lb )

/*--------------------------------------------------------------------------*/

void SFDCRBRelaxationSolver::prepare_branch( double refine_seconds )
{
 // A short proof attempt can avoid an entire subtree even when the
 // Lagrangian already offers two distinct paths. Reuse the screening heap.
 refine_lb( false , 24 , refine_seconds );
 const OFValue refined_ub = get_ub();
 const OFValue refined_tol = std::max( f_AbsAcc ,
  f_RelAcc * std::max( OFValue( 1 ) , std::abs( refined_ub ) ) );
 if( ( std::isfinite( f_UpCutOff ) && get_lb() >= f_UpCutOff ) ||
     ( std::isfinite( refined_ub ) && refined_ub - get_lb() <= refined_tol ) )
  return;

 prepare_branch_result();
 invalidate_solution_cache();
 const OFValue prepared_ub = get_ub();
 const OFValue prepared_tol = std::max( f_AbsAcc ,
  f_RelAcc * std::max( OFValue( 1 ) , std::abs( prepared_ub ) ) );
 if( ( std::isfinite( f_UpCutOff ) && get_lb() >= f_UpCutOff ) ||
     ( std::isfinite( prepared_ub ) && prepared_ub - get_lb() <= prepared_tol ) )
  return;
 const auto & result = branch_result();
 if( result.available ) {
  ++f_profile.cache_hits;
  const bool equal = result.plus == result.minus;
  if( f_log )
   *f_log << "[prepare_branch] solver=" << this << " run=" << f_compute_id
          << " source=benbound_master rmin=" << result.rmin
          << " value=" << result.value << " equal=" << equal
          << " local_converged=" << result.local_converged
          << " dual_converged=" << result.dual_converged << '\n';
  if( ! equal ) {
   v_branch_plus = result.plus;
   v_branch_minus = result.minus;
   return;
   }
  }

 // Equality at a fixed rate is not a global certificate if the master
 // stopped early. Retain its actual bound; the fallback can strengthen
 // it, but must not manufacture an unrelated branching pair.
 refine_lb();
 const OFValue ub = get_ub();
 const double tolerance = std::max( f_AbsAcc ,
  f_RelAcc * std::max( OFValue( 1 ) , std::abs( ub ) ) );
 if( ( std::isfinite( f_UpCutOff ) && get_lb() >= f_UpCutOff ) ||
     ( std::isfinite( ub ) && ub - get_lb() <= tolerance ) )
  return;

 // The interval master now selects [f_za,b]. Its bound is evaluated
 // with rates >= f_za, so use that same rate for the original flow.
 // This replaces the old scans over unrelated cached witnesses.
 const auto probe_start = f_log ? ProfileClock::now() : ProfileClock::time_point{};
 evaluate_branch_candidate( f_za , get_lb() );
 invalidate_solution_cache();
 const auto & refined = branch_result();
 const bool distinct = refined.available && refined.plus != refined.minus;
 if( distinct ) {
  v_branch_plus = refined.plus;
  v_branch_minus = refined.minus;
  }
 if( f_log ) {
  ++f_profile.probes;
  f_profile.probe_seconds += profile_seconds( probe_start );
  *f_log << "[prepare_branch] solver=" << this << " run=" << f_compute_id
         << " source=refined_master rmin=" << refined.rmin
         << " value=" << refined.value << " available=" << refined.available
         << " equal=" << ( refined.available && refined.plus == refined.minus )
         << " found=" << distinct << '\n';
  }
}

/*--------------------------------------------------------------------------*/

std::vector< Change * > SFDCRBRelaxationSolver::branch( void )
{
 if( f_log ) {
  f_log->precision( 17 );
  *f_log << "[branch] solver=" << this << " run=" << f_compute_id
         << " stat=" << BenBound::getStat()
         << " r*=" << branch_result().rmin
         << " primal_rmin=" << BenBound::getr_min()
         << " za=" << f_za << " zstar=" << f_zstar
         << " LB=" << get_lb() << " UB=" << get_ub()
         << " heur=" << BenBound::getHeurVal()
         << " approx=" << BenBound::getApproxVal()
         << " objval=" << BenBound::getObjVal()
         << " benIt=" << BenBound::getNumIterationBender()
         << " narcs=" << v_links.size() << std::endl
         << "[true] has=" << has_true_var_solution()
         << " tub=" << get_true_ub() << " tlb=" << get_true_lb()
         << std::endl;
  }

 if( f_infeasible || ( BenBound::getStat() != BenBound::OK ) )
  throw( std::logic_error( "SFDCRBRelaxationSolver::branch: the node has "
                           "no relaxation" ) );

 const auto & Pp = v_branch_plus;
 const auto & Pm = v_branch_minus;
 if( Pp.empty() || Pm.empty() )
  throw( std::logic_error( "SFDCRBRelaxationSolver::branch: no pair of "
                           "distinct paths at the final master candidate" ) );

 // the arcs of the two paths after their common prefix leave the same node
 std::size_t k = 0;
 while( ( k < Pp.size() ) && ( k < Pm.size() ) && ( Pp[ k ] == Pm[ k ] ) )
  ++k;
 if( ( k >= Pp.size() ) || ( k >= Pm.size() ) )
  throw( std::logic_error( "SFDCRBRelaxationSolver::branch: one path is a "
                           "prefix of the other" ) );

 if( f_log )
  *f_log << "[branch] close arc " << Pm[ k ] << " (path not meeting the "
         << "deadline) / arc " << Pp[ k ] << " (path meeting it), common "
         << "prefix " << k << " arcs" << std::endl;

 // the child that keeps the path meeting the deadline comes first
 return( std::vector< Change * >{
          new SingleFlowDCRBlockChange( SingleFlowDCRBlockChange::eCloseArc ,
                                        Index( Pm[ k ] ) ) ,
          new SingleFlowDCRBlockChange( SingleFlowDCRBlockChange::eCloseArc ,
                                        Index( Pp[ k ] ) ) } );

 }  // end( SFDCRBRelaxationSolver::branch )
/*--------------------------------------------------------------------------*/

Change * SFDCRBRelaxationSolver::apply( Change * chg , bool doUndo )
{
 auto dchg = dynamic_cast< SingleFlowDCRBlockChange * >( chg );
 if( ! dchg )  // anything else is for the Block
  return( chg->apply( f_Block , doUndo ) );

 const auto a = dchg->arc();
 if( a >= v_closed.size() )
  throw( std::invalid_argument( "SFDCRBRelaxationSolver::apply: invalid "
                                "arc " + std::to_string( a ) ) );

 Change * undo = nullptr;
 if( doUndo )
  undo = new SingleFlowDCRBlockChange( v_closed[ a ] ?
                                       SingleFlowDCRBlockChange::eCloseArc :
                                       SingleFlowDCRBlockChange::eOpenArc ,
                                       a );

 v_closed[ a ] = ( dchg->type() == SingleFlowDCRBlockChange::eCloseArc );
 f_instance_fresh = false;
 f_bound_only = true;
 f_cert = - SMSpp_di_unipi_it::Inf< OFValue >();
 invalidate_solution_cache();

 return( undo );

 }  // end( SFDCRBRelaxationSolver::apply )

/*--------------------------------------------------------------------------*/

bool SFDCRBRelaxationSolver::arc_closed( Index i )
{
 if( ( i < v_closed.size() ) && v_closed[ i ] )
  return( true );
 auto DCRB = static_cast< SingleFlowDCRBlock * >( f_Block );
 return( DCRB->is_closed( i ) || DCRB->is_deleted( i ) );

 }  // end( SFDCRBRelaxationSolver::arc_closed )

/*--------------------------------------------------------------------------*/

bool SFDCRBRelaxationSolver::sink_reachable( void )
{
 const auto nn = v_nodes.size();
 std::vector< std::vector< Index > > adj( nn );
 for( Index i = 0 ; i < v_links.size() ; ++i )
  if( ! arc_closed( i ) )
   adj[ v_links[ i ].startnode ].push_back( v_links[ i ].endnode );

 std::vector< bool > seen( nn , false );
 std::queue< Index > q;
 q.push( f_flow.sourcenode );
 seen[ f_flow.sourcenode ] = true;
 while( ! q.empty() ) {
  const auto u = q.front();
  q.pop();
  for( auto v : adj[ u ] )
   if( ! seen[ v ] ) {
    seen[ v ] = true;
    q.push( v );
    }
  }

 return( seen[ f_flow.sinknode ] );

 }  // end( SFDCRBRelaxationSolver::sink_reachable )

/*--------------------------------------------------------------------------*/

double SFDCRBRelaxationSolver::interval_lb( double a , double b )
{
 if( f_log )
  ++f_profile.interval_solves;
 DCR::DCRFlow fl = f_flow;
 fl.burst = f_flow.burst * a / b;           // burst' / a == burst / b

 DCRLagrangianSolver pr;
 pr.LoadProblem( int( v_nodes.size() ) , int( v_links.size() ) , fl ,
                 v_links.data() , v_nodes.data() ,
                 static_cast< SingleFlowDCRBlock * >( f_Block )->get_MTU() ,
                 a );
 pr.Solve();

 if( pr.getStatus() == DCRLagrangianSolver::Infeasible )
  return( SMSpp_di_unipi_it::Inf< OFValue >() );     // no solution here
 if( pr.getStatus() != DCRLagrangianSolver::OK )
  return( - SMSpp_di_unipi_it::Inf< OFValue >() );   // claim nothing
 return( pr.getObjVal() );
 }//end (SFDCRBRelaxationSolver::interval_lb )

void SFDCRBRelaxationSolver::refine_lb( bool cutoff_only , int max_splits ,
                                        double max_seconds )
{
 const auto refine_start = ( f_log || max_seconds > 0 ) ?
  ProfileClock::now() : ProfileClock::time_point{};
 const OFValue lb_before = get_lb();
 const char * reason = "iteration_limit";
 int splits = 0;
 if( f_log )
  ++f_profile.refine_calls;
 auto cmp = []( const BoundInterval & x , const BoundInterval & y ) {
  return( x.lb > y.lb );
  };
 auto & Q = v_bound_intervals;
 int interval_solves = 0;
 auto push = [ & ]( double a , double b ) {
  const double lb = interval_lb( a , b );
  ++interval_solves;
  Q.push_back( { a , b , lb } );
  std::push_heap( Q.begin() , Q.end() , cmp );
  };

 if( Q.empty() ) {
  const double L = std::max( f_flow.rate , f_flow.burst / f_flow.deadline );
  double U = L;
  for( const auto & l : v_links )
   U = std::max( U , double( l.capacity ) );
  push( L , U );
  }

 // Before BenBound, only the driver's incumbent is meaningful. After
 // it, the verified local UB can further tighten the stopping target.
 OFValue target = f_UpCutOff;
 if( ! cutoff_only ) {
  const OFValue ub = get_ub();
  if( std::isfinite( ub ) )
   target = std::min( target , ub - std::max( f_AbsAcc ,
    f_RelAcc * std::max( OFValue( 1 ) , std::abs( ub ) ) ) );
  }
 const int remaining = std::max( 0 , f_MaxIv - int( Q.size() - 1 ) );
 const int requested = max_splits >= 0 ? max_splits :
                       ( cutoff_only ? 16 : remaining );
 int budget = std::min( requested , remaining );
 const int extended_budget = cutoff_only ?
  std::min( remaining , budget + 8 ) : budget;
 const double initial_lb = Q.front().lb;
 for( int it = 0 ; it < budget ; ++it ) {
  const BoundInterval top = Q.front();
  if( std::isfinite( target ) && ( top.lb >= target ) ) {
   reason = "target_reached";
   break;
   }
  if( max_seconds > 0 && profile_seconds( refine_start ) >= max_seconds ) {
   reason = "time_budget";
   break;
   }
  if( top.b - top.a <= 1e-1 * f_RelAcc * top.b ) {
   reason = "interval_narrow";
   break;
   }
  const double m = top.a + 0.5 * ( top.b - top.a );
  if( m <= top.a || m >= top.b ) {
   reason = "interval_narrow";
   break;
   }
  std::pop_heap( Q.begin() , Q.end() , cmp );
  Q.pop_back();
  push( top.a , m );
  push( m , top.b );
  ++splits;
  // Give a nearly successful incumbent check eight more splits, only
  // while its certificate is improving. Proximity alone never prunes.
  if( cutoff_only && ( it + 1 == budget ) && budget < extended_budget ) {
   const double lb = Q.front().lb;
   if( std::isfinite( target ) && std::isfinite( lb ) && lb > initial_lb &&
       target - lb > 0 &&
       target - lb <= 1e-4 * std::max( 1.0 , std::abs( target ) ) )
    budget = extended_budget;
   }
  }

 f_cert = std::max( f_cert , Q.front().lb );          // valid: min over leaves
 f_za = Q.front().a;
 f_zstar = 0.5 * ( Q.front().a + Q.front().b );         // where the bound is lowest
 if( f_log ) {
  const double elapsed = profile_seconds( refine_start );
  f_profile.refine_seconds += elapsed;
  f_profile.refine_splits += splits;
  *f_log << "[refine_profile] solver=" << this << " run=" << f_compute_id
         << " elapsed_s=" << elapsed << " splits=" << splits
         << " interval_solves=" << interval_solves << " leaves=" << Q.size()
         << " phase=" << ( cutoff_only ? "before_benbound" :
                              ( max_splits >= 0 ? "before_branch" : "fallback" ) )
         << " max_s=" << max_seconds
         << " reason=" << reason << " target=" << target
         << " LB_before=" << lb_before << " LB=" << get_lb()
         << " a=" << Q.front().a << " b=" << Q.front().b << '\n';
  }
 }//end (SFDCRBRelaxationSolver::refine_lb )

/*--------------------------------------------------------------------------*/
/*------------------- End File SFDCRBRelaxationSolver.cpp ------------------*/
/*--------------------------------------------------------------------------*/
