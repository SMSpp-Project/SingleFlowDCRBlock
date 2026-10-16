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

/// a DCRLagrangianSolver whose state is read at the end of the solve
/** Gives access to the protected state of DCRLagrangianSolver, to find the
 * two paths that define the optimal multiplier [see two_paths()] and to
 * tell whether the multiplier 0 is optimal [see zero_opt()]. */

class SFDCRBProbe : public DCRLagrangianSolver
{
 public:

/*--------------------------------------------------------------------------*/
 /// the two paths that define the optimal multiplier at the loaded r_min
 /** Solves the Lagrangian relaxation and writes in \p Pplus and \p Pminus
  * the paths, as indices of the arcs of the whole graph in the order from
  * the source to the sink, optimal just above and just below the optimal
  * multiplier \f$ \lambda^* \f$, i.e., the one meeting the deadline and
  * the one not meeting it; if the two are the same, the paths at the two
  * cuts of the line search are taken instead. Returns true if the two
  * paths are there and are different. What is found is written on \p log,
  * if any. */

 bool two_paths( std::vector< int > & Pplus , std::vector< int > & Pminus ,
                 std::ostream * log = nullptr ) {
  Solve();

  if( log ){
    log->precision( 12 );
    *log << "[probe] status=" << lagstat << " lambda=" << lambda
    << " pLambda=" << pLambda << " mLambda=" << mLambda
    << " rmin=" << rmin << " d(r*)=" << getObjVal() 
    << " heur=" << getHeurVal() 
    << " pCut(q,m)=(" << pCut.q << "," << pCut.m << ")"
    << " mCut(q,m)=(" << mCut.q << "," << mCut.m << ")" 
    << std::endl;
  }

  if( ( lagstat != DCRLagrangianSolver::OK ) || ( lambda <= 0 ) )
   return( false );

  const double ls = lambda;  // path_at() overwrites lambda
  const double d = 1e-6 * ls;

  Pminus = path_at( ls - d );
  Pplus = path_at( ls + d );

  if( Pminus == Pplus ) {    // the paths at the two cuts
   Pminus = path_at( pLambda );
   Pplus = path_at( mLambda );
   }

  if( log )
   *log << "[probe] |P-|=" << Pminus.size() << " |P+|=" << Pplus.size()
        << " equal=" << ( Pminus == Pplus ) << std::endl;

  return( ( ! Pminus.empty() ) && ( ! Pplus.empty() ) &&
          ( Pminus != Pplus ) );
  }

/*--------------------------------------------------------------------------*/
 /// true if the multiplier 0 is optimal and its path meets the deadline

 bool zero_opt( void ) {
  if( isFeasible() )  // true if the relaxation is infeasible
   return( false );
  // Inf is that of DCRLagrangianSolver, by which getHeurVal() says that
  // there is no feasible path
  return( ( getLambda() == 0 ) && ( getHeurVal() < Inf< double >() ) );
  }

/*--------------------------------------------------------------------------*/

 private:

 /// the shortest path at multiplier lam, as indices of the arcs of the
 /// whole graph, empty if the shortest path tree fails

 std::vector< int > path_at( double lam ) {
  lambda = lam;
  setSPTcosts();  // the costs the cuts are generated with
  spt.updCosts( Linksp );
  spt.Solve();
  if( spt.getStatus() != SPT::OK )
   return( std::vector< int >() );

  auto p = spt.getPath();
  std::vector< int > out( p.size() );
  for( std::size_t k = 0 ; k < p.size() ; ++k )
   out[ k ] = RedGraPos[ p[ k ] ];  // reduced graph -> whole graph
  return( out );
  }

 };  // end( class( SFDCRBProbe ) )

}  // end( anonymous namespace )

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS OF SFDCRBRelaxationSolver -----------------*/
/*--------------------------------------------------------------------------*/

void SFDCRBRelaxationSolver::set_Block( Block * block )
{
 if( block && ( ! dynamic_cast< SingleFlowDCRBlock * >( block ) ) )
  throw( std::invalid_argument( "SFDCRBRelaxationSolver::set_Block: block "
                                "must be a SingleFlowDCRBlock" ) );

 v_closed.clear();  // before the base loads the instance [see arc_closed()]
 f_infeasible = false;
 f_cert = - SMSpp_di_unipi_it::Inf< OFValue >();
 f_zstar = 0;
 f_za = 0;

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
 return( SingleFlowDCRBendersSolver::get_dbl_par( par ) );

 }  // end( SFDCRBRelaxationSolver::get_dbl_par )

/*--------------------------------------------------------------------------*/

int SFDCRBRelaxationSolver::compute( bool changedvars )
{
 if( ! f_Block )
  return( kError );

 f_cert = - SMSpp_di_unipi_it::Inf< OFValue >();
 f_zstar = 0;
 f_za = 0;

 const bool owned = f_Block->is_owned_by( f_id );
 if( ( ! owned ) && ( ! f_Block->read_lock() ) )
  return( kBlockLocked );

 // a flow needs a source and a sink [see load_BenBound()]
 auto DCRB = static_cast< SingleFlowDCRBlock * >( f_Block );
 bool source = false;
 bool sink = false;
 for( Index i = 0 ; i < DCRB->get_NNodes() ; ++i ) {
  source |= ( DCRB->get_B( i ) < 0 );
  sink |= ( DCRB->get_B( i ) > 0 );
  }
 const bool routable = source && sink;

 // BenBound cannot solve the same data twice, so the instance, with the
 // arcs closed by the Solver, is loaded at every compute(), after the
 // Modification of the Block
 if( routable ) {
  process_outstanding_Modification();
  load_BenBound();
  f_infeasible = ! sink_reachable();
  }

 if( ! owned )
  f_Block->read_unlock();

 if( ! routable )
  return( kError );

 if( f_log )
  f_log->precision( 17 );
  *f_log << "[compute] rho=" << f_flow.rate << " burst=" << f_flow.burst
         << " deadline=" << f_flow.deadline << " infeasible="
         << f_infeasible << std::endl;

 if( f_infeasible )    // no path left: nothing to solve
  return( kInfeasible );

 const int rc = SingleFlowDCRBendersSolver::compute( changedvars );
 if( rc != kOK )
  return( rc );

 switch( BenBound::getStat() ) {
  case( BenBound::OK ):
   break;
  case( BenBound::Infeasible ):
   return( kInfeasible );
  default:
   return( kError );
  }

 const OFValue ub = get_ub();
 if( ub - SingleFlowDCRBendersSolver::get_lb() >
     f_RelAcc * std::max( OFValue( 1 ) , std::abs( ub ) ) )
  refine_lb();
 // closed within the tolerance: report LB == UB exactly, so that a driver
 // with zero tolerance does not branch on roundoff
 if (ub - get_lb() <= f_RelAcc * std::max(OFValue(1), std::abs(ub)))
   f_cert = std::max(f_cert, ub);

 return( kOK );

 }  // end( SFDCRBRelaxationSolver::compute )

/*--------------------------------------------------------------------------*/

Solver::OFValue SFDCRBRelaxationSolver::get_lb( void )
{
 if( f_infeasible || ( BenBound::getStat() == BenBound::Infeasible ) )
  return( SMSpp_di_unipi_it::Inf< OFValue >() );
 if( BenBound::getStat() != BenBound::OK )
  return( - SMSpp_di_unipi_it::Inf< OFValue >() );
 return( std::max( SingleFlowDCRBendersSolver::get_lb() , f_cert ) );

 }  // end( SFDCRBRelaxationSolver::get_lb )

/*--------------------------------------------------------------------------*/

std::vector< Change * > SFDCRBRelaxationSolver::branch( void )
{
 if( f_log ) {
  f_log->precision( 17 );
  *f_log << "[branch] stat=" << BenBound::getStat()
         << " r*=" << BenBound::getr_min()
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

 // the values of r_min at which to look for the two paths, in order: the
 // left end of the interval that gives the bound (it is in the domain by
 // construction), the point BenBound found, the middle of the interval
 // (which can fall just past the end of the domain)
 std::vector< double > cands;
 if( f_za > 0 )
  cands.push_back( f_za );
 cands.push_back( BenBound::getr_min() );
 if( f_zstar > 0 )
  cands.push_back( f_zstar );

 const double mtu =
  static_cast< SingleFlowDCRBlock * >( f_Block )->get_MTU();

 std::vector< int > Pp , Pm;
 bool found = false;
 for( const double z : cands ) {
  SFDCRBProbe pr;
  pr.LoadProblem( int( v_nodes.size() ) , int( v_links.size() ) , f_flow ,
                  v_links.data() , v_nodes.data() , mtu , z );

  Pp.clear();
  Pm.clear();
  if( pr.two_paths( Pp , Pm , f_log ) ) {
   found = true;
   break;
   }
  }

 if( ! found )
  throw( std::logic_error( "SFDCRBRelaxationSolver::branch: no pair of "
                           "distinct paths at any candidate r_min" ) );

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
 DCR::DCRFlow fl = f_flow;
 fl.burst = f_flow.burst * a / b;           // burst' / a == burst / b

 SFDCRBProbe pr;
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

void SFDCRBRelaxationSolver::refine_lb( void )
{
 struct Iv { double a , b , lb; };
 auto cmp = []( const Iv & x , const Iv & y ) { return( x.lb > y.lb ); };
 std::priority_queue< Iv , std::vector< Iv > , decltype( cmp ) > Q( cmp );

 const double L = std::max( f_flow.rate , f_flow.burst / f_flow.deadline );
 double U = L;
 for( auto & l : v_links )
  U = std::max( U , double( l.capacity ) );

 Q.push( { L , U , interval_lb( L , U ) } );
 const OFValue ub = get_ub();

 for( int it = 0 ; it < f_MaxIv ; ++it ) {
  const Iv top = Q.top();
  if (top.lb >= ub - f_RelAcc * std::max(1.0, std::abs(ub)) )
    break; // node closed
  if( top.b - top.a <= 1e-1*f_RelAcc * top.b )
   break;                                           // interval narrow enough
  Q.pop();
  const double m = 0.5 * ( top.a + top.b );
  Q.push( { top.a , m , interval_lb( top.a , m ) } );
  Q.push( { m , top.b , interval_lb( m , top.b ) } );
  }

 f_cert = std::max( f_cert , Q.top().lb );          // valid: min over leaves
 f_za = Q.top().a;
 f_zstar = 0.5 * ( Q.top().a + Q.top().b );         // where the bound is lowest
 }//end (SFDCRBRelaxationSolver::refine_lb )

/* void SFDCRBRelaxationSolver::certify( void )
{
 const double rlow = std::max( f_flow.rate ,
                               f_flow.burst / f_flow.deadline );
 SFDCRBProbe pr;
 pr.LoadProblem( int( v_nodes.size() ) , int( v_links.size() ) , f_flow ,
                 v_links.data() , v_nodes.data() ,
                 static_cast< SingleFlowDCRBlock * >( f_Block )->get_MTU() ,
                 rlow );
 if( pr.zero_opt() )
  f_cert = pr.getObjVal();

 if( f_log )
  *f_log << "[cert] rlow=" << rlow << " f_cert=" << f_cert << std::endl;

 }  // end( SFDCRBRelaxationSolver::certify ) */

/*--------------------------------------------------------------------------*/
/*------------------- End File SFDCRBRelaxationSolver.cpp ------------------*/
/*--------------------------------------------------------------------------*/
