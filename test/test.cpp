/*--------------------------------------------------------------------------*/
/*-------------------------- File test.cpp ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Tester of the SingleFlowDCRBlock module that needs nothing but the core
 * SMS++ library: all the instances are tiny graphs built in memory, whose
 * s-t paths are few enough to be enumerated, so that every answer of the
 * solvers of the module is compared with the one the enumeration gives.
 *
 * - SPT, on random graphs with non-negative costs and on random acyclic
 *   ones with negative costs, is compared with the shortest of the simple
 *   paths, both in the path it returns and in the label of every node.
 *
 * - DCR_SPT, with both its heuristics, is compared with the paths that meet
 *   the deadline: ERA-I says an instance is infeasible exactly when no path
 *   does with the whole capacity reserved on every arc, and whatever either
 *   heuristic returns is a path meeting the deadline, whose cost is not
 *   below the continuous optimum.
 *
 * - DCRLagrangianSolver, for a given r_min, has to give a Lagrangian value
 *   not above the optimum of the problem it relaxes and a heuristic value
 *   not below the optimum of the DCR problem on the same reduced graph.
 *
 * - SingleFlowDCRBlock is checked in the feasibility of given solutions,
 *   the default data of load(), the netCDF round trip, the feasibility of
 *   the instance against the enumeration, and the Modification its changes
 *   issue, which a Solver attached to it sees when it solves again.
 *
 * The edge cases are checked on their own: source equal to sink, no path at
 * all, every path violating the deadline, a deadline met with equality, ties
 * between paths, zero-cost arcs, parallel arcs and closed arcs.
 *
 * The exit code is the number of failed checks.
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "SingleFlowDCRBlock.h"

#include "DCRLagrangianSolver.h"

#include "DCR_SPT.h"

#include "SPT.h"

#include "Solver.h"

#include <cmath>

#include <cstdio>

#include <functional>

#include <iostream>

#include <limits>

#include <random>

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*-------------------------------- GLOBALS ---------------------------------*/
/*--------------------------------------------------------------------------*/

static int failures = 0;  // number of failed checks

static const double INF = std::numeric_limits< double >::max();

using Path = std::vector< int >;  // arc indices, from source to sink

/*--------------------------------------------------------------------------*/
/*------------------------------ FUNCTIONS ---------------------------------*/
/*--------------------------------------------------------------------------*/

static void check( bool ok , const std::string & what )
{
 if( ! ok ) {
  ++failures;
  std::cout << "FAILED: " << what << std::endl;
  }
 }

/*--------------------------------------------------------------------------*/
/// true if a and b agree up to a relative tolerance

static bool near( double a , double b , double tol = 1e-9 )
{
 if( a == b )
  return( true );
 return( std::abs( a - b ) <= tol * std::max( { 1.0 , std::abs( a ) ,
                                                  std::abs( b ) } ) );
 }

/*--------------------------------------------------------------------------*/
/*------------------------------ THE NETWORK -------------------------------*/
/*--------------------------------------------------------------------------*/
/* A network with 0-based nodes; the speed of an arc is its capacity, as in
 * SingleFlowDCRBlock::is_feasible_instance(), so that the fixed delay of an
 * arc is MTU / cap + ldel + ndel of its tail and the one depending on the
 * reserved rate r is MTU / r, plus burst / min r for the path as a whole. */

struct Net
{
 int n = 0;
 std::vector< int > sn , en;
 std::vector< double > cap , cost , ldel;
 std::vector< double > ndel;
 double MTU = 1;
 double burst = 0;
 double deadline = 0;
 double rate = 0;
 int s = 0 , t = 0;

 int m( void ) const { return( int( sn.size() ) ); }

 void arc( int i , int j , double u , double c , double ld = 0 ) {
  sn.push_back( i ); en.push_back( j );
  cap.push_back( u ); cost.push_back( c ); ldel.push_back( ld );
  }
 };

/*--------------------------------------------------------------------------*/
/// all the simple paths from u to t whose arcs pass allowed

static void paths_dfs( const Net & g , int u , int t ,
                       const std::function< bool( int ) > & allowed ,
                       std::vector< bool > & vis , Path & cur ,
                       std::vector< Path > & out )
{
 if( u == t ) {
  out.push_back( cur );
  return;
  }
 for( int a = 0 ; a < g.m() ; ++a )
  if( ( g.sn[ a ] == u ) && ( ! vis[ g.en[ a ] ] ) && allowed( a ) ) {
   vis[ g.en[ a ] ] = true;
   cur.push_back( a );
   paths_dfs( g , g.en[ a ] , t , allowed , vis , cur , out );
   cur.pop_back();
   vis[ g.en[ a ] ] = false;
   }
 }

static std::vector< Path > all_paths( const Net & g , int s , int t ,
                                      const std::function< bool( int ) > &
                                      allowed = []( int ) { return true; } )
{
 std::vector< Path > out;
 std::vector< bool > vis( g.n , false );
 vis[ s ] = true;
 Path cur;
 paths_dfs( g , s , t , allowed , vis , cur , out );
 return( out );
 }

/*--------------------------------------------------------------------------*/
/// true if P is a path from s to t in g, its arcs in source-to-sink order

static bool is_path( const Net & g , const Path & P , int s , int t )
{
 int u = s;
 for( int a : P ) {
  if( ( a < 0 ) || ( a >= g.m() ) || ( g.sn[ a ] != u ) )
   return( false );
  u = g.en[ a ];
  }
 return( u == t );
 }

/*--------------------------------------------------------------------------*/
/// the delay of P, with R[ k ] the rate reserved on its k-th arc

static double path_delay( const Net & g , const Path & P ,
                          const std::vector< double > & R )
{
 double d = 0;
 double rmin = INF;
 for( size_t k = 0 ; k < P.size() ; ++k ) {
  const int a = P[ k ];
  d = g.MTU / R[ k ] + g.MTU / g.cap[ a ] + g.ndel[ g.sn[ a ] ] +
      g.ldel[ a ] + d;
  rmin = std::min( rmin , R[ k ] );
  }
 return( P.empty() ? d : d + g.burst / rmin );
 }

/*--------------------------------------------------------------------------*/
/// the minimum capacity on P (+INF for the empty path)

static double min_cap( const Net & g , const Path & P )
{
 double u = INF;
 for( int a : P )
  u = std::min( u , g.cap[ a ] );
 return( u );
 }

/*--------------------------------------------------------------------------*/
/// true if P meets the deadline with the whole capacity reserved

static bool full_feasible( const Net & g , const Path & P )
{
 if( min_cap( g , P ) < g.rate )
  return( false );
 std::vector< double > R;
 for( int a : P )
  R.push_back( g.cap[ a ] );
 return( path_delay( g , P , R ) <= g.deadline );
 }

/*--------------------------------------------------------------------------*/
/* The cheapest rates on P, every one in [ m , cap ], with the burst term
 * taken as burst / m: min sum c r s.t. sum MTU / r <= slack, which is
 * separable and convex, hence solved by bisection on the multiplier of the
 * delay constraint. +INF if not even the whole capacity meets the deadline.
 * With rates of 1 on the empty path, since it has no arc. */

static double path_opt_lb( const Net & g , const Path & P , double m )
{
 if( P.empty() )
  return( g.burst / m <= g.deadline ? 0 : INF );

 double fixed = 0;
 for( int a : P ) {
  if( g.cap[ a ] < m )
   return( INF );
  fixed += g.MTU / g.cap[ a ] + g.ndel[ g.sn[ a ] ] + g.ldel[ a ];
  }
 const double slack = g.deadline - g.burst / m - fixed;

 auto eval = [ & ]( double mu , double & cst ) {
  double dl = 0;
  cst = 0;
  for( int a : P ) {
   double r = g.cap[ a ];
   if( g.cost[ a ] > 0 )
    r = std::min( g.cap[ a ] ,
                  std::max( m , std::sqrt( mu * g.MTU / g.cost[ a ] ) ) );
   dl += g.MTU / r;
   cst += g.cost[ a ] * r;
   }
  return( dl );
  };

 double cst;
 const double dmin = eval( INF , cst );   // every rate at its capacity
 if( dmin > slack * ( 1 + 1e-12 ) )
  return( INF );
 if( dmin >= slack )                      // no slack left
  return( cst );
 if( eval( 0 , cst ) <= slack )
  return( cst );

 double lo = 0 , hi = 1;
 while( eval( hi , cst ) > slack )
  hi *= 2;
 for( int it = 0 ; it < 200 ; ++it ) {
  const double mid = ( lo + hi ) / 2;
  if( eval( mid , cst ) > slack )
   lo = mid;
  else
   hi = mid;
  }
 eval( hi , cst );
 return( cst );
 }

/*--------------------------------------------------------------------------*/
/* The cheapest rates on P, every one >= rmin, with the burst term taken as
 * burst / min r: the minimum over m in [ rmin , min cap ] of path_opt_lb(),
 * a convex function of m, which is finite on a right interval. */

static double path_opt_true( const Net & g , const Path & P , double rmin )
{
 if( P.empty() )
  return( path_opt_lb( g , P , rmin ) );
 double lo = rmin , hi = min_cap( g , P );
 if( hi < lo )
  return( INF );
 for( int it = 0 ; it < 100 ; ++it ) {
  const double m1 = lo + ( hi - lo ) / 3;
  const double m2 = hi - ( hi - lo ) / 3;
  // infeasible values lie on the left, a larger m being never worse
  const double g1 = path_opt_lb( g , P , m1 );
  if( ( g1 == INF ) || ( g1 > path_opt_lb( g , P , m2 ) ) )
   lo = m1;
  else
   hi = m2;
  }
 return( std::min( { path_opt_lb( g , P , lo ) , path_opt_lb( g , P , hi ) ,
                     path_opt_lb( g , P , min_cap( g , P ) ) } ) );
 }

/*--------------------------------------------------------------------------*/
/*------------------------ THE DATA OF THE DCR SOLVERS ---------------------*/
/*--------------------------------------------------------------------------*/

struct DCRData
{
 std::vector< DCR::DCRLink > links;
 std::vector< DCR::DCRNode > nodes;
 std::vector< double > costs , caps;
 DCR::DCRFlow flow{};
 };

static void fill( DCRData & d , const Net & g )
{
 d.links.assign( g.m() , DCR::DCRLink{} );
 for( int a = 0 ; a < g.m() ; ++a ) {
  d.links[ a ].startnode = g.sn[ a ];
  d.links[ a ].endnode = g.en[ a ];
  d.links[ a ].speed = g.cap[ a ];
  d.links[ a ].capacity = g.cap[ a ];
  d.links[ a ].delay = g.ldel[ a ];
  d.links[ a ].cost = g.cost[ a ];
  }
 d.nodes.assign( g.n , DCR::DCRNode{} );
 for( int i = 0 ; i < g.n ; ++i )
  d.nodes[ i ].delay = g.ndel[ i ];
 d.costs = g.cost;
 d.caps = g.cap;
 d.flow.sourcenode = g.s;
 d.flow.sinknode = g.t;
 d.flow.burst = g.burst;
 d.flow.rate = g.rate;
 d.flow.deadline = g.deadline;
 d.flow.costs = d.costs.data();
 d.flow.caps = d.caps.data();
 }

/*--------------------------------------------------------------------------*/
/// a random network: n nodes, up to maxm arcs, parallel ones included

static Net random_net( std::mt19937 & gen , int maxn , int maxm )
{
 std::uniform_int_distribution< int > dn( 2 , maxn );
 Net g;
 g.n = dn( gen );
 std::uniform_int_distribution< int > dm( 0 , maxm );
 std::uniform_int_distribution< int > dv( 0 , g.n - 1 );
 std::uniform_int_distribution< int > dc( 0 , 5 );
 std::uniform_int_distribution< int > du( 0 , 3 );
 std::uniform_int_distribution< int > dl( 0 , 2 );
 const int m = dm( gen );
 for( int k = 0 ; k < m ; ++k ) {
  const int i = dv( gen ) , j = dv( gen );
  if( i == j )
   continue;
  g.arc( i , j , double( 1 << du( gen ) ) , dc( gen ) , 0.25 * dl( gen ) );
  }
 g.ndel.assign( g.n , 0 );
 for( int i = 0 ; i < g.n ; ++i )
  g.ndel[ i ] = 0.125 * dl( gen );
 g.s = 0;
 g.t = g.n - 1;
 std::uniform_real_distribution< double > dd( 0.5 , 6 );
 g.deadline = dd( gen );
 g.burst = dl( gen );
 g.rate = du( gen ) ? 0 : 2;
 return( g );
 }

/*--------------------------------------------------------------------------*/
/*---------------------------------- SPT -----------------------------------*/
/*--------------------------------------------------------------------------*/
/// solves g with SPT, costs g.cost, and compares it with the enumeration

static void check_spt( const Net & g , const std::string & what )
{
 std::vector< SPT::SPTLink > links( g.m() );
 for( int a = 0 ; a < g.m() ; ++a )
  links[ a ] = { g.sn[ a ] , g.en[ a ] , g.cost[ a ] , 0 };

 SPT spt;
 spt.LoadProblem( g.n , g.m() , links , g.s , g.t );
 spt.Solve();

 // the shortest simple path to every node
 std::vector< double > best( g.n , INF );
 best[ g.s ] = 0;
 for( int v = 0 ; v < g.n ; ++v )
  if( v != g.s )
   for( const auto & P : all_paths( g , g.s , v ) ) {
    double c = 0;
    for( int a : P )
     c += g.cost[ a ];
    best[ v ] = std::min( best[ v ] , c );
    }

 if( best[ g.t ] == INF ) {
  check( spt.getStatus() == SPT::Error , what + ": no path, Error expected" );
  check( spt.getPath().empty() && ( spt.getNHops() == 0 ) ,
         what + ": no path, yet a path is returned" );
  return;
  }

 check( spt.getStatus() == SPT::OK , what + ": status not OK" );
 const auto P = spt.getPath();
 check( is_path( g , P , g.s , g.t ) , what + ": not an s-t path" );
 check( int( P.size() ) == spt.getNHops() , what + ": getNHops() differs" );
 double c = 0;
 for( int a : P )
  if( ( a >= 0 ) && ( a < g.m() ) )
   c += g.cost[ a ];
 check( near( c , best[ g.t ] ) , what + ": path cost " + std::to_string( c )
        + " instead of " + std::to_string( best[ g.t ] ) );

 const auto L = spt.getLabel();
 check( int( L.size() ) == g.n , what + ": wrong number of labels" );
 for( int v = 0 ; v < std::min( g.n , int( L.size() ) ) ; ++v )
  check( near( L[ v ] , best[ v ] ) , what + ": label of node " +
         std::to_string( v ) + " is " + std::to_string( L[ v ] ) +
         " instead of " + std::to_string( best[ v ] ) );
 }

/*--------------------------------------------------------------------------*/
/// random graphs with non-negative costs (cycles, parallel arcs, zero costs)
/// and random acyclic ones with negative costs

static void test_spt_random( void )
{
 std::mt19937 gen( 20260926 );
 for( int k = 0 ; k < 300 ; ++k ) {
  Net g = random_net( gen , 6 , 12 );
  check_spt( g , "SPT random #" + std::to_string( k ) );
  }

 std::uniform_int_distribution< int > dc( -5 , 5 );
 for( int k = 0 ; k < 200 ; ++k ) {
  Net g = random_net( gen , 6 , 12 );
  for( int a = 0 ; a < g.m() ; ++a ) {   // arcs from lower to higher nodes
   if( g.sn[ a ] > g.en[ a ] )
    std::swap( g.sn[ a ] , g.en[ a ] );
   g.cost[ a ] = dc( gen );
   }
  check_spt( g , "SPT acyclic, negative costs #" + std::to_string( k ) );
  }
 }

/*--------------------------------------------------------------------------*/
/// the edge cases of SPT: source equal to sink, no path, parallel arcs,
/// ties, a zero-cost cycle, and a change of the costs

static void test_spt_edges( void )
{
 // source equal to sink: the empty path
 {
  Net g;
  g.n = 3;
  g.arc( 0 , 1 , 1 , 2 ); g.arc( 1 , 2 , 1 , 3 );
  g.ndel.assign( 3 , 0 );
  g.s = g.t = 1;
  check_spt( g , "SPT source == sink" );
  }

 // no path at all, not even an arc
 {
  Net g;
  g.n = 2;
  g.ndel.assign( 2 , 0 );
  g.s = 0; g.t = 1;
  check_spt( g , "SPT no arc" );
  g.arc( 1 , 0 , 1 , 1 );  // the wrong way
  check_spt( g , "SPT arc the wrong way" );
  }

 // parallel arcs: the cheaper one, whatever its position
 {
  Net g;
  g.n = 2;
  g.arc( 0 , 1 , 1 , 5 ); g.arc( 0 , 1 , 1 , 2 ); g.arc( 0 , 1 , 1 , 7 );
  g.ndel.assign( 2 , 0 );
  g.s = 0; g.t = 1;
  check_spt( g , "SPT parallel arcs" );
  std::vector< SPT::SPTLink > links( 3 );
  for( int a = 0 ; a < 3 ; ++a )
   links[ a ] = { 0 , 1 , g.cost[ a ] , 0 };
  SPT spt;
  spt.LoadProblem( 2 , 3 , links , 0 , 1 );
  spt.Solve();
  check( spt.getPath() == Path{ 1 } , "SPT parallel arcs: not arc 1" );
  }

 // ties and a zero-cost cycle 1 -> 2 -> 1, which must not loop forever
 {
  Net g;
  g.n = 4;
  g.arc( 0 , 1 , 1 , 1 ); g.arc( 0 , 2 , 1 , 1 );
  g.arc( 1 , 3 , 1 , 2 ); g.arc( 2 , 3 , 1 , 2 );
  g.arc( 1 , 2 , 1 , 0 ); g.arc( 2 , 1 , 1 , 0 );
  g.ndel.assign( 4 , 0 );
  g.s = 0; g.t = 3;
  check_spt( g , "SPT ties and zero-cost cycle" );
  }

 // a change of the costs, then of the whole graph after a failure
 {
  std::vector< SPT::SPTLink > links = {
   { 0 , 1 , 1 , 0 } , { 1 , 2 , 1 , 0 } , { 0 , 2 , 3 , 0 } };
  SPT spt;
  spt.LoadProblem( 3 , 3 , links , 0 , 2 );
  spt.Solve();
  check( spt.getPath() == Path( { 0 , 1 } ) , "SPT before updCosts()" );
  links[ 2 ].cost = 1;
  spt.updCosts( links );
  spt.Solve();
  check( spt.getPath() == Path{ 2 } , "SPT after updCosts()" );
  check( near( spt.getLabel()[ 2 ] , 1 ) , "SPT label after updCosts()" );

  std::vector< SPT::SPTLink > back = { { 2 , 0 , 1 , 0 } };
  spt.LoadProblem( 3 , 1 , back , 0 , 2 );
  spt.Solve();
  check( spt.getStatus() == SPT::Error , "SPT reloaded with no path" );
  spt.LoadProblem( 3 , 3 , links , 0 , 2 );
  spt.Solve();
  check( ( spt.getStatus() == SPT::OK ) && ( spt.getPath() == Path{ 2 } ) ,
         "SPT reloaded after a failure" );
  }
 }

/*--------------------------------------------------------------------------*/
/*-------------------------------- DCR_SPT ---------------------------------*/
/*--------------------------------------------------------------------------*/
/* Solves g with one heuristic of DCR_SPT and checks what it returns: the
 * status of ERA-I against the enumeration, and for both a path meeting the
 * deadline at the rates returned, whose cost is the objective and is not
 * below the continuous optimum opt; returns the status. */

static DCR::DCRStatus check_dcr_spt( const Net & g , char heur , double opt ,
                                     const std::string & what )
{
 DCRData d;
 fill( d , g );
 DCR_SPT slv;
 slv.DCRloadProblem( g.n , g.m() , 1 , &d.flow , d.links.data() ,
                     d.nodes.data() , g.MTU , DCR::SRP );
 slv.DCRsetHeur( heur );
 const auto st = slv.DCRsolve();
 const std::string who = what + ", ERA-" + ( heur == '1' ? "I" : "H" );

 bool feasible = false;
 double cheapest = INF;   // cheapest path meeting the deadline at full rate
 for( const auto & P : all_paths( g , g.s , g.t ) )
  if( full_feasible( g , P ) ) {
   feasible = true;
   double c = 0;
   for( int a : P )
    c += g.cost[ a ] * g.cap[ a ];
   cheapest = std::min( cheapest , c );
   }

 if( heur == '1' )
  check( ( st == DCR::OK ) == feasible , who + ": status " +
         std::to_string( st ) + ", while a path meeting the deadline " +
         ( feasible ? "exists" : "does not exist" ) );
 else
  check( ( st != DCR::OK ) || feasible , who + ": OK with no feasible path" );

 if( st != DCR::OK )
  return( st );

 check( slv.DCRgetUBSolNPaths( 0 ) == 1 , who + ": not a single path" );
 std::vector< int > X( std::max( g.n - 1 , 1 ) );
 std::vector< double > R( std::max( g.n - 1 , 1 ) );
 const int nh = slv.DCRgetUBSolPath( 0 , 0 , X.data() , R.data() );
 // the arcs come from the sink back to the source
 Path P( X.rbegin() + ( X.size() - nh ) , X.rend() );
 std::vector< double > RP( R.rbegin() + ( R.size() - nh ) , R.rend() );
 check( is_path( g , P , g.s , g.t ) , who + ": not an s-t path" );
 if( ! is_path( g , P , g.s , g.t ) )
  return( st );

 double c = 0;
 for( size_t k = 0 ; k < P.size() ; ++k ) {
  const int a = P[ k ];
  c += g.cost[ a ] * RP[ k ];
  check( ( RP[ k ] >= g.rate ) && ( RP[ k ] <= g.cap[ a ] * ( 1 + 1e-12 ) ) ,
         who + ": rate " + std::to_string( RP[ k ] ) + " out of [ rate , "
         "cap ]" );
  if( heur == '1' )
   check( RP[ k ] == g.cap[ a ] , who + ": ERA-I rate is not the capacity" );
  else
   check( RP[ k ] == RP[ 0 ] , who + ": ERA-H rates are not all equal" );
  }
 check( path_delay( g , P , RP ) <= g.deadline * ( 1 + 1e-9 ) ,
        who + ": delay " + std::to_string( path_delay( g , P , RP ) ) +
        " over the deadline " + std::to_string( g.deadline ) );
 check( near( c , slv.DCRgetObj() ) , who + ": objective " +
        std::to_string( slv.DCRgetObj() ) + " is not the cost " +
        std::to_string( c ) + " of the path" );
 check( slv.DCRgetObj() >= opt * ( 1 - 1e-7 ) - 1e-9 , who + ": objective " +
        std::to_string( slv.DCRgetObj() ) + " below the optimum " +
        std::to_string( opt ) );
 if( heur == '1' )
  check( slv.DCRgetObj() >= cheapest * ( 1 - 1e-12 ) , who + ": objective "
         "below the cheapest path at full rate" );

 bool threw = false;
 try { slv.DCRgetPSolNPaths( 0 ); }
 catch( DCR::DCRException & ) { threw = true; }
 check( threw , who + ": DCRgetPSolNPaths() does not throw" );

 return( st );
 }

/*--------------------------------------------------------------------------*/
/// the continuous optimum of g over all its paths, rates >= g.rate

static double dcr_optimum( const Net & g )
{
 double opt = INF;
 for( const auto & P : all_paths( g , g.s , g.t ) )
  opt = std::min( opt , path_opt_true( g , P , std::max( g.rate , 1e-12 ) ) );
 return( opt );
 }

/*--------------------------------------------------------------------------*/
/// random instances, both heuristics

static void test_dcr_spt_random( void )
{
 std::mt19937 gen( 42 );
 int nok[ 2 ] = { 0 , 0 };
 for( int k = 0 ; k < 300 ; ++k ) {
  const Net g = random_net( gen , 5 , 9 );
  const double opt = dcr_optimum( g );
  const std::string what = "DCR_SPT random #" + std::to_string( k );
  const auto s1 = check_dcr_spt( g , '1' , opt , what );
  const auto s2 = check_dcr_spt( g , '2' , opt , what );
  nok[ 0 ] += ( s1 == DCR::OK );
  nok[ 1 ] += ( s2 == DCR::OK );
  check( ( s2 != DCR::OK ) || ( s1 == DCR::OK ) , what + ": ERA-H finds a "
         "path, ERA-I does not" );
  }
 // the instances are not all trivially infeasible, nor all feasible
 check( ( nok[ 0 ] > 50 ) && ( nok[ 0 ] < 250 ) && ( nok[ 1 ] > 20 ) ,
        "DCR_SPT random: " + std::to_string( nok[ 0 ] ) + " and " +
        std::to_string( nok[ 1 ] ) + " feasible out of 300" );
 }

/*--------------------------------------------------------------------------*/
/* The edge cases of DCR_SPT. The instance of reference is the path 0 -> 1
 * -> 2 with capacities 2 and MTU 1, so that each arc has a fixed delay of
 * 1 / 2 and one of 1 / 2 at full rate, and burst 1, which adds 1 / 2 at
 * full rate: 2.5 in all, a number exactly representable. */

static Net base_path( void )
{
 Net g;
 g.n = 3;
 g.arc( 0 , 1 , 2 , 1 ); g.arc( 1 , 2 , 2 , 3 );
 g.ndel.assign( 3 , 0 );
 g.s = 0; g.t = 2;
 g.burst = 1;
 g.rate = 0;
 g.deadline = 10;
 return( g );
 }

static void test_dcr_spt_edges( void )
{
 const char heurs[] = { '1' , '2' };

 // source equal to sink: the empty path, of cost 0
 {
  Net g = base_path();
  g.s = g.t = 1;
  for( char h : heurs ) {
   const auto st = check_dcr_spt( g , h , 0 , "DCR_SPT source == sink" );
   check( st == DCR::OK , "DCR_SPT source == sink: not OK" );
   }
  }

 // no path at all
 {
  Net g = base_path();
  g.sn[ 1 ] = 2; g.en[ 1 ] = 1;  // 2 -> 1 instead of 1 -> 2
  for( char h : heurs )
   check( check_dcr_spt( g , h , INF , "DCR_SPT no path" ) ==
          DCR::Infeasible , "DCR_SPT no path: not Infeasible" );
  }

 // the only path violates the deadline, even at full rate
 {
  Net g = base_path();
  g.deadline = 2.4;
  for( char h : heurs )
   check( check_dcr_spt( g , h , INF , "DCR_SPT deadline violated" ) ==
          DCR::Infeasible , "DCR_SPT deadline violated: not Infeasible" );
  }

 // the deadline is met with equality at full rate: the path is feasible,
 // with no slack left for either heuristic to reserve less
 {
  Net g = base_path();
  g.deadline = 2.5;
  const double opt = 2 * 1 + 2 * 3;
  for( char h : heurs ) {
   const auto st = check_dcr_spt( g , h , opt , "DCR_SPT deadline met "
                                  "with equality" );
   check( st == DCR::OK , std::string( "DCR_SPT deadline met with "
          "equality: ERA-" ) + ( h == '1' ? "I" : "H" ) + " says " +
          std::to_string( st ) + " instead of OK" );
   }
  }

 // a minimum rate above every capacity
 {
  Net g = base_path();
  g.rate = 3;
  for( char h : heurs )
   check( check_dcr_spt( g , h , INF , "DCR_SPT rate above capacity" ) ==
          DCR::Infeasible , "DCR_SPT rate above capacity: not Infeasible" );
  }

 // parallel arcs, a tie in the delay, and zero-cost arcs: 0 -> 1 twice with
 // the same capacity, one of cost 0 and one of cost 4, then 1 -> 2 of cost
 // 0; ERA-I prices both at full rate, and the cheaper one costs nothing
 {
  Net g;
  g.n = 3;
  g.arc( 0 , 1 , 2 , 4 ); g.arc( 0 , 1 , 2 , 0 ); g.arc( 1 , 2 , 4 , 0 );
  g.ndel.assign( 3 , 0 );
  g.s = 0; g.t = 2;
  g.burst = 1;
  g.deadline = 10;
  const double opt = dcr_optimum( g );
  check( opt == 0 , "DCR_SPT parallel arcs: optimum not 0" );
  for( char h : heurs )
   check_dcr_spt( g , h , opt , "DCR_SPT parallel arcs" );
  }

 // arcs closed with DCRcloseArcs(): the shortest-delay path 0 -> 2, of
 // capacity 8, is closed, and whatever is returned has to avoid it
 {
  Net g;
  g.n = 3;
  g.arc( 0 , 2 , 8 , 1 ); g.arc( 0 , 1 , 4 , 1 ); g.arc( 1 , 2 , 4 , 1 );
  g.ndel.assign( 3 , 0 );
  g.s = 0; g.t = 2;
  g.burst = 1;
  g.deadline = 10;
  for( char h : heurs ) {
   DCRData d;
   fill( d , g );
   DCR_SPT slv;
   slv.DCRloadProblem( 3 , 3 , 1 , &d.flow , d.links.data() ,
                       d.nodes.data() , g.MTU , DCR::SRP );
   slv.DCRsetHeur( h );
   int closed = 0;
   slv.DCRcloseArcs( &closed , 1 );
   const auto st = slv.DCRsolve();
   const std::string who = std::string( "DCR_SPT closed arc, ERA-" ) +
                           ( h == '1' ? "I" : "H" );
   check( st == DCR::OK , who + ": not OK" );
   if( st == DCR::OK ) {
    int X[ 2 ];
    double R[ 2 ];
    const int nh = slv.DCRgetUBSolPath( 0 , 0 , X , R );
    check( ( nh == 2 ) && ( X[ 0 ] == 2 ) && ( X[ 1 ] == 1 ) ,
           who + ": the path uses the closed arc" );
    check( slv.DCRgetObj() < INF , who + ": objective " +
           std::to_string( slv.DCRgetObj() ) );
    }
   }
  }

 // reloading: a feasible instance after an infeasible one on the same object
 {
  Net bad = base_path();
  bad.deadline = 1;
  Net good = base_path();
  DCRData db , dg;
  fill( db , bad );
  fill( dg , good );
  DCR_SPT slv;
  slv.DCRloadProblem( 3 , 2 , 1 , &db.flow , db.links.data() ,
                      db.nodes.data() , 1 , DCR::SRP );
  check( slv.DCRsolve() == DCR::Infeasible , "DCR_SPT reload: first" );
  slv.DCRloadProblem( 3 , 2 , 1 , &dg.flow , dg.links.data() ,
                      dg.nodes.data() , 1 , DCR::SRP );
  check( ( slv.DCRsolve() == DCR::OK ) && near( slv.DCRgetObj() , 8 ) ,
         "DCR_SPT reload: second" );
  }

 // what DCR_SPT does not support throws
 {
  Net g = base_path();
  DCRData d;
  fill( d , g );
  DCR_SPT slv;
  bool threw = false;
  try { slv.DCRloadProblem( 3 , 2 , 1 , &d.flow , d.links.data() ,
                            d.nodes.data() , 1 , DCR::WRT ); }
  catch( DCR::DCRException & ) { threw = true; }
  check( threw , "DCR_SPT WRT does not throw" );
  threw = false;
  try { slv.DCRsetHeur( '3' ); }
  catch( DCR::DCRException & ) { threw = true; }
  check( threw , "DCR_SPT heuristic 3 does not throw" );
  }
 }

/*--------------------------------------------------------------------------*/
/*-------------------------- DCRLagrangianSolver ---------------------------*/
/*--------------------------------------------------------------------------*/
/* Solves g with DCRLagrangianSolver for the given r_min and checks the two
 * sides: the Lagrangian value is not above the optimum of the problem it
 * relaxes (rates in [ r_min , cap ], burst / r_min), and the heuristic value
 * is the cost of the rates it returns, which meet the deadline and are not
 * below the optimum of the DCR problem on the same reduced graph (burst /
 * min r). If withupd, r_min is set by updrmin() before Solve(). Returns the
 * status. */

static DCRLagrangianSolver::LAGStatus check_lagr( const Net & g , double rmin ,
                                                  const std::string & what ,
                                                  bool withupd = false )
{
 DCRData d;
 fill( d , g );
 DCRLagrangianSolver slv;
 slv.LoadProblem( g.n , g.m() , d.flow , d.links.data() , d.nodes.data() ,
                  g.MTU , withupd ? rmin * 2 : rmin );
 if( withupd )
  slv.updrmin( rmin );
 slv.Solve();
 const auto st = slv.getStatus();

 double opt_lb = INF , opt_true = INF;
 for( const auto & P : all_paths( g , g.s , g.t , [ & ]( int a ) {
       return( g.cap[ a ] >= rmin ); } ) ) {
  opt_lb = std::min( opt_lb , path_opt_lb( g , P , rmin ) );
  opt_true = std::min( opt_true , path_opt_true( g , P , rmin ) );
  }

 if( opt_lb == INF ) {
  check( st == DCRLagrangianSolver::Infeasible , what + ": status " +
         std::to_string( st ) + " on an infeasible instance" );
  return( st );
  }

 check( st == DCRLagrangianSolver::OK , what + ": status " +
        std::to_string( st ) + " on a feasible instance" );
 if( st != DCRLagrangianSolver::OK )
  return( st );

 check( slv.getObjVal() <= opt_lb * ( 1 + 1e-7 ) + 1e-9 , what +
        ": Lagrangian value " + std::to_string( slv.getObjVal() ) +
        " above the optimum " + std::to_string( opt_lb ) );

 const double hv = slv.getHeurVal();
 check( hv < INF , what + ": no heuristic solution" );
 if( hv < INF ) {
  check( hv >= opt_true * ( 1 - 1e-7 ) - 1e-9 , what + ": heuristic value "
         + std::to_string( hv ) + " below the optimum " +
         std::to_string( opt_true ) );
  // the heuristic takes the burst term as burst / min r, which is not above
  // the burst / r_min of the relaxed problem: the two values are ordered
  // only if there is no burst
  if( g.burst == 0 )
   check( hv >= slv.getObjVal() * ( 1 - 1e-7 ) - 1e-9 , what +
          ": heuristic value below the Lagrangian one" );
  // the rates of the heuristic solution: a path meeting the deadline
  const auto HR = slv.getHeurRSOLS();
  Path P;
  std::vector< double > RP;
  int u = g.s;
  for( int k = 0 ; ( k < g.n ) && ( u != g.t ) ; ++k )
   for( int a = 0 ; a < g.m() ; ++a )
    if( ( g.sn[ a ] == u ) && ( HR[ a ] > 0 ) ) {
     P.push_back( a );
     RP.push_back( HR[ a ] );
     u = g.en[ a ];
     break;
     }
  check( is_path( g , P , g.s , g.t ) , what + ": heuristic rates are not "
         "an s-t path" );
  double c = 0;
  for( size_t k = 0 ; k < P.size() ; ++k )
   c += g.cost[ P[ k ] ] * RP[ k ];
  check( near( c , hv , 1e-7 ) , what + ": heuristic value " +
         std::to_string( hv ) + " is not the cost " + std::to_string( c ) +
         " of its rates" );
  check( path_delay( g , P , RP ) <= g.deadline * ( 1 + 1e-6 ) , what +
         ": heuristic rates violate the deadline" );
  }

 return( st );
 }

/*--------------------------------------------------------------------------*/
/// random instances, r_min at every capacity of the instance

static void test_lagrangian_random( void )
{
 std::mt19937 gen( 7 );
 int nok = 0 , ntot = 0;
 for( int k = 0 ; k < 150 ; ++k ) {
  Net g = random_net( gen , 5 , 9 );
  for( int a = 0 ; a < g.m() ; ++a )
   g.cost[ a ] += 1;  // positive costs: see test_lagrangian_edges()
  for( double rmin : { 1.0 , 2.0 , 4.0 } ) {
   const auto st = check_lagr( g , rmin , "Lagrangian random #" +
                               std::to_string( k ) + ", r_min " +
                               std::to_string( rmin ) , k % 2 );
   nok += ( st == DCRLagrangianSolver::OK );
   ++ntot;
   }
  }
 check( ( nok > ntot / 10 ) && ( nok < ntot ) , "Lagrangian random: " +
        std::to_string( nok ) + " feasible out of " + std::to_string( ntot ) );
 }

/*--------------------------------------------------------------------------*/
/// the edge cases of DCRLagrangianSolver

static void test_lagrangian_edges( void )
{
 // source equal to sink: the empty path, if burst / r_min meets the deadline
 {
  Net g = base_path();
  g.s = g.t = 1;
  check( check_lagr( g , 1 , "Lagrangian source == sink" ) ==
         DCRLagrangianSolver::OK , "Lagrangian source == sink: not OK" );
  }

 // no path at all
 {
  Net g = base_path();
  g.sn[ 1 ] = 2; g.en[ 1 ] = 1;
  check_lagr( g , 1 , "Lagrangian no path" );
  }

 // no path in the reduced graph: every capacity below r_min
 check_lagr( base_path() , 3 , "Lagrangian r_min above capacity" );

 // the only path violates the deadline
 {
  Net g = base_path();
  g.deadline = 2.4;
  check_lagr( g , 2 , "Lagrangian deadline violated" );
  }

 // the deadline met with equality at full rate, r_min = capacity
 {
  Net g = base_path();
  g.deadline = 2.5;
  check_lagr( g , 2 , "Lagrangian deadline met with equality" );
  }

 // a loose deadline: lambda = 0 is optimal, every rate at r_min
 {
  Net g = base_path();
  g.deadline = 100;
  DCRData d;
  fill( d , g );
  DCRLagrangianSolver slv;
  slv.LoadProblem( 3 , 2 , d.flow , d.links.data() , d.nodes.data() , 1 ,
                   0.5 );
  slv.Solve();
  check( ( slv.getStatus() == DCRLagrangianSolver::OK ) &&
         near( slv.getObjVal() , 2 ) && near( slv.getHeurVal() , 2 ) &&
         ( slv.getLambda() == 0 ) , "Lagrangian loose deadline" );
  }

 // a deadline that binds, so that lambda > 0 at the optimum
 {
  Net g = base_path();
  g.deadline = 4;
  check_lagr( g , 1 , "Lagrangian binding deadline" );
  check_lagr( g , 1 , "Lagrangian binding deadline, updrmin()" , true );
  }

 // zero-cost arcs: at lambda = 0 any rate is as cheap on them, and the
 // path through them has to be found all the same
 {
  Net g = base_path();
  g.cost[ 0 ] = 0;
  check_lagr( g , 1 , "Lagrangian zero-cost arc, loose deadline" );
  g.deadline = 4;
  check_lagr( g , 1 , "Lagrangian zero-cost arc, binding deadline" );
  g.cost[ 1 ] = 0;
  check_lagr( g , 1 , "Lagrangian all arcs of zero cost" );
  }

 // parallel arcs and a tie: two arcs 0 -> 1 of the same data
 {
  Net g = base_path();
  g.arc( 0 , 1 , 2 , 1 );
  g.deadline = 4;
  check_lagr( g , 1 , "Lagrangian parallel arcs" );
  }

 // closing and reopening arcs, and changing r_min, then solving again
 {
  Net g;
  g.n = 3;
  g.arc( 0 , 2 , 8 , 3 ); g.arc( 0 , 1 , 4 , 1 ); g.arc( 1 , 2 , 4 , 1 );
  g.ndel.assign( 3 , 0 );
  g.s = 0; g.t = 2;
  g.burst = 1;
  g.deadline = 10;
  DCRData d;
  fill( d , g );
  DCRLagrangianSolver slv;
  slv.LoadProblem( 3 , 3 , d.flow , d.links.data() , d.nodes.data() , 1 ,
                   1 );
  slv.Solve();
  check( near( slv.getHeurVal() , 2 ) , "Lagrangian before closing: " +
         std::to_string( slv.getHeurVal() ) );
  slv.closeArcs( { 1 } , 1 );             // 0 -> 2 alone, at r_min 2
  slv.updrmin( 2 );
  slv.Solve();
  check( ( slv.getStatus() == DCRLagrangianSolver::OK ) &&
         near( slv.getHeurVal() , 6 ) , "Lagrangian after closing: " +
         std::to_string( slv.getHeurVal() ) );
  int arc = 1;
  slv.openArcs( &arc , 1 );
  slv.updrmin( 1 );
  slv.Solve();
  check( near( slv.getHeurVal() , 2 ) , "Lagrangian after reopening: " +
         std::to_string( slv.getHeurVal() ) );
  slv.updrmin( 5 );                       // 0 -> 2 alone again
  slv.Solve();
  check( ( slv.getStatus() == DCRLagrangianSolver::OK ) &&
         near( slv.getHeurVal() , 15 ) , "Lagrangian r_min 5: " +
         std::to_string( slv.getHeurVal() ) );
  }
 }

/*--------------------------------------------------------------------------*/
/*--------------------------- SingleFlowDCRBlock ---------------------------*/
/*--------------------------------------------------------------------------*/
/// loads g into blk, with its source and sink

static void load( SingleFlowDCRBlock & blk , const Net & g )
{
 Block::Subset sn( g.m() ) , en( g.m() );
 for( int a = 0 ; a < g.m() ; ++a ) {
  sn[ a ] = g.sn[ a ] + 1;
  en[ a ] = g.en[ a ] + 1;
  }
 blk.load( g.n , g.m() , sn , en , g.cap , g.cost , g.ndel , g.ldel ,
           g.burst , g.deadline , g.MTU , g.rate );
 blk.chg_st( g.s , g.t );
 }

/*--------------------------------------------------------------------------*/
/// the feasibility of given solutions, the deadline met with equality

static void test_block_feasibility( void )
{
 // 0 -> 1 -> 2 and 0 -> 2
 Net g = base_path();
 g.arc( 0 , 2 , 1 , 1 );
 g.rate = 1;
 g.deadline = 2.5;  // met with equality by 0 -> 1 -> 2 at full rate
 SingleFlowDCRBlock blk;
 load( blk , g );

 const std::vector< double > X = { 1 , 1 , 0 } , R = { 2 , 2 , 0 };
 check( blk.flow_feasible( 0 , X ) , "Block flow_feasible() of a path" );
 check( blk.bound_feasible( 0 , X , R ) , "Block bound_feasible()" );
 check( blk.link_feasible( 0 , X , R ) , "Block link_feasible()" );
 check( blk.delay_feasible( 0 , X , R ) , "Block delay_feasible() with "
        "equality" );

 blk.set_x( X.begin() );
 blk.set_r( R.begin() );
 check( blk.is_feasible() , "Block is_feasible() of the path" );
 auto sol = blk.get_Solution( nullptr , false );
 check( blk.is_sol_feasible( sol ) , "Block is_sol_feasible() of the path" );
 delete sol;

 // a little less than the capacity on one arc: the deadline is violated
 const std::vector< double > R2 = { 2 , 1.99 , 0 };
 check( ! blk.delay_feasible( 0 , X , R2 ) , "Block delay_feasible() over "
        "the deadline" );
 check( blk.delay_feasible( 0.01 , X , R2 ) , "Block delay_feasible() "
        "within the tolerance" );

 // over the capacity, under rho, not a path, a rate on an unused arc
 check( ! blk.bound_feasible( 0 , X , { 2 , 3 , 0 } ) ,
        "Block bound_feasible() over the capacity" );
 check( ! blk.link_feasible( 0 , X , { 2 , 0.5 , 0 } ) ,
        "Block link_feasible() under rho" );
 check( ! blk.flow_feasible( 0 , { 1 , 0 , 0 } ) ,
        "Block flow_feasible() of half a path" );
 check( ! blk.link_feasible( 0 , X , { 2 , 2 , 1 } ) ,
        "Block link_feasible() with a rate on an unused arc" );

 // the other path, of one hop of capacity 1: 1 + 1 + burst 1 = 3 > 2.5
 check( ! blk.delay_feasible( 0 , { 0 , 0 , 1 } , { 0 , 0 , 1 } ) ,
        "Block delay_feasible() of the slow path" );
 }

/*--------------------------------------------------------------------------*/
/// load() with the optional data left empty, as its documentation allows

static void test_block_defaults( void )
{
 SingleFlowDCRBlock blk;
 const Block::Subset sn = { 1 , 2 } , en = { 2 , 3 };
 blk.load( 3 , 2 , sn , en );
 check( ( blk.get_NNodes() == 3 ) && ( blk.get_NArcs() == 2 ) ,
        "Block defaults: sizes" );
 check( blk.get_U().empty() && ( blk.get_U( 1 ) == Inf< double >() ) ,
        "Block defaults: capacities not infinite" );
 check( ( blk.get_C( 0 ) == 0 ) && ( blk.get_C( 1 ) == 0 ) ,
        "Block defaults: costs not 0" );
 check( blk.get_NodeDelays().empty() && blk.get_LinkDelays().empty() ,
        "Block defaults: delays not empty" );

 // a capacity given to an arc of a Block whose capacities are all infinite
 blk.chg_ucap( 3 , 0 );
 check( ( blk.get_U( 0 ) == 3 ) && ( blk.get_U( 1 ) == Inf< double >() ) ,
        "Block chg_ucap() on infinite capacities: U[ 0 ] = " +
        std::to_string( blk.get_U( 0 ) ) );
 }

/*--------------------------------------------------------------------------*/
/// entry i of a vector of delays, empty meaning all 0

static double delay_of( const std::vector< double > & d , Block::Index i )
{
 return( d.empty() ? 0 : d[ i ] );
 }

/*--------------------------------------------------------------------------*/
/// the data of two Blocks are the same

static bool same_data( const SingleFlowDCRBlock & a ,
                       const SingleFlowDCRBlock & b )
{
 if( ( a.get_NNodes() != b.get_NNodes() ) ||
     ( a.get_NArcs() != b.get_NArcs() ) )
  return( false );
 for( Block::Index i = 0 ; i < a.get_NArcs() ; ++i )
  if( ( a.get_SN( i ) != b.get_SN( i ) ) || ( a.get_EN( i ) != b.get_EN( i ) )
      || ( a.get_U( i ) != b.get_U( i ) ) || ( a.get_C( i ) != b.get_C( i ) )
      || ( delay_of( a.get_LinkDelays() , i ) !=
           delay_of( b.get_LinkDelays() , i ) ) )
   return( false );
 for( Block::Index i = 0 ; i < a.get_NNodes() ; ++i )
  if( ( a.get_B( i ) != b.get_B( i ) ) ||
      ( delay_of( a.get_NodeDelays() , i ) !=
        delay_of( b.get_NodeDelays() , i ) ) )
   return( false );
 return( ( a.get_FlowBurst() == b.get_FlowBurst() ) &&
         ( a.get_FlowDeadline() == b.get_FlowDeadline() ) &&
         ( a.get_MTU() == b.get_MTU() ) && ( a.get_rho() == b.get_rho() ) );
 }

/*--------------------------------------------------------------------------*/
/// the netCDF round trip, and a deserialize() over a loaded Block

static void test_block_serialize( void )
{
 const std::string file = "SingleFlowDCRBlock_unit_test.nc4";

 std::mt19937 gen( 11 );
 for( int k = 0 ; k < 5 ; ++k ) {
  Net g = random_net( gen , 6 , 10 );
  if( g.m() == 0 )
   g.arc( 0 , g.n - 1 , 1 , 1 );
  if( k == 1 )   // all capacities infinite: U is not written
   g.cap.assign( g.m() , Inf< double >() );
  if( k == 2 ) { // no delays: they are not written
   g.ldel.assign( g.m() , 0 );
   g.ndel.assign( g.n , 0 );
   }
  SingleFlowDCRBlock blk;
  load( blk , g );
  blk.Block::serialize( file , eBlockFile );
  auto copy = dynamic_cast< SingleFlowDCRBlock * >(
   Block::deserialize( file ) );
  check( copy != nullptr , "Block round trip #" + std::to_string( k ) +
         ": not a SingleFlowDCRBlock" );
  if( copy ) {
   check( same_data( blk , *copy ) , "Block round trip #" +
          std::to_string( k ) + ": data differ" );
   check( blk.is_feasible_instance() == copy->is_feasible_instance() ,
          "Block round trip #" + std::to_string( k ) + ": feasibility" );
   delete copy;
   }
  }

 // a file with no delays, read into a Block that had them
 {
  Net g = base_path();
  g.ldel.clear();
  g.ndel.clear();
  SingleFlowDCRBlock blk;
  load( blk , g );
  blk.Block::serialize( file , eBlockFile );

  g.ldel = { 0.5 , 0.25 };
  g.ndel = { 0.125 , 0 , 0 };
  SingleFlowDCRBlock old;
  load( old , g );
  {
   netCDF::NcFile f( file , netCDF::NcFile::read );
   old.deserialize( f.getGroup( "Block_0" ) );
   }
  check( same_data( blk , old ) , "Block deserialize() over a loaded Block "
         "keeps its delays" );
  }

 std::remove( file.c_str() );
 }

/*--------------------------------------------------------------------------*/
/// is_feasible_instance() against the enumeration, on random instances

static void test_block_feasible_instance( void )
{
 std::mt19937 gen( 5 );
 for( int k = 0 ; k < 200 ; ++k ) {
  const Net g = random_net( gen , 5 , 9 );
  SingleFlowDCRBlock blk;
  load( blk , g );
  bool feasible = false;
  for( const auto & P : all_paths( g , g.s , g.t ) )
   feasible = feasible || full_feasible( g , P );
  check( blk.is_feasible_instance() == feasible , "Block "
         "is_feasible_instance() #" + std::to_string( k ) + " says " +
         ( feasible ? "infeasible" : "feasible" ) );
  }

 // the only path meeting the deadline goes through a closed arc
 Net g;
 g.n = 3;
 g.arc( 0 , 2 , 8 , 1 ); g.arc( 0 , 1 , 1 , 1 ); g.arc( 1 , 2 , 1 , 1 );
 g.ndel.assign( 3 , 0 );
 g.s = 0; g.t = 2;
 g.deadline = 1;   // 0 -> 2 takes 0.25, 0 -> 1 -> 2 takes 4
 SingleFlowDCRBlock blk;
 load( blk , g );
 check( blk.is_feasible_instance() , "Block is_feasible_instance() before "
        "closing" );
 blk.close_arc( 0 );
 check( ! blk.is_feasible_instance() , "Block is_feasible_instance() with "
        "the only feasible path closed" );
 blk.open_arc( 0 );
 check( blk.is_feasible_instance() , "Block is_feasible_instance() after "
        "reopening" );
 }

/*--------------------------------------------------------------------------*/
/* A Solver that counts the physical Modification it receives, by type, and
 * solves the SingleFlowDCRBlock again with ERA-I on its current data, the
 * closed arcs left out. */

class ERAITestSolver : public Solver
{
 public:

 int count[ 5 ] = { 0 , 0 , 0 , 0 , 0 };  // per SingleFlowDCRBlockMod type
 double value = INF;                      // the cost ERA-I finds

 int compute( bool changedvars = true ) override {
  for( ; ! v_mod.empty() ; v_mod.pop_front() )
   if( auto dm = std::dynamic_pointer_cast< const SingleFlowDCRBlockMod >(
        v_mod.front() ) )
    ++count[ dm->type() ];

  auto blk = static_cast< SingleFlowDCRBlock * >( f_Block );
  Net g;
  g.n = blk->get_NNodes();
  for( Block::Index a = 0 ; a < blk->get_NArcs() ; ++a )
   if( ! blk->is_closed( a ) )
    g.arc( blk->get_SN( a ) - 1 , blk->get_EN( a ) - 1 , blk->get_U( a ) ,
           blk->get_C( a ) );
  g.ndel.assign( g.n , 0 );
  for( int i = 0 ; i < g.n ; ++i ) {
   if( blk->get_B( i ) < 0 ) g.s = i;
   if( blk->get_B( i ) > 0 ) g.t = i;
   }
  g.MTU = blk->get_MTU();
  g.burst = blk->get_FlowBurst();
  g.deadline = blk->get_FlowDeadline();
  g.rate = blk->get_rho();
  DCRData d;
  fill( d , g );
  DCR_SPT slv;
  slv.DCRloadProblem( g.n , g.m() , 1 , &d.flow , d.links.data() ,
                      d.nodes.data() , g.MTU , DCR::SRP );
  slv.DCRsetHeur( '1' );
  const bool ok = ( slv.DCRsolve() == DCR::OK );
  value = ok ? slv.DCRgetObj() : INF;
  return( ok ? kOK : kInfeasible );
  }

 void get_var_solution( Configuration * solc = nullptr ) override {}

 private:

 const std::string & private_name( void ) const override {
  static const std::string name = "ERAITestSolver";
  return( name );
  }
 };

/*--------------------------------------------------------------------------*/
/* The changes of the SingleFlowDCRBlock, each issuing its physical
 * Modification, seen by the Solver attached to it when it solves again.
 * The instance: 0 -> 3 (capacity 8, cost 1), the fastest, and 0 -> 1 -> 3
 * and 0 -> 2 -> 3 (capacity 4, cost 1 each), which ERA-I prices at full
 * rate: 8 the first, 8 each of the other two. ERA-I takes, among the paths
 * of least delay for each candidate r_min, the cheapest one. */

static void test_block_modifications( void )
{
 Net g;
 g.n = 4;
 g.arc( 0 , 3 , 8 , 1 );
 g.arc( 0 , 1 , 4 , 1 ); g.arc( 1 , 3 , 4 , 1 );
 g.arc( 0 , 2 , 4 , 1 ); g.arc( 2 , 3 , 4 , 1 );
 g.ndel.assign( 4 , 0 );
 g.s = 0; g.t = 3;
 g.burst = 1;
 g.deadline = 10;
 SingleFlowDCRBlock blk;
 load( blk , g );

 ERAITestSolver slv;
 blk.register_Solver( &slv );
 slv.compute();
 check( near( slv.value , 8 ) , "Block Mod: initial value " +
        std::to_string( slv.value ) );

 blk.chg_cost( 3 , 0 );                   // 0 -> 3 costs 24
 slv.compute();
 check( slv.count[ SingleFlowDCRBlockMod::eChgCost ] == 1 ,
        "Block Mod: chg_cost() issues no eChgCost" );
 check( near( slv.value , 24 ) , "Block Mod: value after chg_cost() " +
        std::to_string( slv.value ) );

 blk.close_arc( 0 );                      // 0 -> 1 -> 3 or 0 -> 2 -> 3
 slv.compute();
 check( slv.count[ SingleFlowDCRBlockMod::eCloseArc ] == 1 ,
        "Block Mod: close_arc() issues no eCloseArc" );
 check( near( slv.value , 8 ) , "Block Mod: value after close_arc() " +
        std::to_string( slv.value ) );

 blk.open_arc( 0 );
 slv.compute();
 check( slv.count[ SingleFlowDCRBlockMod::eOpenArc ] == 1 ,
        "Block Mod: open_arc() issues no eOpenArc" );
 check( near( slv.value , 24 ) , "Block Mod: value after open_arc() " +
        std::to_string( slv.value ) );

 // 0 -> 1 -> 3 at capacity 32 becomes the fastest, at 64 in all
 const std::vector< double > caps = { 32 , 32 };
 blk.chg_ucaps( caps.begin() , Block::Range( 1 , 3 ) );
 slv.compute();
 check( slv.count[ SingleFlowDCRBlockMod::eChgCaps ] == 1 ,
        "Block Mod: chg_ucaps() issues no eChgCaps" );
 check( near( slv.value , 64 ) , "Block Mod: value after chg_ucaps() " +
        std::to_string( slv.value ) );

 blk.chg_st( 1 , 3 );                     // from node 1: 1 -> 3 alone
 slv.compute();
 check( slv.count[ SingleFlowDCRBlockMod::eChgDfct ] > 0 ,
        "Block Mod: chg_st() issues no eChgDfct" );
 check( near( slv.value , 32 ) , "Block Mod: value after chg_st() " +
        std::to_string( slv.value ) );

 blk.chg_st( 3 , 0 );                     // no path back
 check( slv.compute() == Solver::kInfeasible , "Block Mod: no path back" );

 blk.unregister_Solver( &slv );
 }

/*--------------------------------------------------------------------------*/
/// fixing a flow Variable is closing the arc, and only that

static void test_block_fixing( void )
{
 Net g;
 g.n = 2;
 g.arc( 0 , 1 , 8 , 1 );
 g.arc( 0 , 1 , 4 , 1 );
 g.ndel.assign( 2 , 0 );
 g.s = 0; g.t = 1;
 g.burst = 1;
 g.deadline = 10;
 SingleFlowDCRBlock blk;
 load( blk , g );
 blk.generate_abstract_variables();
 auto x0 = blk.i2p_x( 0 );

 x0->set_value( 0 );
 x0->is_fixed( true );
 check( blk.is_closed( 0 ) , "Block fixing: fixed at 0, closed" );
 x0->is_fixed( false );
 check( ! blk.is_closed( 0 ) , "Block fixing: unfixed, open" );

 bool refused = false;
 x0->set_value( 1 );
 try { x0->is_fixed( true ); }
 catch( std::invalid_argument & ) { refused = true; }
 check( refused , "Block fixing: fixed at 1, refused" );
 x0->is_fixed( false );
 }

/*--------------------------------------------------------------------------*/
/*--------------------------------- MAIN -----------------------------------*/
/*--------------------------------------------------------------------------*/

int main( void )
{
 test_spt_random();
 test_spt_edges();
 test_dcr_spt_random();
 test_dcr_spt_edges();
 test_lagrangian_random();
 test_lagrangian_edges();
 test_block_feasibility();
 test_block_defaults();
 test_block_serialize();
 test_block_feasible_instance();
 test_block_modifications();
 test_block_fixing();

 if( failures )
  std::cout << failures << " checks failed" << std::endl;
 else
  std::cout << "All tests passed!!" << std::endl;

 return( failures );
 }

/*--------------------------------------------------------------------------*/
/*------------------------- End File test.cpp ------------------------------*/
/*--------------------------------------------------------------------------*/
