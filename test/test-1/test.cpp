/*--------------------------------------------------------------------------*/
/*-------------------------- File test.cpp ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Main for testing SingleFlowDCRBlock and SingleFlowDCRBendersSolver.
 *
 * Reads an instance of a MCF from a file (in either DIMACS or netCDF format)
 * in an SingleFlowDCRBlock, and from there in an object of a class MCFC derived from
 * MCFClass, as decided by the macro WHICH_MCF. Then, a SingleFlowDCRBendersSolver< MCFC > is
 * attached to the SingleFlowDCRBlock. The MCF problem is then repeatedly solved with
 * several changes in costs / capacities / deficits, arcs openings / closures
 * and arcs additions / deletions. The same operations are performed on the
 * two solvers, and the results are compared. This mostly tests SingleFlowDCRBlock and
 * SingleFlowDCRBendersSolver, since the actual MCFClass solved is the same, and so it can
 * easily be wrong in the same way for both the objects.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni
 */
/*--------------------------------------------------------------------------*/
/*------------------------------ DEFINES -----------------------------------*/
/*--------------------------------------------------------------------------*/

#define LOG_LEVEL 1
// 0 = only pass/fail
// 1 = result of each test

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/* Defines which :MCFClass solver is included and the corresponding version
 * of SingleFlowDCRBendersSolver< :MCFClass > is tested:
 *
 * - 0      for the CS2 class
 *
 * - 1      for the MCFCplex class
 *
 * - 2      for the MCFSimplex class
 *
 * - 3      for the MCFZIB class
 *
 * - 4      for the RelaxIV class
 *
 * - 5      for the SPTree class; note that SPTree cannot solve
 *          most MCF instances, except those with SPT structure */

#define WHICH_MCF 1

#if( LOG_LEVEL >= 1 )
#define LOG1( x ) cout << x
#define CLOG1( y , x ) if( y ) cout << x
#else
#define LOG1( x )
#define CLOG1( y , x )
#endif

#define USECOLORS 1
#if( USECOLORS )
#define RED( x ) "\x1B[31m" #x "\033[0m"
#define GREEN( x ) "\x1B[32m" #x "\033[0m"
#else
#define RED( x ) #x
#define GREEN( x ) #x
#endif

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <fstream>
#include <iomanip>

#include <random>

#if WHICH_MCF == 0

 #include "CS2.h"
 #define MCFC CS2

#elif WHICH_MCF == 1

 #include "MCFCplex.h"
 #define MCFC MCFCplex

#elif WHICH_MCF == 2

 #include "MCFSimplex.h"
 #define MCFC MCFSimplex

#elif WHICH_MCF == 3

 #include "MCFZIB.h"
 #define MCFC MCFZIB

#elif WHICH_MCF == 4

 #include "RelaxIV.h"
 #define MCFC RelaxIV

#elif WHICH_MCF == 5

 #include "SPTree.h"
 #define MCFC SPTree

#endif

#include "SingleFlowDCRBendersSolver.h"
#include "UpdateSolver.h"

/*--------------------------------------------------------------------------*/
/*-------------------------------- MACROS ----------------------------------*/
/*--------------------------------------------------------------------------*/

#define arg_2_str( s ) #s
#define solver_name( snm ) "SingleFlowDCRBendersSolver<" arg_2_str( snm ) ">"

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace MCFClass_di_unipi_it;
using namespace SMSpp_di_unipi_it;

// FIXME: Avoid these declarations
template<> const std::vector< int > SingleFlowDCRBendersSolver< MCFC >::Solver_2_MCFClass_int;
template<> const std::vector< int > SingleFlowDCRBendersSolver< MCFC >::Solver_2_MCFClass_dbl;

/*--------------------------------------------------------------------------*/
/*------------------------- TYPES & CONSTEXPRS -----------------------------*/
/*--------------------------------------------------------------------------*/
// check that MCFClass types and SingleFlowDCRBlock types agree

using Index = Block::Index;
using c_Index = Block::c_Index;
static constexpr Index IInf = SMSpp_di_unipi_it::Inf< Index >();

static_assert( std::is_same< Index , MCFClass::Index >::value );

using FNumber = SingleFlowDCRBlock::FNumber;
static constexpr FNumber FInf = SMSpp_di_unipi_it::Inf< FNumber >();

static_assert( std::is_same< FNumber , MCFClass::FNumber >::value );

using CNumber = SingleFlowDCRBlock::CNumber;
static constexpr CNumber CInf = SMSpp_di_unipi_it::Inf< CNumber >();

static_assert( std::is_same< CNumber , MCFClass::CNumber >::value );

using Range = Block::Range;
using c_Range = Block::c_Range;

using Subset = Block::Subset;
using c_Subset = Block::c_Subset;

static constexpr double BA = 1e-12;  // base accuracy of the SingleFlowDCRBendersSolver

/*--------------------------------------------------------------------------*/
/*------------------------------- GLOBALS ----------------------------------*/
/*--------------------------------------------------------------------------*/

unsigned int mode = 0;         // what is modified, what is solved

SingleFlowDCRBlock * oMCFB = nullptr;    // original SingleFlowDCRBlock
SingleFlowDCRBlock * dMCFB = nullptr;    // "derived" (i.e., R3B) SingleFlowDCRBlock
SingleFlowDCRBlock * sMCFB = nullptr;    // SingleFlowDCRBlock that is solved
SingleFlowDCRBlock * mMCFB = nullptr;    // SingleFlowDCRBlock that is modified

MCFClass * mcf;                // the MCFClass object

bool isnc4 = false;            // true if the file is a ntCDF one

std::mt19937 rg;               // base random generator
std::uniform_real_distribution<> dis( 0.0 , 1.0 );

FNumber MaxC = 0;              // max absolute value of costs
FNumber MaxU = 0;              // max absolute value of capacities / deficits

/*--------------------------------------------------------------------------*/
/*------------------------------ FUNCTIONS ---------------------------------*/
/*--------------------------------------------------------------------------*/

template< class T >
static void Str2Sthg( const char* const str , T &sthg )
{
 istringstream( str ) >> sthg;
 }

/*--------------------------------------------------------------------------*/
// returns a random number between 0.5 and 2, with 50% probability of
// being < 1

static double rndfctr( void )
{
 double fctr = dis( rg ) - 0.5;
 return( fctr < 0 ? - fctr : fctr * 4 );
 }

/*--------------------------------------------------------------------------*/
// generates a random k-vector of unique integers in 0 ... m - 1

static Subset GenerateRand( Index m , Index k )
{
 Subset rnd( m );
 std::iota( rnd.begin() , rnd.end() , 0 );
 std::shuffle( rnd.begin() , rnd.end() , rg );
 rnd.resize( k );
 sort( rnd.begin() , rnd.end() );

 return( rnd );
 }

/*--------------------------------------------------------------------------*/

static void CreateProb( unsigned int Optns )
{
 bool reoptmz = Optns & 1u;
 Optns /= 2;

 mcf = nullptr;  // unknown solver, or the required solver is not
                 // available due to the macros settings

 #if WHICH_MCF == 0  //- - - - - - - - - - - - - - - - - - - - - - - - - - -

  mcf = new CS2();
  LOG1( "CS2" );
  assert( false );  // CS2 not fully supported yet

 #elif WHICH_MCF == 1  //- - - - - - - - - - - - - - - - - - - - - - - - - -

  MCFCplex *cpx = new MCFCplex();

  if( Optns >= 0 )
   cpx->SetPar( CPX_PARAM_NETPPRIIND , int( Optns ) );
  mcf = cpx;
  LOG1( "MCFCplex" );

 #elif WHICH_MCF == 2  //- - - - - - - - - - - - - - - - - - - - - - - - - -

  auto *mcfs = new MCFSimplex();
  bool PrmlSmplx = Optns & 1u;
  char Prcng = 0;
  switch( Optns / 2 ) {
   case( 0 ): Prcng = char( MCFSimplex::kDantzig ); break;
   case( 1 ): Prcng = char( MCFSimplex::kFirstEligibleArc ); break;
   default:   Prcng = char( MCFSimplex::kCandidateListPivot );
   }
  if( ( ! PrmlSmplx ) && ( Prcng == MCFSimplex::kDantzig ) )
   Prcng = char( MCFSimplex::kCandidateListPivot );
  mcf = mcfs;
  LOG1( "MCFSimplex" );

 #elif WHICH_MCF == 3  //- - - - - - - - - - - - - - - - - - - - - - - - - -

  MCFZIB *zib = new MCFZIB();
  bool PrmlSmplx = Optns & 1;
  char Prcng;
  switch( Optns / 2 ) {
   case( 0 ): Prcng = char( MCFZIB::kDantzig ); break;
   case( 1 ): Prcng = char( MCFZIB::kFrstElA ); break;
   default:   Prcng = char( MCFZIB::kMltPrPr );
   }
  if( ( ! PrmlSmplx ) && ( Prcng == MCFZIB::kDantzig ) )
   Prcng = char( MCFZIB::kMltPrPr );
  zib->SetAlg( PrmlSmplx , Prcng );
  mcf = zib;
  LOG1( "MCFZIB" );
  assert( false );  // MCFZIB not fully supported yet

 #elif WHICH_MCF == 4  //- - - - - - - - - - - - - - - - - - - - - - - - - -

  RelaxIV *rlx = new RelaxIV();
  #if( AUCTION )
   if( Optns )
    rlx->SetPar( RelaxIV::kAuction , MCFClass::kYes );
  #endif
  mcf = rlx;
  LOG1( "RelaxIV" );

 #elif WHICH_MCF == 5  //- - - - - - - - - - - - - - - - - - - - - - - - - -

  mcf = new SPTree();
  LOG1( "SPTree" );
  assert( false );  // SPTree not fully supported yet

 #endif

 if( ! reoptmz )
  mcf->SetPar( MCFClass::kReopt , MCFClass::kNo );

 if( ! mcf )
  throw( std::logic_error( "No MCFClass defined" ) );
 }

/*--------------------------------------------------------------------------*/

static void load( char * fn )
{
 try {
  // note that usually the "original" SingleFlowDCRBlock is loaded, unless mode == 2
  // (==> original solved, R3 modified) *and* the R3 as been constructed
  // already, in which case the R3 is loaded; this is perhaps stretching
  // the concept of "R3Block" close to the breaking point, but in this case
  // it works because the R3B is a copy *and* always the same instance is
  // loaded

  if( isnc4 ) {
   netCDF::NcFile f( fn , netCDF::NcFile::read );
   if( f.isNull() ) {
    std::cerr << "cannot open nc4 file " << fn << std::endl;
    exit( 1 );
    }

   netCDF::NcGroupAtt gtype = f.getAtt( "SMS++_file_type" );
   if( gtype.isNull() ) {
    std::cerr << fn << " is not an SNS++ nc4 file" << std::endl;
    exit( 1 );
    }

   int type = 0;
   gtype.getValues( & type );

   if( type != eBlockFile ) {
    std::cerr << fn << " is not an SNS++ nc4 Block file" << std::endl;
    exit( 1 );
    }

   netCDF::NcGroup bg = f.getGroup( "Block_0" );
   if( bg.isNull() ) {
    std::cerr << "Block_0 empty or undefined in " << fn << std::endl;
    exit( 1 );
    }
   
   if( ( ( mode & 3u ) == 2 ) && dMCFB ) {
    dMCFB->deserialize( bg );  // load the (derived) SingleFlowDCRBlock

    // load the MCFClass out of the SingleFlowDCRBlock using the in-memory interface
    mcf->LoadNet( dMCFB->get_MaxNNodes() , dMCFB->get_MaxNArcs() ,
		  dMCFB->get_NNodes() , dMCFB->get_NArcs() ,
		  dMCFB->get_U().empty() ? nullptr : dMCFB->get_U().data() ,
		  dMCFB->get_C().empty() ? nullptr : dMCFB->get_C().data() ,
		  dMCFB->get_B().empty() ? nullptr : dMCFB->get_B().data() ,
		  dMCFB->get_SN().data() , dMCFB->get_EN().data() );
    }
   else {
    oMCFB->deserialize( bg );  // load the (original) SingleFlowDCRBlock

    // load the MCFClass out of the SingleFlowDCRBlock using the in-memory interface
    mcf->LoadNet( oMCFB->get_MaxNNodes() , oMCFB->get_MaxNArcs() ,
		  oMCFB->get_NNodes() , oMCFB->get_NArcs() ,
		  oMCFB->get_U().empty() ? nullptr : oMCFB->get_U().data() ,
		  oMCFB->get_C().empty() ? nullptr : oMCFB->get_C().data() ,
		  oMCFB->get_B().empty() ? nullptr : oMCFB->get_B().data() ,
		  oMCFB->get_SN().data() , oMCFB->get_EN().data() );
    }
   }
  else {
   ifstream iFile( fn );
   if( ! iFile ) {
    cerr << "Can't open dmx file " << fn << endl;
    exit( 1 );
    }

   mcf->LoadDMX( iFile );  // load the MCFClass

   iFile.clear();
   iFile.seekg( 0 );       // rewind the file
/*
   if( ( ( mode & 3u ) == 2 ) && dMCFB )
    iFile >> *dMCFB;       // load the (derived) SingleFlowDCRBlock
   else
    iFile >> *oMCFB;       // load the (original) SingleFlowDCRBlock
*/
  string fn1 = fn;
  ifstream iFile1( fn1.substr(0,fn1.find_last_of('.'))+".dcr" );
  if( ! iFile1 ) {
    cerr << "Can't open dcr file " << fn << endl;
    exit( 1 );
    }
  if( ( ( mode & 3u ) == 2 ) && dMCFB ) {
   dMCFB->load( iFile );
   Index NNodes = dMCFB->get_NNodes();
   Index NArcs = dMCFB->get_NArcs();
   try {
    dMCFB->load_dcr( iFile1 , NNodes , NArcs );
    }
   catch(...) {
    cerr << "Error: dcr file error!" << endl;
    exit( 1 );
   }
  }
  else {
   oMCFB->load( iFile );
   Index NNodes = oMCFB->get_NNodes();
   Index NArcs = oMCFB->get_NArcs();
   try {  
    oMCFB->load_dcr( iFile1 , NNodes , NArcs );
    }
   catch(...) {
    cerr << "Error: dcr file error!" << endl;
    exit( 1 );
   }
   }
  }
  //if( mode & 4u ) {
   // if so instructed, generate abstract representation for oMCFB
   oMCFB->generate_abstract_constraints();
   oMCFB->generate_dynamic_constraints();
   oMCFB->generate_objective();
   //}

  //if( ( mode & 8u ) && dMCFB ) {
  if( dMCFB ) {
   // if so instructed, generate abstract representation for dMCFB (if any)
   dMCFB->generate_abstract_constraints();
   dMCFB->generate_dynamic_constraints();
   dMCFB->generate_objective();
   }
  }
 catch( exception &e ) {
  cerr << "MCFClass: " << e.what() << endl;
  exit( 1 );
  }
 catch(...) {
  cerr << "Error: unknown exception thrown" << endl;
  exit( 1 );
  }

 }  // end( load )

/*--------------------------------------------------------------------------*/

static bool SolveMCF( void ) 
{
 try {
  // solve the MCFClass- - - - - - - - - - - - - - - - - - - - - - - - - - - -
  mcf->SolveMCF();
  auto stat = mcf->MCFGetStatus();

  // solve the SingleFlowDCRBlock- - - - - - - - - - - - - - - - - - - - - - - - - - - -
  Solver * slvr = (sMCFB->get_registered_solvers()).front();
  int rtrn = slvr->compute( false );

  if( ( stat == MCFClass::kOK ) &&
      ( rtrn >= Solver::kOK ) && ( rtrn < Solver::kError ) ) {
   auto fo1 = mcf->MCFGetFO();
   auto fo2 = slvr->get_ub();
   if( abs( fo1 - fo2 )
       <= 1e-9 *  max( double( 1 ) , abs( max( fo1 , fo2 ) ) ) ) {
    LOG1( "OK(f)" << endl );
    return( true );
    }
   }

  if( ( stat == MCFClass::kUnfeasible ) &&
      ( rtrn == Solver::kInfeasible ) ) {
   LOG1( "OK(e)" << endl );
   return( false );
   }

  if( ( stat == MCFClass::kUnbounded ) &&
      ( rtrn == Solver::kUnbounded ) ) {
   LOG1( "OK(u)" << endl );
   return( false );
   }

  #if( LOG_LEVEL >= 1 )
   cout << "MCFClass = ";
   switch( mcf->MCFGetStatus() ) {
    case( MCFClass::kOK ):         cout << mcf->MCFGetFO();
                                   break;
    case( MCFClass::kUnfeasible ): cout << "        +INF";
                                   break;
    case( MCFClass::kUnbounded ):  cout << "        -INF";
                                   break;
    default:                       cout << "      Error!";
    }

   cout << " ~ SingleFlowDCRBlock = ";
  #endif

  if( ( rtrn >= Solver::kOK ) && ( rtrn < Solver::kError ) ) {
   LOG1( slvr->get_ub() << endl );
   return( false );
   }

  #if( LOG_LEVEL >= 1 )
   switch( rtrn ) {
    case( Solver::kInfeasible ):   cout << "        +INF";
                                   break;
    case( Solver::kUnbounded ):    cout << "        -INF";
                                   break;
    default:                       cout << "      Error!";
    }

   cout << endl;
  #endif

  return( false );
  }
 catch( exception &e ) {
  cerr << e.what() << endl;
  exit( 1 );
  }
 catch(...) {
  cerr << "Error: unknown exception thrown" << endl;
  exit( 1 );
  }
 }

/*--------------------------------------------------------------------------*/

int main( int argc , char **argv )
{
 // reading command line parameters - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 long int seed = 1;
 unsigned int wchg = 127;
 double p_change = 0.5;
 Index n_change = 10;
 Index n_repeat = 40;
 int optns = 1;

 switch( argc ) {
  case( 9 ): Str2Sthg( argv[ 8 ] , optns );
  case( 8 ): Str2Sthg( argv[ 7 ] , p_change );
  case( 7 ): Str2Sthg( argv[ 6 ] , n_change );
  case( 6 ): Str2Sthg( argv[ 5 ] , n_repeat );
  case( 5 ): Str2Sthg( argv[ 4 ] , wchg );
  case( 4 ): Str2Sthg( argv[ 3 ] , mode );
  case( 3 ): Str2Sthg( argv[ 2 ] , seed );
  case( 2 ): break;
  default: cerr << "Usage: " << argv[ 0 ] <<
	   " <dmx file> [seed mode wchg #rounds #chng %chng optns]"
		<< endl <<
	   "       mode: 0 = only one, 1 = orig -> solve, 2 = solve -> orig"
		<< endl <<
	   "             +4 = abstract orig, +8 = abstract solve,"
		<< endl <<
 	   "             +16 = also change abstract representation"
		<< endl <<
           "       wchg: what to change, coded bit-wise "
		<< endl <<
           "             0 = cost, 1 = cap, 2 = dfct, 3 = o.arc, 4 = c.arc"
		<< endl <<
           "             5 = delete arc, 6 = add arc"
		<< endl <<
	   "       optns: bit 0 = re-optimize, other bits MCF-specific" 
		<< endl <<
	   "              Relax   : > 0 uses Auction"
		<< endl <<
	   "              Cplex   : network pricing parameter"
		<< endl <<
	   "              ZIB     : 1st bit == 1 ==> primal +"
		<< endl <<
	   "                        0 = Dantzig, 2 = First Eligible, 4 = MPP"
	        << endl <<
	   "              Simplex : 1st bit == 1 ==> primal +"
		<< endl <<
	   "                        0 = Dantzig, 2 = First Eligible, 4 = MPP"
	        << endl;
	   return( 1 );
  }

 if( ( mode & 3u ) == 3 ) {
  std::cerr << "Wrong mode (must be 0, 1 or 2)" << std::endl;
  exit( 1 );
  }
  
 if( mode & 16u )  // if the abstract representations are changed
  mode |= 12u;     // ensure they exist in the first place

 // construction and loading of the objects - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // construct MCFClass- - - - - - - - - - - - - - - - - - - - - - - - - - - -

 CreateProb( optns );

 // construct the "original" SingleFlowDCRBlock - - - - - - - - - - - - - - - - - - - -
 // ... in a bit of a roundabout way by un-necessarily using the factory

 oMCFB = dynamic_cast< SingleFlowDCRBlock * >( Block::new_Block( "SingleFlowDCRBlock" ) );
 assert( oMCFB );

 // load the instance - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // first of check if the file is a .dmx (default) or a .nc4 one

 std::string name( argv[ 1 ] );
 if( name.size() > 4 ) {
  std::string sffx = name.substr( name.size() - 4 , 4 );
  std::string nc4( ".nc4" );

  isnc4 = std::equal( sffx.begin() , sffx.end() , nc4.begin() ,
		      []( auto a , auto b ) {
		       return( std::tolower( a ) == std::tolower( b ) ); } );
  }

 load( argv[ 1 ] );

 // if so instructed, construct the R3 SingleFlowDCRBlock = copy- - - - - - - - - - - -

 if( mode & 3u ) {
  dMCFB = dynamic_cast< SingleFlowDCRBlock * >( oMCFB->get_R3_Block() );
  assert( dMCFB );           // excess of caution (we know it is)

  if( ( mode & 8u ) && dMCFB ) {
   // if so instructed, also generate abstract representation for dMCFB
   dMCFB->generate_abstract_constraints();
   dMCFB->generate_objective();
   }

  if( ( mode & 3u ) == 1 ) {  // modify the original, solve the R3
   mMCFB = oMCFB;
   sMCFB = dMCFB;
   }
  else {                     // modify the R3, solve the original
   mMCFB = dMCFB;
   sMCFB = oMCFB;
   }

  // attach an UpdateSolver to the modified: if the modified is the original
  // then map_forward, otherwise (i.e., it is the copy) map_back
  mMCFB->register_Solver(
	      new UpdateSolver( sMCFB , nullptr , mMCFB != oMCFB ? 1 : 0 ) );
  }
 else                        // just use one SingleFlowDCRBlock
  sMCFB = mMCFB = oMCFB;

 // compute min/max cost & max deficit- - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Index n = mcf->MCFn();
 Index m = mcf->MCFm();

 LOG1( ", n = " << n << ", m = " << m << endl );
 if( n_change > m )
  n_change = m;

 CNumber c_max = - CInf;    // max cost
 CNumber c_abs = 0;         // max absolute value of cost
 CNumber c_min = - c_max;   // min cost
 FNumber u_max = 0;         // max (finite) capacity
 FNumber u_avg = 0;         // average capacity
 FNumber u_min = FInf;      // min capacity
 FNumber b_abs = 0;         // max absolute value of deficit

 for( Index i = 0 ; i < m ; ++i ) {
  auto ci = mcf->MCFCost( i );
  if( std::abs( ci ) > c_abs )
   c_abs = std::abs( ci );
  if( ci < CInf ) {
   if( ci < c_min )
    c_min = ci;
   if( ci > c_max )
    c_max = ci;
   }

  auto ui = mcf->MCFUCap( i );
  if( ui < FInf ) {
   if( ui > u_max )
    u_max = ui;
   u_avg += ui;
   if( ui < u_min )
    u_min = ui;
   }
  }

 u_avg /= m;

 for( Index i = 0 ; i < n ; )
  if( auto bi = mcf->MCFDfct( i++ ) )
   if( std::abs( bi ) > b_abs )
    b_abs = std::abs( bi );

 bool nzdfct = b_abs > 0;

 // set epsilons in MCFClass and SingleFlowDCRBendersSolver- - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // set epsilons in MCFClass
 mcf->SetPar( MCFClass::kEpsFlw ,
	      BA * std::max( u_max , std::max( b_abs , FNumber( 1 ) ) ) );
 mcf->SetPar( MCFClass::kEpsCst , BA * std::max( c_abs , CNumber( 1 ) ) );

 // attach a "true" SingleFlowDCRBendersSolver to the one that is actually solved
 // auto MCFS = Solver::new_Solver( solver_name( MCFC ) );
 auto MCFS = new SingleFlowDCRBendersSolver< MCFC >();

 sMCFB->register_Solver( MCFS );

 // set epsilons in SingleFlowDCRBendersSolver
 MCFS->set_par( Solver::dblAbsAcc ,
		BA * std::max( u_max , std::max( b_abs , FNumber( 1 ) ) ) );
 MCFS->set_par( CDASolver::dblAAccDSol ,
		BA * std::max( c_abs , CNumber( 1 ) ) ); 

 // first solver call - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 #if( LOG_LEVEL >= 1 )
  cout << "First call: ";
  cout.setf( ios::scientific, ios::floatfield );
  cout << setprecision( 6 );
 #endif

 auto OK = SolveMCF();
 
 // main loop - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // now, for n_repeat times:
 // - up tp n_change costs are changed, then the two problems are re-solved;
 // - up to n_change capacities are changed, then the two problems are
 //   re-solved, then the original capacities are restored;
 // - if the problem is not a circulation problem, 2 deficits are modified
 //   (adding and subtracting the same number), then the two problems are
 //   re-solved, then the original deficits are restored;
 // - up to n_change arcs are closed, then the two problems are re-solved;
 //   the same arcs are re-opened, then the two problems are re-solved

 rg.seed( seed );  // seed the pseudo-random number generator

 bool diffarcs = false;  // whether added arcs ended up with different names

 for( Index iter = 0 ; iter < n_repeat ; ++iter ) {

  LOG1( iter << ": " );

  // before making changes, lock() mMCFB: this is of course useless since
  // nothing else has it, but there you go. Use "mcf" as the "owner", since
  // it clearly cannot be a reserved address

  bool owned = mMCFB->is_owned_by( mcf );
  if( ( ! owned ) && ( ! mMCFB->lock( mcf ) ) )
   throw( std::logic_error( "can't lock mMCFB" ) );

  // change costs - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( ( wchg & 1u ) && ( dis( rg ) <= p_change ) ) {
   Index tochange = max( double( 1 ) , dis( rg ) * n_change );

   LOG1( tochange << " cost" );

   LinearFunction * lf = nullptr;
   if( ( mode & 16u ) && ( dis( rg ) < 0.5 ) )
    // change via abstract representation
    lf = static_cast< LinearFunction * >( static_cast< FRealObjective * >(
				  mMCFB->get_objective() )->get_function() );

   if( tochange == 1 ) {
    CNumber newcst = c_min + CNumber( dis( rg ) * ( c_max - c_min ) );
    Index arc = Index( dis( rg ) * ( m - 1 ) );

    mcf->ChgCost( arc , newcst );

    if( lf ) {  // change via abstract representation
     LOG1( "(a)" );
     lf->modify_coefficient( arc , newcst );
     }
    else  // change via call to chg_* method
     mMCFB->chg_cost( newcst , arc );

    LOG1( " - " );
    }
   else {
    SingleFlowDCRBlock::Vec_CNumber newcsts( tochange );
    for( Index i = 0 ; i < tochange ; ++i )
     newcsts[ i ] = c_min + CNumber( dis( rg ) * ( c_max - c_min ) );

    // in 50% of the cases do a ranged change, in the others a sparse change
    if( dis( rg ) <= 0.5 ) {
     Index strt = dis( rg ) * ( m - tochange );
     Index stp = strt + tochange;
     mcf->ChgCosts( newcsts.data() , nullptr , strt , stp );

     if( lf ) {  // change via abstract representation
      LOG1( "s(r,a) - " );
      lf->modify_coefficients( std::move( newcsts ) , Range( strt , stp ) );
      }
     else {  // change via call to chg_* method
      // in 50% of the cases a direct call, otherwise use the methods factory
      //!! if( dis( rg ) <= 0.5 ) {
       mMCFB->chg_costs( newcsts.begin() , Range( strt , stp ) );
       LOG1( "s(r) - " );
      //!!  }
      //!! else {
      //!!  std::string mthd_name = "SingleFlowDCRBlock::chg_costs";
      //!!  const auto * mthd = Block::get_method_fs( mthd_name ,
      //!!                                      Block::MS_dbl_rngd::args() );
      //!!  assert( mthd_name == Block::get_method_name_fs
      //!!          ( mthd , Block::MS_dbl_rngd::args() ) );
      //!!  std::invoke( *mthd , mMCFB , newcsts.begin() ,
       //!!              Range( strt , stp ) , eNoBlck , eNoBlck );
      //!!  LOG1( "s(r-mf) - " );
      //!!  }
      }
     }
    else {
     Subset nms( GenerateRand( m , tochange ) );
     nms.push_back( OPTtypes_di_unipi_it::Inf< MCFClass::Index >() );

     mcf->ChgCosts( newcsts.data() , nms.data() );
     nms.resize( tochange );

     if( lf ) {  // change via abstract representation
      LOG1( "s(s,a) - " );
      lf->modify_coefficients( std::move( newcsts ) , std::move( nms ) ,
			       true );
      }
     else {  // change via call to chg_* method
      // in 50% of the cases a direct call, otherwise use the methods factory
      //!! if( dis( rg ) <= 0.5 ) {
       mMCFB->chg_costs( newcsts.begin() , std::move( nms ) , true );
       LOG1( "s(s) - " );
      //!! }
      //!! else {
      //!!  std::string mthd_name = "SingleFlowDCRBlock::chg_costs";
      //!!  const auto * mthd = Block::get_method_fs( mthd_name,
      //!!                                      Block::MS_dbl_sbst::args() );
      //!!  assert( mthd_name == Block::get_method_name_fs
      //!!          ( mthd , Block::MS_dbl_sbst::args() ) );
      //!!  std::invoke( *mthd , mMCFB , newcsts.begin() , std::move( nms ) ,
       //!!              true , eNoBlck , eNoBlck );
      //!!  LOG1( "s(s-mf) - " );
      //!!  }
      }
     }
    }
   }  // end( if( change costs ) )

  // change capacities- - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( ( wchg & 2u ) && ( dis( rg ) <= p_change ) ) {
   Index tochange = max( double( 1 ) , dis( rg ) * n_change );
   LOG1( tochange << " capacit" );

   if( tochange == 1 ) {
    auto arc = Index( dis( rg ) * ( m - 1 ) );
    CNumber newcap = mcf->MCFUCap( arc ) * rndfctr();
    mcf->ChgUCap( arc , newcap );

    if( ( mode & 16u ) && ( dis( rg ) < 0.5 ) ) {
     // change via abstract representation
     LOG1( "y(a) - " );
     mMCFB->i2p_ub( arc )->set_rhs( newcap );
     }
    else {  // change via call to chg_* method
     mMCFB->chg_ucap( newcap , arc );
     LOG1( "y - " );
     }
    }
   else {
    SingleFlowDCRBlock::Vec_FNumber newcaps( tochange );

    // in 50% of the cases do a ranged change, in the others a sparse change
    if( dis( rg ) <= 0.5 ) {
     Index strt = dis( rg ) * ( m - tochange );
     Index stp = strt + tochange;
     for( Index i = 0 ; i < tochange ; ++i )
      newcaps[ i ] = mcf->MCFUCap( i + strt ) * rndfctr();
     mcf->ChgUCaps( newcaps.data() , nullptr , strt , stp );

     if( ( mode & 16u ) && ( dis( rg ) < 0.5 ) ) {
      // change via abstract representation, sending to a new channel
      LOG1( "ies(a,r) - " );
      auto chnl = mMCFB->open_channel();
      auto modpar = Observer::make_par( eModBlck , chnl );
      for( Index i = 0 ; i < tochange ; ++i )
       mMCFB->i2p_ub( i + strt )->set_rhs( newcaps[ i ] , modpar );
      mMCFB->close_channel( chnl );
      }
     else {  // change via call to chg_* method
      // in 50% of the cases a direct call, otherwise use the methods factory
      //!! if( dis( rg ) <= 0.5 ) {
       mMCFB->chg_ucaps( newcaps.begin() , Range( strt , stp ) );
       LOG1( "ies(r) - " );
      //!!  }
      //!! else {
      //!!  std::string mthd_name = "SingleFlowDCRBlock::chg_ucaps";
      //!!  const auto * mthd = Block::get_method_fs( mthd_name,
      //!!                                     Block::MS_dbl_rngd::args() );
      //!!  assert( mthd_name == Block::get_method_name_fs
      //!!          ( mthd , Block::MS_dbl_rngd::args() ) );
      //!!  std::invoke( *mthd , mMCFB , newcaps.begin() ,
       //!!              Range( strt , stp ) , eNoBlck , eNoBlck );
      //!!  LOG1( "ies(r-mf) - " );
      //!!  }
      }
     }
    else {
     Subset nms( GenerateRand( m , tochange ) );
     auto ncit = newcaps.begin();
     for( auto i : nms )
      *(ncit++) = mcf->MCFUCap( i ) * rndfctr();

     nms.push_back( OPTtypes_di_unipi_it::Inf< MCFClass::Index >() );
     mcf->ChgUCaps( newcaps.data() , nms.data() );
     nms.resize( tochange );

     if( ( mode & 16u ) && ( dis( rg ) < 0.5 ) ) {
      // change via abstract representation, sending to a new channel
      LOG1( "ies(a,s) - " );
      auto chnl = mMCFB->open_channel();
      auto modpar = Observer::make_par( eModBlck , chnl );
      for( Index i = 0 ; i < tochange ; ++i )
       mMCFB->i2p_ub( nms[ i ] )->set_rhs( newcaps[ i ] , modpar );
      mMCFB->close_channel( chnl );
      }
     else {  // change via call to chg_* method
      // in 50% of the cases a direct call, otherwise use the methods factory
      //!! if( dis( rg ) <= 0.5 ) {
       mMCFB->chg_ucaps( newcaps.begin() , std::move( nms ) , true );
       LOG1( "ies(s) - " );
      //!!  }
      //!! else {
      //!!  std::string mthd_name = "SingleFlowDCRBlock::chg_ucaps";
      //!!  const auto * mthd = Block::get_method_fs( mthd_name,
      //!!                                      Block::MS_dbl_sbst::args() );
      //!!  assert( mthd_name == Block::get_method_name_fs
      //!!          ( mthd , Block::MS_dbl_sbst::args() ) );
      //!!  std::invoke( *mthd , mMCFB , newcaps.begin() , std::move( nms ) ,
       //!!              true , eNoBlck , eNoBlck );
      //!!  LOG1( "ies(s-mf) - " );
      //!!  }
      }
     }
    }
   }  // end( if( change capacities ) )

  // change deficits- - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( ( wchg & 4u ) && ( dis( rg ) <= p_change ) ) {
 /*  
   LOG1( "2 deficits" );

   Index posn = 0;
   Index negn = 0;
   FNumber posd = NAN;
   FNumber negd = NAN;

   if( nzdfct ) {  // if there are nonzero deficits
    SingleFlowDCRBlock::Vec_FNumber dfcts( n );
    mcf->MCFDfcts( dfcts.data() );

    do
     posn = Index( dis( rg ) * n );  // select node with positive
    while( dfcts[ posn ] <= 0 );     // deficit (one must exist)
    posd = dfcts[ posn ];

    do
     negn = Index( dis( rg ) * n );  // select node with negative
    while( dfcts[ negn ] >= 0 );     // deficit (one must exist)
    negd = dfcts[ negn ];
    }
   else {
    posn = Index( dis( rg ) * n );   // just select at random
    negn = Index( dis( rg ) * n );
    posd = negd = 0;
    }

   FNumber Dlt = u_avg * 2 * dis( rg );
   if( dis( rg ) <= 0.5 ) {  // in 50% of cases up, in 50% of cases down
    posd += Dlt;
    negd -= Dlt;
    }
   else {
    Dlt = min( Dlt , max( max( posd , - negd ) / 2 , double( 1 ) ) );
    posd -= Dlt;
    negd += Dlt;
    }

   mcf->ChgDfct( posn , posd );
   mcf->ChgDfct( negn , negd );

   // pack the two Modification into a new channel
   auto chnl = mMCFB->open_channel();
   auto modpar = Observer::make_par( eModBlck , chnl );

   if( ( mode & 16u ) && ( dis( rg ) < 0.5 ) ) {
    // change via abstract representation
    LOG1( "(a)" );
    mMCFB->i2p_e( posn )->set_both( posd , modpar );
    mMCFB->i2p_e( negn )->set_both( negd , modpar );
    }
   else {  // change via call to chg_* method
    // note that eModBlck makes no sense for a physical Modification,
    // but SingleFlowDCRBlock is supposed to take care of this
    mMCFB->chg_dfct( posd , posn , modpar , modpar );
    mMCFB->chg_dfct( negd , negn , modpar , modpar );
    }

   mMCFB->close_channel( chnl );

   LOG1( " - " );
*/
   }  // end( change deficits )

  // closing arcs- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( ( wchg & 8u ) && ( dis( rg ) <= p_change ) ) {
   Index changed = 0;

   Subset nms( n_change );
   for( auto i = mMCFB->get_NStaticArcs() ; i < mMCFB->get_NArcs() ; ++i ) {
    if( mcf->IsDeletedArc( i ) )
     continue;
    if( mcf->IsClosedArc( i ) )
     continue;
    if( dis( rg ) <= 0.5 )
     continue;
    
    nms[ changed++ ] = i;
    mcf->CloseArc( i );

    if( changed >= n_change )
     break;
    }

   if( changed ) {
    nms.resize( changed );
    LOG1( changed << " close" );

    if( ( mode & 16u ) && ( dis( rg ) < 0.5 ) ) {
     // change via abstract representation
     LOG1( "(a)" );
     auto modpar = mMCFB->open_if_needed( eModBlck , changed );
     for( auto i : nms ) {
      auto *x = mMCFB->i2p_x( i );
      x->set_value( 0 );
      x->is_fixed( true , modpar );
      }
     mMCFB->close_if_needed( modpar , changed );
     }
    else {  // change via call to chg_* method
     // in 50% of the cases a direct call, otherwise use the methods factory
     //!! if( dis( rg ) <= 0.5 )
      mMCFB->close_arcs( std::move( nms ) );
     //!! else {
     //!!  std::string mthd_name = "SingleFlowDCRBlock::close_arcs";
     //!!  const auto * mthd = Block::get_method_fs( mthd_name,
     //!!                                          Block::MS_sbst::args() );
     //!!  assert( mthd_name == Block::get_method_name_fs
     //!!          ( mthd , Block::MS_sbst::args() ) );
     //!!  std::invoke( *mthd , mMCFB , std::move( nms ) , false ,
      //!!              eNoBlck , eNoBlck );
     //!!  LOG1( "(mf)" );
     //!!  }
     }

    LOG1( " - " );
    }
   }

  // re-opening arcs - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( ( wchg & 16u ) && ( dis( rg ) <= p_change ) ) {
   Index changed = 0;

   Subset nms( n_change );
   for( auto i = mMCFB->get_NStaticArcs() ; i < mMCFB->get_NArcs() ; ++i ) {
    if( ( mcf->IsDeletedArc( i ) ) || ( ! mcf->IsClosedArc( i ) ) )
     continue;
    if( dis( rg ) <= 0.5 )
     continue;

    nms[ changed++ ] = i;
    mcf->OpenArc( i );

    if( changed >= n_change )
     break;
    }

   if( changed ) {
    nms.resize( changed );
    LOG1( changed << " open" );

    if( ( mode & 16u ) && ( dis( rg ) < 0.5 ) ) {
     // change via abstract representation
     LOG1( "(a)" );
     auto modpar = mMCFB->open_if_needed( eModBlck , changed );
     for( auto i : nms )
      mMCFB->i2p_x( i )->is_fixed( false , modpar );
     mMCFB->close_if_needed( modpar , changed );
     }
    else {  // change via call to chg_* method
     // in 50% of the cases a direct call, otherwise use the methods factory
     //!! if( dis( rg ) <= 0.5 )
      mMCFB->open_arcs( std::move( nms ) );
     //!! else {
     //!!  std::string mthd_name = "SingleFlowDCRBlock::open_arcs";
     //!!  const auto * mthd = Block::get_method_fs( mthd_name,
     //!!                                            Block::MS_sbst::args() );
     //!!  assert( mthd_name == Block::get_method_name_fs
     //!!          ( mthd , Block::MS_sbst::args() ) );
     //!!  std::invoke( *mthd , mMCFB , std::move( nms ) , false ,
      //!!              eNoBlck , eNoBlck );
     //!!  LOG1( "(mf)" );
     //!!  }
     }

    LOG1( " - " );
    }
   }

  // deleting arcs - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( ( wchg & 32u ) && ( dis( rg ) <= p_change ) ) {
   Index changed = 0;

   if( dis( rg ) < 0.5 ) {
    // delete somewhere in the middle

    for( auto i = mMCFB->get_NStaticArcs() ; i < mMCFB->get_NArcs() ; ++i ) {
     if( mcf->IsDeletedArc( i ) )
      continue;
     if( dis( rg ) <= 0.75 )
      continue;

     mcf->DelArc( i );
     mMCFB->remove_arc( i );
     if( ++changed >= n_change )
      break;
     }

    CLOG1( changed , changed << " delete(m) - " );
    }
   else {
    for( auto i =  mMCFB->get_NArcs() ; --i >= mMCFB->get_NStaticArcs() ; ) {
     if( mcf->IsDeletedArc( i ) )
      continue;
     if( dis( rg ) <= 0.13 )
      break;

     mcf->DelArc( i );
     mMCFB->remove_arc( i );
     if( ++changed >= n_change )
      break;
     }

    CLOG1( changed , changed << " delete(e) - " );
    }
   }

  // creating new arcs - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( ( wchg & 64u ) && ( dis( rg ) <= p_change ) ) {
   Index changed = 0;
   Index afterend = 0;
   while( changed < n_change ) {
    if( dis( rg ) <= 0.13 )
     break;

    ++changed;

    // random sn != en
    Index sn = 0;
    Index en = 0;
    do {
     sn = dis( rg ) * mMCFB->get_NNodes() + 1;
     en = dis( rg ) * mMCFB->get_NNodes() + 1;
     } while( sn == en );

    // random cost in [ - c_max , c_max ]
    auto cst = c_max * ( 1 - 2 * dis( rg ) );

    // random capacity <= 0.75 u_avg
    auto cap = 1.5 * ( u_avg - u_min ) * dis( rg ) + u_min;

    auto arc = mMCFB->add_arc( sn , en , cst , cap );
    if( arc != mcf->AddArc( sn , en , cap , cst ) )
     diffarcs = true;

    if( arc >= m )
     ++afterend;

    if( mMCFB->get_NArcs() >= mMCFB->get_MaxNArcs() )
     break;
    }

   #if( LOG_LEVEL >= 1 )
    if( changed ) {
     cout << "create " << changed << "(" << afterend << ")";
     if( diffarcs )
      cout << "[d]";
     cout << " - ";
     }
   #endif
   }

  // check that the status of the arcs is the same- - - - - - - - - - - - - -

  if( wchg & 120u ) {  // ... if it can ever change

   n = mcf->MCFn();
   if( n != mMCFB->get_NNodes() ) {
    cerr << endl << "error: different number of nodes" << endl;
    exit( 1 );
    }

   m = mcf->MCFm();
   if( m != mMCFB->get_NArcs() ) {
    cerr << endl << "error: different number of arcs" << endl;
    exit( 1 );
    }

   if( ! diffarcs )
    for( Index i = 0 ; i < m ; ++i ) {
     if( mcf->IsDeletedArc( i ) != mMCFB->is_deleted( i ) ) {
      std::cerr << "inconsistent del status for arc " << i << std::endl;
      exit( 1 );
      }

     if( mcf->IsClosedArc( i ) != mMCFB->is_closed( i ) ) {
      std::cerr << "inconsistent cls status for arc " << i << std::endl;
      exit( 1 );
      }
     }
   }

  // since all changes are done, unlock mMCFB
  if( ! owned )
   mMCFB->unlock( mcf );

  // finally, re-solve the problems- - - - - - - - - - - - - - - - - - - - -
  // yet, if the problem is either unfeasible or unbounded, or something has
  // gone awry with the arcs names, re-load it in both MCFClass and SingleFlowDCRBlock

  auto OKi = SolveMCF();
  OK |= OKi;
  if( ( ! OKi ) || diffarcs ) {
   load( argv[ 1 ] );
   n = mcf->MCFn();
   m = mcf->MCFm();
   diffarcs = false;
   }
  }  // end( main loop )- - - - - - - - - - - - - - - - - - - - - - - - - - -
     // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( OK )
  cout << GREEN( All test passed!! ) << endl;
 else
  cout << RED( Shit happened!! ) << endl;

 // destroy objects - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // unregister (and delete) all Solvers attached to the SingleFlowDCRBlocks
 mMCFB->unregister_Solvers();
 sMCFB->unregister_Solvers();

 // delete the SingleFlowDCRBlocks
 delete dMCFB;
 delete oMCFB;

 // delete the :MCFClass object
 delete mcf;

 // terminate - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 return( 0 );

 }  // end( main )

/*--------------------------------------------------------------------------*/
/*------------------------ End File test.cpp -------------------------------*/
/*--------------------------------------------------------------------------*/
