/*--------------------------------------------------------------------------*/
/*-------------------------- File test.cpp ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Main for testing 
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni
 */
/*--------------------------------------------------------------------------*/
/*-------------------------------- MACROS ----------------------------------*/
/*--------------------------------------------------------------------------*/

#define LOG_LEVEL 0
// 0 = only pass/fail
// 1 = result of each test
// 2 = + solver log

#if( LOG_LEVEL >= 1 )
 #define LOG1( x ) cout << x
 #define CLOG1( y , x ) if( y ) cout << x

 #if( LOG_LEVEL >= 2 )
  #define LOG_ON_COUT 1
  // if nonzero, the NDO Solver log is sent on cout rather than on a file
 #endif
#else
 #define LOG1( x )
 #define CLOG1( y , x )
#endif

/*--------------------------------------------------------------------------*/
// if nonzero, the :MILPSolver attached to the NCoCubeBlock is detached and
// re-attached to it at all iterations

#define DETACH_LP 0

/*--------------------------------------------------------------------------*/

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
#include <sstream>
#include <iomanip>

#include <random>

#include "MCFCplex.h"
#include "BlockSolverConfig.h"
#include "SingleFlowDCRBlock.h"
#include "CPXMILPSolver.h"
#include "MILPSolver.h"
#include "UpdateSolver.h"

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace std;
using namespace MCFClass_di_unipi_it;
using namespace SMSpp_di_unipi_it;

SingleFlowDCRBlock * oMCFB = nullptr;    // original SingleFlowDCRBlock
MCFClass * mcf;                          // the MCFClass object

static void load( char * fn ){

 ifstream iFile( fn );
 if( ! iFile ) {
   cerr << "Can't open dmx file " << fn << endl;
   exit( 1 );
  }

 //mcf->LoadDMX( iFile );  // load the MCFClass

 iFile.clear();
 iFile.seekg( 0 );       // rewind the file

 string fn1 = fn;
 ifstream iFile1( fn1.substr(0,fn1.find_last_of('.'))+".dcr" );
 if( ! iFile1 ) {
   cerr << "Can't open dcr file " << fn << endl;
   exit( 1 );
  }
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

/*--------------------------------------------------------------------------*/

int main( int argc , char **argv )
{

// load the problem- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

 // attach the Solver to the Block- - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // do it by using a single a BlockSolverConfig, read from file
 
 auto bsc = dynamic_cast< BlockSolverConfig * >(
		               Configuration::deserialize( "MILPPar.txt" ) );
 if( ! bsc ) {
  cerr << "Error: configuration file not a BlockSolverConfig" << endl;
  exit( 1 );    
  }

 oMCFB = dynamic_cast< SingleFlowDCRBlock * >( Block::new_Block( "SingleFlowDCRBlock" ) );
 assert( oMCFB );

 load( argv[ 1 ] );

 oMCFB->generate_abstract_constraints();
 oMCFB->generate_dynamic_constraints();
 oMCFB->generate_objective();
 
 bsc->apply( oMCFB );
 bsc->clear();  // keep the clear()-ed BlockSolverConfig for final cleanup

 // check Solvers - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( oMCFB->get_registered_solvers().empty() ) {
  cerr << "Error: BlockSolverConfig did not register any Solver" << endl;
  exit( 1 );    
  }

 Solver * slvr = oMCFB->get_registered_solvers().front();
 int rtrn = slvr->compute( false );
 
 }  // end( main )

/*--------------------------------------------------------------------------*/
/*------------------------ End File test.cpp -------------------------------*/
/*--------------------------------------------------------------------------*/
