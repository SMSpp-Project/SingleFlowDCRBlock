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
 * \author Laura Galli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Luca Mencarelli \n
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
#include "MultiFlowDCRBlock.h"
#include "CPXMILPSolver.h"
#include "MILPSolver.h"
#include "UpdateSolver.h"
#include "LagrangianDualSolver.h"
#include "BundleSolver.h"

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace std;
//using namespace MCFClass_di_unipi_it;
using namespace SMSpp_di_unipi_it;

MultiFlowDCRBlock * oMCFB = nullptr;    // original MultiFlowDCRBlock
MultiFlowDCRBlock * oMCFB1 = nullptr;    // original MultiFlowDCRBlock

/*--------------------------------------------------------------------------*/

int main( int argc , char **argv )
{
  std::string file = argv[ 1 ];
  std::ifstream iNode( file.substr(0,file.find_last_of('.'))+"-1.nod" );
  if( ! iNode.is_open() )
   throw( std::invalid_argument( "can't open file .nod" ) );
  int commodity, node, arc;
  iNode >> commodity;
  iNode >> node;
  iNode >> arc;
  iNode >> commodity;

  int j;

  for(j=1; j<commodity; ++j){
    std::ofstream fout( file );
    fout << j << " " << node << " " << arc << " " << arc;  
    fout.close();

    auto bsc = dynamic_cast< BlockSolverConfig * >(Configuration::deserialize( "MILPPar.txt" ) );
    if( ! bsc ) {
    cerr << "Error: configuration file not a BlockSolverConfig" << endl;
    exit( 1 );    
    }

    auto bsc1 = dynamic_cast< BlockSolverConfig * >(Configuration::deserialize( "BSPar.txt" ) );
    if( ! bsc1 ) {
    cerr << "Error: configuration file not a BlockSolverConfig" << endl;
    exit( 1 );    
    }

    oMCFB = dynamic_cast< MultiFlowDCRBlock * >( Block::new_Block( "MultiFlowDCRBlock" ) );
    assert( oMCFB );

    oMCFB1 = dynamic_cast< MultiFlowDCRBlock * >( Block::new_Block( "MultiFlowDCRBlock" ) );
    assert( oMCFB1 );

    oMCFB->load( argv[ 1 ] );
    bsc->apply( oMCFB );
    bsc->clear();  // keep the clear()-ed BlockSolverConfig for final cleanup

    // check Solvers - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    if( oMCFB->get_registered_solvers().empty() ) {
      cerr << "Error: BlockSolverConfig did not register any Solver" << endl;
      exit( 1 );    
    }
      
    Solver * slvr = oMCFB->get_registered_solvers().front();
    int rtrn = slvr->compute();
      
    if( slvr->has_var_solution() ){

      oMCFB1->load( argv[ 1 ] );

      //oMCFB1->generate_abstract_constraints();
      //oMCFB1->generate_abstract_variables();
      bsc1->apply( oMCFB1 );
      bsc1->clear(); 

      Solver * slvrLD = oMCFB1->get_registered_solvers().front();
      //std::ofstream lagrangian_log("lagrangian_log.txt");
      //slvrLD->set_log( &lagrangian_log );

      int rtrn1 = slvrLD->compute();
      //std::cout << slvrLD->get_ub() << "\n";

      double CPX_LB = slvr->get_lb();
      double CPX_UB = slvr->get_lb();
      
      double LD_sol = slvrLD->get_lb();

      oMCFB->print(std::cout);
      std::cout << std::abs( CPX_LB-LD_sol ) << "\n";
      std::cout << std::abs( CPX_UB-LD_sol ) << "\n";
    } else {
      break;
    }
  } 
 }  // end( main )

/*--------------------------------------------------------------------------*/
/*------------------------ End File test.cpp -------------------------------*/
/*--------------------------------------------------------------------------*/
