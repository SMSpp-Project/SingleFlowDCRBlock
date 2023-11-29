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
#include <filesystem>  

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
MultiFlowDCRBlock * oMCFB2 = nullptr;    // original MultiFlowDCRBlock

/*--------------------------------------------------------------------------*/

int main( int argc , char **argv )
{
  namespace stdfs = std::filesystem;
  std::string file = argv[ 1 ];
  std::ifstream iNode( file.substr(0,file.find_last_of('.'))+".nod" );
  int commodity, node, arc;
  
  if( ! iNode.is_open() )
   throw( std::invalid_argument( "can't open file .nod" ) );

  std::ifstream iNode1( file.substr(0,file.find_last_of('.'))+"-1.nod" );
  if( ! iNode1.good() ){
    stdfs::path p = stdfs::current_path();
    stdfs::rename(file.substr(0,file.find_last_of('.'))+".nod", file.substr(0,file.find_last_of('.'))+"-1.nod");
    std::ifstream iNode2( file.substr(0,file.find_last_of('.'))+"-1.nod" );
    iNode2 >> commodity;
    iNode2 >> node >> arc;
    iNode2.close();
  } else {
    iNode1 >> commodity;
    iNode1 >> node >> arc;
    iNode1.close();
  }

  int index = 0;
  int k = 0;
  int mtu;
  string buffer1[commodity];

  std::ifstream iDCR1( file.substr(0,file.find_last_of('.'))+"-1.dcr" );
  if( ! iDCR1.good() ){
    stdfs::path p = stdfs::current_path();
    stdfs::rename(file.substr(0,file.find_last_of('.'))+".dcr", file.substr(0,file.find_last_of('.'))+"-1.dcr");
    std::ofstream iDCR( file.substr(0,file.find_last_of('.'))+".dcr" );
    std::ifstream iDCR1( file.substr(0,file.find_last_of('.'))+"-1.dcr" );
    std::string buffer;
    getline(iDCR1, buffer);
    iDCR << buffer << '\n';
    while(index<arc+node) {  
      getline(iDCR1, buffer);
      iDCR << buffer << '\n';
      index +=1 ;
    }
    while(index<arc+node+commodity) {  
      getline(iDCR1, buffer1[k]);
      k += 1;
      index += 1;
    }
    iDCR1 >> mtu;
  } else {
    stdfs::path p = stdfs::current_path();
    std::ofstream iDCR( file.substr(0,file.find_last_of('.'))+".dcr" );
    std::ifstream iDCR1( file.substr(0,file.find_last_of('.'))+"-1.dcr" );
    std::string buffer;
    getline(iDCR1, buffer);
    iDCR << buffer << '\n';
    while(index<arc+node) {  
      getline(iDCR1, buffer);
      iDCR << buffer << '\n';
      index +=1 ;
    }
    while(index<arc+node+commodity) {  
      getline(iDCR1, buffer1[k]);
      k += 1;
      index += 1;
    }
    iDCR1 >> mtu;
  }

  std::ofstream iParam( file.substr(0,file.find_last_of('.'))+".param" );
  iParam << mtu << "\n";
  for(index=0; index<commodity; index++) {  
    iParam << buffer1[index] << '\n';
  }
  iParam.close();
  int j;

  for(j=1; j<commodity+1; ++j){
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

    auto bsc2 = dynamic_cast< BlockSolverConfig * >(Configuration::deserialize( "MILPPar2.txt" ) );
    if( ! bsc ) {
    cerr << "Error: configuration file not a BlockSolverConfig" << endl;
    exit( 1 );    
    }

    oMCFB = dynamic_cast< MultiFlowDCRBlock * >( Block::new_Block( "MultiFlowDCRBlock" ) );
    assert( oMCFB );

    oMCFB1 = dynamic_cast< MultiFlowDCRBlock * >( Block::new_Block( "MultiFlowDCRBlock" ) );
    assert( oMCFB1 );

    oMCFB2 = dynamic_cast< MultiFlowDCRBlock * >( Block::new_Block( "MultiFlowDCRBlock" ) );
    assert( oMCFB2 );

    oMCFB->load( argv[ 1 ] );
    bsc->apply( oMCFB );
    bsc->clear();  // keep the clear()-ed BlockSolverConfig for final cleanup

    // check Solvers - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    if( oMCFB->get_registered_solvers().empty() ) {
      cerr << "Error: BlockSolverConfig did not register any Solver" << endl;
      exit( 1 );    
    }
      
    Solver * slvrCONT = oMCFB->get_registered_solvers().front();
    int rtrn = slvrCONT->compute();
      
    if( slvrCONT->has_var_solution() ){

      oMCFB1->load( argv[ 1 ] );

      //oMCFB1->generate_abstract_constraints();
      //oMCFB1->generate_abstract_variables();
      bsc1->apply( oMCFB1 );
      bsc1->clear(); 

      Solver * slvrLD = oMCFB1->get_registered_solvers().front();
      //std::ofstream lagrangian_log("lagrangian_log.txt");
      //slvrLD->set_log( &lagrangian_log );

      int rtrn1 = slvrLD->compute();

      oMCFB2->load( argv[ 1 ] );
      bsc2->apply( oMCFB2 );
      bsc2->clear();
      Solver * slvrEXACT = oMCFB2->get_registered_solvers().front();
      int rtrn2 = slvrEXACT->compute();

      oMCFB->print(std::cout);
      std::cout << "CONTINUOUS RELAX: " << slvrCONT->get_lb() << "\n";
      std::cout << "LAGRANGIAN RELAX: " << slvrLD->get_lb() << "\n";
      std::cout << "OPTIMAL VALUE: " << slvrEXACT->get_ub() << "\n\n";
    } else {
      break;
    }
  } 
 }  // end( main )

/*--------------------------------------------------------------------------*/
/*------------------------ End File test.cpp -------------------------------*/
/*--------------------------------------------------------------------------*/
