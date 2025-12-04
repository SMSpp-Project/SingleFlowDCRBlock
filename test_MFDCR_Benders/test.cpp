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

#include <chrono>

#include "MCFCplex.h"
#include "BlockSolverConfig.h"
#include "SingleFlowDCRBlock.h"
#include "MultiFlowDCRBlock.h"
#include "CPXMILPSolver.h"
#include "GRBMILPSolver.h"
#include "MILPSolver.h"
#include "UpdateSolver.h"
#include "LagrangianDualSolver.h"
#include "BundleSolver.h"
#include "SingleFlowDCRBendersSolver.h"

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

  //std::ofstream output_cpu( file.substr(0,file.find_last_of('.'))+"_cpu.txt" );
  //std::ofstream output_gap( file.substr(0,file.find_last_of('.'))+"_gap.txt" );
  //std::ofstream output( file.substr(0,file.find_last_of('.'))+".txt" );
  
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

  int no_continuous = 0;

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

  //output_cpu << "NC" << "\t\t" << "exact" << "\t\t" << "lagrangian" << "\t\t" << "continuous" << "\n";
  //output_gap << "NC" << "\t\t" << "exact" << "\t\t" << "lagrangian" << "\t\t" << "continuous" << "\n";

  double max_GAP = 0.0;
  double max_ERROR = 0.0;

  //for(j=52; j<53; ++j){ //./MFDCRBlock_test ods24/Garr200109.nod
  //for(j=44; j<=50; ++j){ //./MFDCRBlock_test ods24/Abilene.nod
  //for(j=1; j<5; ++j){ //./MFDCRBlock_test ods24/Garr200109.nod
  for(j=1; j<=100; ++j){
  //for(j=1; j<std::min(500,commodity+1); j+=std::min(500,commodity+1)/10){
    //std::cout << "commodities: " << j << std::endl;
    std::cout << "commodities: " << j << std::endl;
    std::ofstream fout( file );
    fout << j << " " << node << " " << arc << " " << arc;  
    fout.close();

    auto bsc = dynamic_cast< BlockSolverConfig * >(Configuration::deserialize( "BSPar1.txt" ) );
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
      
    std::clock_t c_start;
    std::clock_t c_end; 

    oMCFB2->load( argv[ 1 ] );
    bsc2->apply( oMCFB2 );
    bsc2->clear();
    Solver * slvrEXACT = oMCFB2->get_registered_solvers().front();
    auto c_start_chrono_exact = std::chrono::high_resolution_clock::now();//std::clock();
    int rtrn2 = slvrEXACT->compute();
    auto c_end_chrono_exact = std::chrono::high_resolution_clock::now();//std::clock();
    //double time_elapsed_exact = 1000.0 * (c_end-c_start) / (double) CLOCKS_PER_SEC;
    double time_elapsed_exact = std::chrono::duration_cast<std::chrono::nanoseconds>(c_end_chrono_exact - c_start_chrono_exact).count()/1e+9;
    std::cout << "OPTIMAL VALUE SOLUTION: " << slvrEXACT->get_ub() << "\n";
    std::cout << "CPU EXACT: " << time_elapsed_exact << "\n";

    Solver * slvrLDMILP = oMCFB->get_registered_solvers().front();
    //slvrLDMILP->set_log(&cout);
    LagrangianDualSolver * LDS = dynamic_cast< LagrangianDualSolver * >( slvrLDMILP );
    LDS->set_par( LDS->str_par_str2idx( "strLogFileName" ) , "LogFileLDS" );

    int index = 0;
/*
    for( const auto & sbi : oMCFB->get_nested_Blocks() ){
      index++;
      MILPSolver* solver = dynamic_cast< MILPSolver * >( sbi->get_registered_solvers().front() );
      solver->set_par( solver->str_par_str2idx( "strLogFileName" ) , "LogFileMILP_" + std::to_string(index));
      //solver->set_par( solver->str_par_str2idx( "strOutputFile" ) , "MILP_" + std::to_string(index));
    }
*/
    auto c_start_chrono_lagrangianMILP = std::chrono::high_resolution_clock::now();//std::clock();
    /////if(slvrEXACT->has_var_solution())
      int rtrn = slvrLDMILP->compute();
    auto c_end_chrono_lagrangianMILP = std::chrono::high_resolution_clock::now();//std::clock();
    //double time_elapsed_continuous = 1000.0 * (c_end-c_start) / (double) CLOCKS_PER_SEC;
    double time_elapsed_lagrangianMILP = std::chrono::duration_cast<std::chrono::nanoseconds>(c_end_chrono_lagrangianMILP - c_start_chrono_lagrangianMILP).count()/1e+9;

    std::cout << "LAGRANGIAN LB MILP: " << slvrLDMILP->get_lb() << "\n";
    std::cout << "LAGRANGIAN UB MILP: " << slvrLDMILP->get_ub() << "\n";
    std::cout << "CPU MILP: " << time_elapsed_lagrangianMILP  << "\n";
    std::cout << (slvrLDMILP->get_lb()-slvrLDMILP->get_ub())/slvrLDMILP->get_lb() << std::endl;
    if((slvrLDMILP->get_lb()-slvrLDMILP->get_ub())/slvrLDMILP->get_lb() > 1e-5){
      std::cout << "slvrLDMILP->get_lb() > slvrLDMILP->get_ub()" << std::endl;
      ///break;
    }

    oMCFB1->load( argv[ 1 ] );
    bsc1->apply( oMCFB1 );
    bsc1->clear(); 
    Solver * slvrLD = oMCFB1->get_registered_solvers().front();
    slvrLD->set_log(&cout);
    LagrangianDualSolver * LDSB = dynamic_cast< LagrangianDualSolver * >( slvrLD );
    LDSB->set_par( LDSB->str_par_str2idx( "strLogFileName" ) , "LogFileLDSBender" );
    auto c_start_chrono_lagrangian = std::chrono::high_resolution_clock::now();//std::clock();
    /////if(slvrEXACT->has_var_solution())
      int rtrn1 = slvrLD->compute();
    auto c_end_chrono_lagrangian = std::chrono::high_resolution_clock::now();//std::clock();
    //double time_elapsed_lagrangian = 1000.0 * (c_end-c_start) / (double) CLOCKS_PER_SEC;
    double time_elapsed_lagrangian = std::chrono::duration_cast<std::chrono::nanoseconds>(c_end_chrono_lagrangian - c_start_chrono_lagrangian).count()/1e+9;
///*
    //if( slvrEXACT->has_var_solution() ){
      //double laggapMILP  = std::abs((slvrEXACT->get_ub()-std::min(slvrLDMILP->get_lb(),slvrLDMILP->get_ub()))/std::min(slvrLDMILP->get_lb(),slvrLDMILP->get_ub()));
      //double laggap  = std::abs((slvrEXACT->get_ub()-std::min(slvrLD->get_lb(),slvrLD->get_ub()))/std::min(slvrLD->get_lb(),slvrLD->get_ub()));
      std::cout << "LAGRANGIAN LB Benders: " << slvrLD->get_lb() << "\n";
      std::cout << "LAGRANGIAN UB Benders: " << slvrLD->get_ub() << "\n";
      std::cout << "CPU Benders: " << time_elapsed_lagrangian  << "\n";
      std::cout << (slvrLD->get_lb()-slvrLD->get_ub())/slvrLD->get_lb() << std::endl;
      if((slvrLD->get_lb()-slvrLD->get_ub())/slvrLD->get_lb() > 1e-2){
        std::cout << "slvrLD->get_lb() > slvrLD->get_ub()" << std::endl;
        ///break;
    }

      std::cout << "\n";
      
      //std::cout << "DIFFERENCE in LAGRANGIAN LBs (abs): " << (slvrLDMILP->get_lb() - slvrLD->get_lb())/slvrLDMILP->get_lb() << "\n";
      std::cout << "DIFFERENCE in LAGRANGIAN LBs: " << (slvrLD->get_lb() - slvrLDMILP->get_lb())/slvrLDMILP->get_lb() << "\n";
    
      if((slvrLD->get_lb() - slvrLDMILP->get_lb())/slvrLDMILP->get_lb() > 0){
        if((slvrLD->get_lb() - slvrLDMILP->get_lb())/slvrLDMILP->get_lb() > max_GAP)
          max_GAP = (slvrLD->get_lb() - slvrLDMILP->get_lb())/slvrLDMILP->get_lb();
      } else {
        if((slvrLD->get_lb() - slvrLDMILP->get_lb())/slvrLDMILP->get_lb() < max_ERROR)
          max_ERROR = (slvrLD->get_lb() - slvrLDMILP->get_lb())/slvrLDMILP->get_lb();
      }
    //} 

    std::cout << "max_GAP = " << max_GAP << std::endl;
    std::cout << "max_ERROR = " << max_ERROR << std::endl << std::endl;
//*/
/*
  if( slvrEXACT->has_var_solution() ) 
    std::cout << file.substr(0,file.find_last_of('.')) << "\t" << commodity << "\t" << node << "\t" << arc << "\t" << j << "\t" << 
          std::abs(slvrEXACT->get_ub()-slvrEXACT->get_lb())/(slvrEXACT->get_lb()) << "\t" << 
              std::abs(slvrEXACT->get_ub()-slvrLD->get_lb())/(slvrLD->get_lb()) << 
                "\t" << std::abs(slvrLD->get_lb()-slvrEXACT->get_lb())/(slvrEXACT->get_lb()) << "\t" <<
                  std::abs(slvrEXACT->get_ub()-slvrLDMILP->get_lb())/(slvrLDMILP->get_lb()) << 
                    "\t" << std::abs(slvrLDMILP->get_lb()-slvrEXACT->get_lb())/(slvrEXACT->get_lb()) << "\t" <<
                      "\t" << round(time_elapsed_exact*100.0)/100.0 << "\t" << round(time_elapsed_lagrangian*100.0)/100.0 << 
                        "\t" << round(time_elapsed_lagrangianMILP*100.0)/100.0 << std::endl;
  else
    std::cout << file.substr(0,file.find_last_of('.')) << "\t" << commodity << "\t" << node << 
          "\t" << arc << "\t" << j << "\t" << "problem is infeasible" << std::endl;
  */
  }
  return 0;
 }  // end( main )

/*--------------------------------------------------------------------------*/
/*------------------------ End File test.cpp -------------------------------*/
/*--------------------------------------------------------------------------*/
