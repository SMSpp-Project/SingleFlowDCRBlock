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

/*--------------------------------------------------------------------------*/

static void load( string fndmx, string fndcr ){

 //iFile.clear();
 //iFile.seekg( 0 );       // rewind the file

 ifstream fndmx1(fndmx);
 ifstream fndcr1(fndcr);
 oMCFB->load( fndmx1 );
 Index NNodes = oMCFB->get_NNodes();
 Index NArcs = oMCFB->get_NArcs();
 try {  
   oMCFB->load_dcr( fndcr1 , NNodes , NArcs );
  } 
  catch(...) {
    cerr << "Error: dcr file error!" << endl;
    exit( 1 );
   }
}

/*--------------------------------------------------------------------------*/

int main( int argc , char **argv )
{

 int i, j, source, destination, commodity;
 int SN, EN, numberArc, numberComm;
 int NComm, NNodes, NArcs;
 float rho, rho1, mtu, FlowBursts, FlowDeadline;
 float capacity, cost;
 string fn = argv[ 1 ];

 ifstream iNode( fn.substr(0,fn.find_last_of('.'))+".nod" );
 ifstream iFile2( fn.substr(0,fn.find_last_of('.'))+".sup" );
 ifstream iFile3( fn.substr(0,fn.find_last_of('.'))+".param" );

 string fn1 = fn;

 iNode >> NComm; 
 iNode >> NNodes;
 iNode >> NArcs;

 iFile3 >> mtu;

 std::ofstream output( fn.substr(0,fn.find_last_of('.'))+"_feasible_flows.txt" );
 std::ofstream output_sup( fn.substr(0,fn.find_last_of('.'))+"_sup.txt" );
 std::ofstream output_param( fn.substr(0,fn.find_last_of('.'))+"_param.txt" );
 
 int index_feas = 0;

 for(j=0; j<std::min(1000,NComm); j++){

  ifstream iFile1( fn1.substr(0,fn1.find_last_of('.'))+".dcr" );
  ifstream iArc( fn.substr(0,fn.find_last_of('.'))+".arc" );

  ofstream foutdmx("output.dmx");
  ofstream foutdcr("output.dcr");

  iFile2 >> source;
  iFile2 >> commodity;
  iFile2 >> rho;
  iFile2 >> destination;
  iFile2 >> commodity;
  iFile2 >> rho1;

  iFile3 >> FlowBursts >> FlowDeadline;

  while (!iFile1.eof()) {
   string buffer;
   getline(iFile1, buffer);
   foutdcr << buffer << '\n';
  }
  
  foutdcr << FlowBursts << '\n';
  foutdcr << FlowDeadline  << '\n';
  foutdcr << mtu  << '\n';
  foutdcr << rho  << '\n';
  iFile1.close();
  foutdcr.close();

  foutdmx << "p min " << NNodes << " " << NArcs << "\n";
  foutdmx << "n " << source << " 1\n";
  foutdmx << "n " << destination << " -1\n";
  
  for(i=0; i<NArcs; i++){
   iArc >> numberArc;
   iArc >> SN;
   iArc >> EN;
   iArc >> numberComm;
   iArc >> cost;
   iArc >> capacity;
   iArc >> numberArc;
   foutdmx << "a " << SN << " " << EN << " -1 " << capacity << " " << cost << "\n";
  }

  foutdmx.close();
  iArc.close();

  auto bsc = dynamic_cast< BlockSolverConfig * >(Configuration::deserialize( "MILPPar.txt" ) );
  if( ! bsc ) {
   cerr << "Error: configuration file not a BlockSolverConfig" << endl;
   exit( 1 );    
  }

  oMCFB = dynamic_cast< SingleFlowDCRBlock * >( Block::new_Block( "SingleFlowDCRBlock" ) );
  assert( oMCFB );

  load( "output.dmx", "output.dcr" );
  
  bsc->apply( oMCFB );
  bsc->clear();  // keep the clear()-ed BlockSolverConfig for final cleanup

  // check Solvers - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( oMCFB->get_registered_solvers().empty() ) {
    cerr << "Error: BlockSolverConfig did not register any Solver" << endl;
    exit( 1 );    
  }

  Solver * slvr = oMCFB->get_registered_solvers().front();
  int rtrn = slvr->compute( false );
  std::cout << "processing flow " << j << ": " << slvr->get_lb() << "\n";
  
  if( slvr->has_var_solution() ){
    index_feas += 1;
    output << j << "\n";
    output_sup << source << "\t" << index_feas << "\t" << rho << "\n";
    output_sup << destination << "\t" << index_feas << "\t" << rho1 << "\n";
    output_param << FlowBursts << "\t" << FlowDeadline <<  "\n";
  }
  
 }

  iFile2.close();
  iFile3.close();
  iNode.close();
  
 }  // end( main )

/*--------------------------------------------------------------------------*/
/*------------------------ End File test.cpp -------------------------------*/
/*--------------------------------------------------------------------------*/
