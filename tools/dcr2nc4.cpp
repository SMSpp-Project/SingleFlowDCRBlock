/*--------------------------------------------------------------------------*/
/*---------------------------- File dcr2nc4.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Small main() for constructing MultiFlowDCRBlock netCDF files out of the
 * multi-file textual format that MultiFlowDCRBlock::load() reads, i.e.,
 * <base>.nod (the number of flows, of nodes and of arcs), <base>.arc (the
 * arcs), <base>.sup (the source, the sink and the rate of each flow),
 * <base>.dcr (the delays of the nodes and of the arcs) and <base>.param (the
 * MTU, then the burst and the deadline of each flow).
 *
 * With -r the instance is rather in the raw format of the data sets it comes
 * from, in which <base>.dcr also holds, after the delays, the burst and the
 * deadline of each flow and, last, the MTU, and there is no <base>.param: the
 * two files that load() reads are written out of it first.
 *
 * With -s i the file is rather the SingleFlowDCRBlock of flow i (from 0) on
 * its own, as the single-flow instances of the data sets are made: each arc
 * with its own capacity and cost out of <base>.arc, not the capacity the
 * flows share and the unit cost of the sub-Block of a MultiFlowDCRBlock, and
 * the flow with its source, sink and rate out of <base>.sup, its burst and
 * deadline out of <base>.param and the MTU.
 *
 * With -k only the first k flows are kept, which is what makes an instance
 * with thousands of flows small enough for a formulation that holds them all
 * to be solved exactly: the flows are read in order, hence it is enough to
 * write k as their number in <base>.nod.
 *
 * All this happens in a temporary directory, which also takes the files that
 * load() writes in the current directory while building each sub-Block, so
 * that the directory the tool is run from, and the instance, are left as they
 * were.
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

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include <MultiFlowDCRBlock.h>

#include <SingleFlowDCRBlock.h>

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

namespace fs = std::filesystem;

/*--------------------------------------------------------------------------*/
/*------------------------------ FUNCTIONS ---------------------------------*/
/*--------------------------------------------------------------------------*/
/// Custom terminate function to print the exception message

static void smspp_terminate( void ) {
 std::cerr << "Uncaught exception in executing SMS++:\n";
 try {
  std::rethrow_exception( std::current_exception() );
 }
 catch( const std::exception & e ) {
  std::cerr << "\tException type: " << typeid( e ).name() << "\n";
  std::cerr << "\tException message: " << e.what() << "\n";
 } catch( ... ) {
  std::cerr << "\tUnknown exception" << std::endl;
 }
 std::abort();
 }

/*--------------------------------------------------------------------------*/
/// the whole content of \p name, which has to exist

static std::string read_file( const fs::path & name ) {
 std::ifstream in( name );
 if( ! in.is_open() ) {
  std::cerr << "Error: cannot open " << name << std::endl;
  exit( 1 );
  }
 std::stringstream ss;
 ss << in.rdbuf();
 return( ss.str() );
 }

/*--------------------------------------------------------------------------*/
/// writes the .dcr and the .param that load() reads out of a raw .dcr
/** The raw .dcr holds a first line, one line per node and per arc, one line
 * per flow with its burst and deadline, and then the MTU: the first part is
 * the .dcr that load() reads, and the MTU followed by the lines of the flows
 * is its .param. */

static void split_raw_dcr( const fs::path & raw , const fs::path & dcr ,
                           const fs::path & param , int nflows , int nnodes ,
                           int narcs ) {
 std::ifstream in( raw );
 if( ! in.is_open() ) {
  std::cerr << "Error: cannot open " << raw << std::endl;
  exit( 1 );
  }

 std::ofstream out_dcr( dcr );
 std::string line;
 for( int i = 0 ; i <= nnodes + narcs ; ++i ) {
  if( ! std::getline( in , line ) ) {
   std::cerr << "Error: " << raw << " ends before its delays" << std::endl;
   exit( 1 );
   }
  out_dcr << line << '\n';
  }

 std::vector< std::string > flows( nflows );
 for( auto & f : flows )
  if( ! std::getline( in , f ) ) {
   std::cerr << "Error: " << raw << " ends before its flows" << std::endl;
   exit( 1 );
   }

 double mtu;
 if( ! ( in >> mtu ) ) {
  std::cerr << "Error: " << raw << " has no MTU at its end" << std::endl;
  exit( 1 );
  }

 std::ofstream out_param( param );
 out_param << mtu << '\n';
 for( auto & f : flows )
  out_param << f << '\n';
 }

/*--------------------------------------------------------------------------*/
/// writes the SingleFlowDCRBlock of flow \p flow of the instance in \p tbase
/** \p tbase names the files that load() reads (the raw .dcr already split),
 * and the Block is built out of the two descriptions that a
 * SingleFlowDCRBlock loads, written in memory. */

static void write_single( const fs::path & tbase , int flow , int nflows ,
                          int nnodes , int narcs , const fs::path & out )
{
 auto name = [ & ]( const char * ext ) {
  return( fs::path( tbase.string() + ext ) );
  };

 // the flow: its source, sink and rate, its burst and deadline, and the MTU
 std::ifstream sup( name( ".sup" ) ) , param( name( ".param" ) );
 double src , dst , cmm , rho = 0 , rho1 , mtu , burst = 0 , deadline = 0;
 if( ! ( param >> mtu ) ) {
  std::cerr << "Error: " << name( ".param" ) << " has no MTU" << std::endl;
  exit( 1 );
  }
 for( int j = 0 ; j <= flow ; ++j )
  if( ! ( ( sup >> src >> cmm >> rho >> dst >> cmm >> rho1 ) &&
          ( param >> burst >> deadline ) ) ) {
   std::cerr << "Error: flow " << flow << " is not in " << tbase << std::endl;
   exit( 1 );
   }

 // the graph, each arc with its own capacity and cost
 std::stringstream dmx;
 dmx << "p min " << nnodes << " " << narcs << "\n"
     << "n " << src << " 1\n" << "n " << dst << " -1\n";
 std::ifstream arc( name( ".arc" ) );
 for( int i = 0 ; i < narcs ; ++i ) {
  double num , sn , en , nc , cost , cap , num2;
  if( ! ( arc >> num >> sn >> en >> nc >> cost >> cap >> num2 ) ) {
   std::cerr << "Error: " << name( ".arc" ) << " ends before its arcs"
             << std::endl;
   exit( 1 );
   }
  dmx << "a " << sn << " " << en << " -1 " << cap << " " << cost << "\n";
  }

 // the delays, then the burst, the deadline, the MTU and the rate
 // (a .dcr whose last line has no end of line would glue it to the burst)
 std::stringstream dcr;
 std::string delays = read_file( name( ".dcr" ) );
 if( ( ! delays.empty() ) && ( delays.back() != '\n' ) )
  delays += '\n';
 dcr << delays << burst << "\n" << deadline << "\n" << mtu << "\n" << rho << "\n";

 auto SFB = dynamic_cast< SingleFlowDCRBlock * >(
                              Block::new_Block( "SingleFlowDCRBlock" ) );
 SFB->load( dmx );
 SFB->load_dcr( dcr , SFB->get_NNodes() , SFB->get_NArcs() );

 netCDF::NcFile f( out.string() , netCDF::NcFile::replace );
 f.putAtt( "SMS++_file_type" , netCDF::NcInt() , eBlockFile );
 SFB->Block::serialize( f , eBlockFile );

 std::cout << out.filename().string() << ": flow " << flow << " of "
           << nflows << ", " << nnodes << " nodes, " << narcs << " arcs"
           << std::endl;
 delete SFB;
 }

/*--------------------------------------------------------------------------*/
/*--------------------------------- Main -----------------------------------*/
/*--------------------------------------------------------------------------*/

int main( int argc , char **argv )
{
 std::set_terminate( smspp_terminate );

 int keep = 0;       // flows to keep, 0 = all of them
 int single = -1;    // the flow written on its own, -1 = none
 bool raw = false;   // the .dcr is in the raw format
 std::vector< std::string > args;

 for( int i = 1 ; i < argc ; ++i ) {
  const std::string a = argv[ i ];
  if( ( a == "-k" ) && ( i + 1 < argc ) )
   keep = std::atoi( argv[ ++i ] );
  else
   if( ( a == "-s" ) && ( i + 1 < argc ) )
    single = std::atoi( argv[ ++i ] );
   else
    if( a == "-r" )
     raw = true;
    else
     args.push_back( a );
  }

 if( args.size() != 2 ) {
  std::cerr << "Usage: " << argv[ 0 ]
            << " [-r] [-k flows | -s flow] base out.nc4"
            << std::endl
            << "  base     the instance, i.e., base.nod, base.arc, base.sup,"
            << std::endl
            << "           base.dcr and base.param (not with -r)" << std::endl
            << "  -r       base.dcr is in the raw format, holding the flows "
            << "and the MTU" << std::endl
            << "  -k flows only the first flows are kept [all]" << std::endl
            << "  -s flow  the SingleFlowDCRBlock of that flow (from 0) "
            << "alone" << std::endl;
  return( 1 );
  }

 const fs::path base = fs::absolute( args[ 0 ] );
 const fs::path out = fs::absolute( args[ 1 ] );

 // the sizes, out of the first line of base.nod - - - - - - - - - - - - - -

 int nflows , nnodes , narcs;
 {
  std::istringstream nod( read_file( base.string() + ".nod" ) );
  if( ! ( nod >> nflows >> nnodes >> narcs ) ) {
   std::cerr << "Error: " << base.string() << ".nod does not start with the "
             << "number of flows, of nodes and of arcs" << std::endl;
   return( 1 );
   }
  }

 if( ( keep < 0 ) || ( keep > nflows ) ) {
  std::cerr << "Error: " << keep << " flows asked of an instance with "
            << nflows << std::endl;
  return( 1 );
  }

 // the files load() reads, in a temporary directory- - - - - - - - - - - - -

 // a directory of this run: create_directory() does not create a name that
 // is already there, so a random one is tried until one is new
 fs::path tmp;
 {
  std::random_device rd;
  std::error_code ec;
  int tries = 0;
  do
   tmp = fs::temp_directory_path() / ( "dcr2nc4." + std::to_string( rd() ) );
  while( ( ! fs::create_directory( tmp , ec ) ) && ( ++tries < 100 ) );
  if( tries == 100 ) {
   std::cerr << "Error: cannot create a temporary directory" << std::endl;
   return( 1 );
   }
  }
 const fs::path tbase = tmp / base.filename();
 auto in_tmp = [ & ]( const char * ext ) {
  return( fs::path( tbase.string() + ext ) );
  };

 for( auto ext : { ".arc" , ".sup" } )
  fs::copy_file( base.string() + ext , in_tmp( ext ) );

 if( raw )
  split_raw_dcr( base.string() + ".dcr" , in_tmp( ".dcr" ) ,
                 in_tmp( ".param" ) , nflows , nnodes , narcs );
 else
  for( auto ext : { ".dcr" , ".param" } )
   fs::copy_file( base.string() + ext , in_tmp( ext ) );

 {
  std::ofstream nod( in_tmp( ".nod" ) );
  nod << ( keep ? keep : nflows ) << "\t" << nnodes << "\t" << narcs << "\t"
      << narcs;
  }

 // one flow on its own, or load it there and write it- - - - - - - - - - - -

 if( single >= 0 ) {
  if( single >= nflows ) {
   std::cerr << "Error: flow " << single << " asked of an instance with "
             << nflows << std::endl;
   fs::remove_all( tmp );
   return( 1 );
   }
  write_single( tbase , single , nflows , nnodes , narcs , out );
  fs::remove_all( tmp );
  return( 0 );
  }

 const fs::path here = fs::current_path();
 fs::current_path( tmp );

 auto MFB = dynamic_cast< MultiFlowDCRBlock * >(
                              Block::new_Block( "MultiFlowDCRBlock" ) );
 MFB->load( in_tmp( ".nod" ).string() );

 fs::current_path( here );
 fs::remove_all( tmp );

 netCDF::NcFile f( out.string() , netCDF::NcFile::replace );
 f.putAtt( "SMS++_file_type" , netCDF::NcInt() , eBlockFile );
 MFB->Block::serialize( f , eBlockFile );

 std::cout << out.filename().string() << ": " << ( keep ? keep : nflows )
           << " flows of " << nflows << ", " << nnodes << " nodes, " << narcs
           << " arcs" << std::endl;

 delete MFB;
 return( 0 );

 }  // end( main )

/*--------------------------------------------------------------------------*/
/*------------------------- End File dcr2nc4.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
