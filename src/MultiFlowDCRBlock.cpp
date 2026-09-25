/*--------------------------------------------------------------------------*/
/*------------------------- File MultiFlowDCRBlock.cpp ---------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the MultiFlowDCRBlock class.
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
 * \copyright &copy; by Antonio Frangioni, Laura Galli, Luca Mencarelli
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "MultiFlowDCRBlock.h"

//#include <math.h>

#include <ctype.h>

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*--------------------------------- TYPES ----------------------------------*/
/*--------------------------------------------------------------------------*/

using Index = Block::Index;

/*--------------------------------------------------------------------------*/
/*----------------------------- FUNCTIONS ----------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------------- STATIC MEMBERS -------------------------------*/
/*--------------------------------------------------------------------------*/

// register MultiFlowDCRBlock to the Block factory
SMSpp_insert_in_factory_cpp_1( MultiFlowDCRBlock );

/*--------------------------------------------------------------------------*/
/*-------------------------------- CONSTANTS -------------------------------*/
/*--------------------------------------------------------------------------*/

// oversubscription/discount factor used in load( const std::string & ) to
// compute the mutual (shared) arc capacity CapTot[] out of the "raw"
// individual arc capacity UTot[] and the number of commodities NC, as
// CapTot[ i ] = floor( FACTOR * NC * UTot[ i ] ): this reflects the fact
// that, typically, not all the commodities need their maximum individual
// capacity at the same time, so the arcs can be "oversubscribed" w.r.t.
// the sum of the individual capacities

static const auto FACTOR = 0.8;

/*--------------------------------------------------------------------------*/
/*------------------------ OTHER INITIALIZATIONS ---------------------------*/
/*--------------------------------------------------------------------------*/

// load the MultiFlowDCRBlock out of an istream: currently unimplemented,
// see the IMPLEMENTATION NOTE in MultiFlowDCRBlock.h

void MultiFlowDCRBlock::load( std::istream & input , char frmt ) {}

/*--------------------------------------------------------------------------*/

void MultiFlowDCRBlock::load( const std::string & input , char frmt )
{
 // ensure starting from clean slate
 guts_of_destructor();

 int i , j , source , destination , commodity;
 int SN , EN , numberArc , numberComm;
 float rho0 , rho1 , mtu , FlowBurstK , FlowDeadlineK;
 float capacity , cost , NC , NA , NN;

 // open the three "global" files: input + ".nod" (sizes and MTU),
 // input + ".sup" (per-commodity source/sink/rate) and input + ".param"
 // (per-commodity burst/deadline); \p frmt is not used, this is the only
 // supported format (see the IMPLEMENTATION NOTE in MultiFlowDCRBlock.h)

 std::ifstream iNode( input.substr( 0 , input.find_last_of( '.' ) ) + ".nod" );
 if( ! iNode.is_open() )
  throw( std::invalid_argument( "can't open file .nod" ) );

 std::ifstream iFile2( input.substr( 0 , input.find_last_of( '.' ) ) +
                       ".sup" );
 if( ! iFile2.is_open() )
  throw( std::invalid_argument( "can't open file .sup" ) );

 std::ifstream iFile3( input.substr( 0 , input.find_last_of( '.' ) ) +
                       ".param" );
 if( ! iFile3.is_open() )
  throw( std::invalid_argument( "can't open file .param" ) );

 // read the global sizes (number of commodities, nodes, arcs) and the MTU
 // (shared by all commodities) - - - - - - - - - - - - - - - - - - - - - -

 iNode >> NC;
 iNode >> NN;
 iNode >> NA;

 NNodes = NN;
 NArcs = NA;
 NComm = NC;
 iFile3 >> mtu;

 NComm = NC;
 NArcs = NA;
 NNodes = NN;
 C.resize( NComm );
 MTU.resize( NComm );
 rho.resize( NComm );
 FlowBursts.resize( NComm );
 FlowDeadlines.resize( NComm );

 NodeDelays.resize( NNodes );
 LinkDelays.resize( NArcs );

 Startn.resize( NArcs );
 Endn.resize( NArcs );
 UTot.resize( NArcs );
 CapTot.resize( NArcs );

 v_Block.resize( NComm );

 // for each commodity, read its own topology/delay/traffic data and build
 // the corresponding SingleFlowDCRBlock sub-Block - - - - - - - - - - - - -

 for( j = 0 ; j < NC ; ++j ) {
  //std::cout << j << "\n";

  C[ j ].resize( NArcs );

  // open the per-commodity delay-related file (".dcr") and the (shared,
  // but re-read for every commodity) arc topology/individual-capacity
  // file (".arc")

  std::ifstream iFile1( input.substr( 0 , input.find_last_of( '.' ) ) +
                        ".dcr" );
  if( ! iFile1.is_open() )
   throw( std::invalid_argument( "can't open file .dcr" ) );

  std::ifstream iArc( input.substr( 0 , input.find_last_of( '.' ) ) + ".arc" );
  if( ! iArc.is_open() )
   throw( std::invalid_argument( "can't open file .arc" ) );

  // read the source/sink node and rate of this commodity from ".sup", and
  // its traffic burst/deadline from ".param"

  iFile2 >> source;
  iFile2 >> commodity;
  iFile2 >> rho0;
  iFile2 >> destination;
  iFile2 >> commodity;
  iFile2 >> rho1;
  iFile3 >> FlowBurstK >> FlowDeadlineK;

  // build, on the fly, a temporary DIMACS-like description of this
  // commodity's DCR sub-problem: "output.dmx" (topology, mutual capacity
  // as arc upper bound, and unit arc cost) and "output.dcr" (the delay
  // data copied verbatim from ".dcr", plus burst/deadline/mtu/rate
  // appended at the end)

  std::ofstream foutdmx( "output.dmx" );
  std::ofstream foutdcr( "output.dcr" );

  while( ! iFile1.eof() ) {
   std::string buffer;
   getline( iFile1 , buffer );
   foutdcr << buffer << '\n';
   }

  foutdcr << FlowBurstK << '\n';
  foutdcr << FlowDeadlineK << '\n';
  foutdcr << mtu << '\n';
  foutdcr << rho0 << '\n';

  iFile1.close();
  foutdcr.close();

  foutdmx << "p min " << NN << " " << NA << "\n";
  foutdmx << "n " << source << " 1\n";
  foutdmx << "n " << destination << " -1\n";

  MTU[ j ] = mtu;
  FlowBursts[ j ] = FlowBurstK;
  FlowDeadlines[ j ] = FlowDeadlineK;
  rho[ j ] = rho0;

  // read the arc topology and individual capacity from ".arc"; the arc
  // cost is hard-wired to 1 (the actual cost read from file is ignored);
  // the mutual capacity CapTot[ i ] is computed here (see FACTOR above)
  // and used as the arc upper bound in the temporary DIMACS file, since
  // the per-commodity sub-Block only sees the shared capacity

  for( i = 0 ; i < NA ; i++ ) {
   iArc >> numberArc;
   iArc >> SN;
   Startn[ i ] = SN;
   iArc >> EN;
   Endn[ i ] = EN;
   iArc >> numberComm;
   iArc >> cost;
   C[ j ][ i ] = 1;
   iArc >> UTot[ i ];
   iArc >> numberArc;
   //UTot[ i ] = (j+1)*UTot[ i ];
   CapTot[ i ] = std::floor( FACTOR * NC * UTot[ i ] );
   foutdmx << "a " << Startn[ i ] << " " << Endn[ i ] << " -1 " << CapTot[ i ]
           << " " << C[ j ][ i ] << "\n"; //0.5*NC*UTot[ i ]
   }

  foutdmx.close();
  iArc.close();

  // create the SingleFlowDCRBlock sub-Block for this commodity (owned by
  // this MultiFlowDCRBlock, i.e., "this" is set as its father) and load
  // it out of the temporary DIMACS/DCR files just written

  auto MCFB = new SingleFlowDCRBlock( this );
  std::ifstream fndmx( "output.dmx" );
  std::ifstream fndcr( "output.dcr" );
  MCFB->load( fndmx );
  Index NNodes = MCFB->get_NNodes();
  Index NArcs = MCFB->get_NArcs();

  try {
   MCFB->load_dcr( fndcr , NNodes , NArcs );

  } catch( ... ) {
   std::cerr << "Error: dcr file error!" << std::endl;
   return;
   }
  v_Block[ j ] = MCFB;
  }

 iFile2.close();
 iFile3.close();
 iNode.close();

 // issue Modification- - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // note: this is a NBModification, the "nuclear option"

 if( anyone_there() )
  add_Modification( std::make_shared< NBModification >( this ) );

 } // end( MultiFlowDCRBlock::load( const std::string & input )

/*--------------------------------------------------------------------------*/
// NOTE: this unconditionally uses the "flow" formulation (each commodity
// is a SingleFlowDCRBlock sub-Block); \p stvv is currently ignored, so the
// alternative "knapsack" formulation documented in MultiFlowDCRBlock.h is
// not actually selectable

void MultiFlowDCRBlock::generate_abstract_variables( Configuration * stvv )
{
 for( auto blck : v_Block )
  blck->generate_abstract_variables();
 }

/*--------------------------------------------------------------------------*/
// first lets each commodity sub-Block construct its own (per-commodity)
// abstract Constraint, then, unless already done (AR == 0, i.e., HasMutual
// not yet set), builds the mutual capacity Constraint MCs: one
// FRowConstraint per arc j, coupling the reserved-rate Variable r^k[ j ] of
// all the NComm commodities via coefficient 1, with RHS = CapTot[ j ] and
// LHS = -Inf (i.e., an inequality \sum_k r^k[ j ] <= CapTot[ j ])

void MultiFlowDCRBlock::generate_abstract_constraints( Configuration * stcc )
{

 // do it in the DCR/BKB respectively
 for( auto blck : v_Block )
  blck->generate_abstract_constraints();

 if( ! ( AR ) ) {

  // initialize the vectors of coefficients, and reset count[]
  std::vector< LinearFunction::v_coeff_pair > coeffs( get_NArcs() );

  for( Index j = 0 ; j < get_NArcs() ; ++j ) {
   coeffs[ j ].resize( NComm );
   }

  // for each commodity k and arc j, the coefficient of the mutual
  // capacity constraint on r^k[ j ] is 1 (a pointer to the ColVariable is
  // obtained from the k-th commodity's SingleFlowDCRBlock sub-Block)

  for( Index k = 0 ; k < get_NComm() ; ++k )
   for( Index j = 0 ; j < get_NArcs() ; ++j )
    coeffs[ j ][ k ] = std::make_pair(
     static_cast< SingleFlowDCRBlock * >( v_Block[ k ] )->i2p_r( j ) ,
     double( 1 ) );

  // generate the mutual capacity constraints  - - - - - - - - - - - - - - -
  // each constraint is an inequality, i.e., RHS = CapTot[ j ]
  MCs.resize( get_NArcs() );
  for( Index j = 0 ; j < get_NArcs() ; ++j ) {
   LinearFunction::v_coeff_pair v_var;
   for( Index k = 0 ; k < get_NComm() ; ++k ) {
    v_var.push_back( coeffs[ j ][ k ] );
    }
   MCs[ j ].set_function( new LinearFunction( std::move( v_var ) ) );
   MCs[ j ].set_rhs( CapTot[ j ] ); //0.5*get_NComm()*UTot[ j ]
   MCs[ j ].set_lhs( -Inf< double >() );
   }
  add_static_constraint( MCs , "Mut" );
  }

 AR |= HasMutual;

 } // end( MultiFlowDCRBlock::generate_abstract_constraints() )

/*--------------------------------------------------------------------------*/
// checks the mutual capacity Constraint MCs, which are the only ones this
// Block has of its own, and then hands the rest over to the sub-Block, each
// of which checks the DCR problem of its own commodity; \p useabstract is
// unused here, as RowConstraint::is_feasible() always works on the abstract
// Constraint, but it is passed down to the sub-Block

bool MultiFlowDCRBlock::is_feasible( bool useabstract , Configuration * fsbc )
{
 double tol = 0;
 bool rel_viol = true;

 // try to extract, from "c", the parameters that determine feasibility; if
 // it succeeds, it sets the values of the parameters and returns true,
 // otherwise it returns false

 auto extract_parameters = [ &tol , &rel_viol ]( Configuration * c ) -> bool {
  if( auto tc = dynamic_cast< SimpleConfiguration< double > * >( c ) ) {
   tol = tc->f_value;
   return( true );
   }
  if( auto tc =
       dynamic_cast< SimpleConfiguration< std::pair< double , int > > * >(
        c ) ) {
   tol = tc->f_value.first;
   rel_viol = tc->f_value.second;
   return( true );
   }
  return( false );
  };

 if( ( ! extract_parameters( fsbc ) ) && f_BlockConfig )
  // if the given Configuration is not valid, try the one from the BlockConfig
  extract_parameters( f_BlockConfig->f_is_feasible_Configuration );

 if( ! RowConstraint::is_feasible( MCs , tol , rel_viol ) )
  return( false );

 // the commodities: Block::is_feasible() checks every sub-Block
 return( Block::is_feasible( useabstract , fsbc ) );

 } // end( MultiFlowDCRBlock::is_feasible )

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/
// currently a no-op: the actual printing code is commented out, so nothing
// is written to output regardless of vlvl

void MultiFlowDCRBlock::print( std::ostream & output , char vlvl ) const
{
 //output << "MultiFlowDCRBlock with " << get_NComm() << " commodities, "
 //     << get_NNodes() << " nodes and " << get_NArcs() << " arcs"
 //     << std::endl;
 }

/*--------------------------------------------------------------------------*/
// writes the data the flows share (topology, individual and mutual
// capacities), the data of each flow and the SingleFlowDCRBlock of each flow,
// each in a group of its own, so that deserialize() rebuilds the whole Block

void MultiFlowDCRBlock::serialize( netCDF::NcGroup & group ) const
{
 // call the method of Block- - - - - - - - - - - - - - - - - - - - - - - - -

 Block::serialize( group );

 // the data the flows share- - - - - - - - - - - - - - - - - - - - - - - - -

 netCDF::NcDim nn = group.addDim( "NNodes" , get_NNodes() );
 netCDF::NcDim na = group.addDim( "NArcs" , get_NArcs() );
 netCDF::NcDim nc = group.addDim( "NComm" , get_NComm() );

 ( group.addVar( "SN" , netCDF::NcUint64() , na ) ).putVar( Startn.data() );

 ( group.addVar( "EN" , netCDF::NcUint64() , na ) ).putVar( Endn.data() );

 ( group.addVar( "Utot" , netCDF::NcDouble() , na ) ).putVar( UTot.data() );

 ( group.addVar( "CapTot" , netCDF::NcDouble() , na ) ).putVar(
                                                           CapTot.data() );

 // the data of each flow - - - - - - - - - - - - - - - - - - - - - - - - - -

 ( group.addVar( "MTU" , netCDF::NcDouble() , nc ) ).putVar( MTU.data() );

 ( group.addVar( "rho" , netCDF::NcDouble() , nc ) ).putVar( rho.data() );

 ( group.addVar( "FlowBursts" , netCDF::NcDouble() , nc ) ).putVar(
                                                       FlowBursts.data() );

 ( group.addVar( "FlowDeadlines" , netCDF::NcDouble() , nc ) ).putVar(
                                                    FlowDeadlines.data() );

 // and the SingleFlowDCRBlock of each flow, which holds all the rest- - - - -

 for( Index k = 0 ; k < get_NComm() ; ++k ) {
  auto sub = group.addGroup( "SingleFlowDCRBlock_" + std::to_string( k ) );
  v_Block[ k ]->serialize( sub );
  }

 } // end( MultiFlowDCRBlock::serialize )

/*--------------------------------------------------------------------------*/
// reads what serialize() writes (see the description of the netCDF format in
// MultiFlowDCRBlock::deserialize() in MultiFlowDCRBlock.h), the
// SingleFlowDCRBlock of each flow included; the matrices "U", "B" and "C"
// that files of an earlier format may hold are read when there

void MultiFlowDCRBlock::deserialize( const netCDF::NcGroup & group )
{
 // erase previous instance, if any- - - - - - - - - - - - - - - - - - - - - -

 guts_of_destructor();

 // read the data the flows share - - - - - - - - - - - - - - - - - - - - - -

 auto nn = group.getDim( "NNodes" );
 if( nn.isNull() )
  throw( std::invalid_argument( "MultiFlowDCRBlock::deserialize: NNodes "
                                "dimension is required" ) );
 NNodes = nn.getSize();

 auto na = group.getDim( "NArcs" );
 if( na.isNull() )
  throw( std::invalid_argument( "MultiFlowDCRBlock::deserialize: NArcs "
                                "dimension is required" ) );
 NArcs = na.getSize();

 auto nc = group.getDim( "NComm" );
 if( nc.isNull() )
  throw( std::invalid_argument( "MultiFlowDCRBlock::deserialize: NComm "
                                "dimension is required" ) );
 NComm = nc.getSize();
 NCnst = NArcs;

 auto sn = group.getVar( "SN" );
 if( sn.isNull() )
  throw( std::invalid_argument( "MultiFlowDCRBlock::deserialize: SN "
                                "(the starting nodes) is required" ) );
 Startn.resize( NArcs );
 sn.getVar( Startn.data() );

 auto en = group.getVar( "EN" );
 if( en.isNull() )
  throw( std::invalid_argument( "MultiFlowDCRBlock::deserialize: EN "
                                "(the ending nodes) is required" ) );
 Endn.resize( NArcs );
 en.getVar( Endn.data() );

 auto ut = group.getVar( "Utot" );
 if( ut.isNull() )
  throw( std::invalid_argument( "MultiFlowDCRBlock::deserialize: Utot "
                                "(the total capacities) is required" ) );
 UTot.resize( NArcs );
 ut.getVar( UTot.data() );

 // the mutual capacities, which a file without them has computed as load()
 // computes them
 CapTot.resize( NArcs );
 auto ct = group.getVar( "CapTot" );
 if( ct.isNull() )
  for( Index i = 0 ; i < NArcs ; ++i )
   CapTot[ i ] = std::floor( FACTOR * NComm * UTot[ i ] );
 else
  ct.getVar( CapTot.data() );

 // the data of each flow, which the SingleFlowDCRBlock hold as well- - - - -

 auto read_flows = [ & ]( const char * name , Vec_double & v ) {
  v.assign( NComm , 0 );
  auto var = group.getVar( name );
  if( ! var.isNull() )
   var.getVar( v.data() );
  };

 read_flows( "MTU" , MTU );
 read_flows( "rho" , rho );
 read_flows( "FlowBursts" , FlowBursts );
 read_flows( "FlowDeadlines" , FlowDeadlines );

 NodeDelays.resize( NNodes );
 LinkDelays.resize( NArcs );

 // the matrices of an earlier format, if there - - - - - - - - - - - - - - -

 if( ! group.getVar( "U" ).isNull() ) {
  U.assign( NComm , Vec_double( NArcs ) );
  ::deserialize( group , "U" , U );
  }
 if( ! group.getVar( "B" ).isNull() ) {
  B.assign( NComm , Vec_double( NNodes ) );
  ::deserialize( group , "B" , B );
  }
 if( ! group.getVar( "C" ).isNull() ) {
  C.assign( NComm , Vec_double( NArcs ) );
  ::deserialize( group , "C" , C );
  }

 // the SingleFlowDCRBlock of each flow - - - - - - - - - - - - - - - - - - -

 v_Block.assign( NComm , nullptr );
 for( Index k = 0 ; k < NComm ; ++k ) {
  const std::string name = "SingleFlowDCRBlock_" + std::to_string( k );
  auto sub = group.getGroup( name );
  if( sub.isNull() )
   throw( std::invalid_argument( "MultiFlowDCRBlock::deserialize: " + name +
                                 " not present" ) );
  v_Block[ k ] = new_Block( sub , this );
  if( ! dynamic_cast< SingleFlowDCRBlock * >( v_Block[ k ] ) )
   throw( std::invalid_argument( "MultiFlowDCRBlock::deserialize: " + name +
                                 " is not a SingleFlowDCRBlock" ) );
  }

 // call the method of Block- - - - - - - - - - - - - - - - - - - - - - - - -
 // inside this the NBModification, the "nuclear option",  is issued

 Block::deserialize( group );

 } // end( MultiFlowDCRBlock::deserialize )

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

void MultiFlowDCRBlock::guts_of_destructor( void )
{
 /* clear() all Constraint to ensure that they do not bother to un-register
    themselves from Variable that are going to be deleted anyway. Then
    deletes all the "abstract representation", if any. */

 for( auto & cnst : MCs )
  cnst.clear();

 MCs.clear();

 for( auto bk : v_Block )
  delete bk;

 v_Block.clear();

 C.clear();
 U.clear();
 B.clear();
 I.clear();

 UTot.clear();
 CapTot.clear();

 Startn.clear();
 Endn.clear();

 NodeDelays.clear();
 FlowDeadlines.clear();
 FlowBursts.clear();
 FlowDeadlines.clear();
 rho.clear();
 MTU.clear();

 StrtNme = 0;
 NNodes = 0;
 NArcs = 0;
 NComm = 0;
 NCnst = 0;

 Constraint::clear( MCs );

 reset_static_constraints();
 reset_static_variables();
 reset_dynamic_constraints();
 reset_dynamic_variables();
 reset_objective();

 // explicitly reset all Constraint and Variable
 // this is done for the case where this method is called prior to re-loading
 // a new instance: if not, the new representation would be added to the
 // (no longer current) one
 reset_static_constraints();
 // not needed, there isn't any - reset_static_variables();
 // not needed, there isn't any - reset_dynamic_constraints();
 // not needed, there isn't any - reset_dynamic_variables();
 // not needed, there isn't any - reset_objective();

 AR = 0;

 } // end( guts_of_destructor )

/*--------------------------------------------------------------------------*/
/*--------------------- End File MultiFlowDCRBlock.cpp ---------------------*/
/*--------------------------------------------------------------------------*/
