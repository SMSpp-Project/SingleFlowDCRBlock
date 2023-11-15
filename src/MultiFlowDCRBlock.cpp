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
 * Copyright &copy by Antonio Frangioni
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
/*------------------------ OTHER INITIALIZATIONS ---------------------------*/
/*--------------------------------------------------------------------------*/

void  MultiFlowDCRBlock::load( std::istream & input , char frmt ){}

void MultiFlowDCRBlock::load( const std::string & input , char frmt )
{
 // ensure starting from clean slate
 guts_of_destructor();

 int i, j, source, destination, commodity;
 int SN, EN, numberArc, numberComm;
 float rho0, rho1, mtu, FlowBurstK, FlowDeadlineK;
 float capacity, cost, NC, NA, NN;

 std::ifstream iNode( input.substr(0,input.find_last_of('.'))+".nod" );
 if( ! iNode.is_open() )
  throw( std::invalid_argument( "can't open file .nod" ) );

 std::ifstream iFile2( input.substr(0,input.find_last_of('.'))+".sup" );
 if( ! iFile2.is_open() )
  throw( std::invalid_argument( "can't open file .sup" ) );

 std::ifstream iFile3( input.substr(0,input.find_last_of('.'))+".param" );
 if( ! iFile3.is_open() )
  throw( std::invalid_argument( "can't open file .param" ) );
 
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

 v_Block.resize( NComm );
 
 for(j=0; j<NC; ++j){
  //std::cout << j << "\n";

  std::ifstream iFile1( input.substr(0,input.find_last_of('.'))+".dcr" );
  if( ! iFile1.is_open() )
    throw( std::invalid_argument( "can't open file .dcr" ) );

  std::ifstream iArc( input.substr(0,input.find_last_of('.'))+".arc" );
  if( ! iArc.is_open() )
    throw( std::invalid_argument( "can't open file .arc" ) );

  iFile2 >> source;
  iFile2 >> commodity;
  iFile2 >> rho0;
  iFile2 >> destination;
  iFile2 >> commodity;
  iFile2 >> rho1;
  iFile3 >> FlowBurstK >> FlowDeadlineK;

  std::ofstream foutdmx("output.dmx");
  std::ofstream foutdcr("output.dcr");

  while (!iFile1.eof()) {
   for(i = 0; i < NNodes; i++){
    iFile1 >> NodeDelays[ i ];
    foutdcr << NodeDelays[ i ]  << '\n';
   }
   for(i = 0; i < NArcs; i++){
    iFile1 >> LinkDelays[ i ];
    foutdcr << LinkDelays[ i ]  << '\n';
   }
  }

  foutdcr << FlowBurstK << '\n';
  foutdcr << FlowDeadlineK  << '\n';
  foutdcr << mtu  << '\n';
  foutdcr << rho0  << '\n';

  iFile1.close();
  foutdcr.close();

  foutdmx << "p min " << NN << " " << NA << "\n";
  foutdmx << "n " << source << " 1\n";
  foutdmx << "n " << destination << " -1\n";

  C[ j ].resize( NArcs );
  rho[ j ].resize( 1 );
  MTU[ j ].resize( 1 );
  FlowBursts[ j ].resize( 1 );
  FlowDeadlines[ j ].resize( 1 );

  MTU[ j ][ 0 ] = mtu;
  FlowBursts[ j ][ 0 ] = FlowBurstK;
  FlowDeadlines[ j ][ 0 ] = FlowDeadlineK;
  rho[ j ][ 0 ] = rho0;

  for(i=0; i<NA; i++){
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
   foutdmx << "a " << SN << " " << EN << " -1 " << UTot[ i ] << " " << cost << "\n";
  }

  foutdmx.close();
  iArc.close();
  
  auto MCFB = new SingleFlowDCRBlock( this );
  std::ifstream fndmx( "output.dmx" );
  std::ifstream fndcr( "output.dcr" );
  MCFB->load( fndmx );
  Index NNodes = MCFB->get_NNodes();
  Index NArcs = MCFB->get_NArcs();

  try {  
   MCFB->load_dcr( fndcr , NNodes , NArcs );
   
  } 
  catch(...) {
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

 }  // end( MultiFlowDCRBlock::load( const std::string & input )

/*--------------------------------------------------------------------------*/

void MultiFlowDCRBlock::generate_abstract_variables( Configuration * stvv )
{/*
 // initialize the children - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ! ( AR & KnapsackRelaxation ) ) {
  v_Block.resize( NComm );

  for( Index k = 0 ; k < NComm ; ++k ) {
   auto MCFb = new SingleFlowDCRBlock( this );
   MCFb->load( NNodes , NArcs , Startn , Endn , UTot , C[ k ] , 0 , 0 , 0 , 0 , 
             NodeDelays , LinkDelays , FlowBursts[ k ] , FlowDeadlines[ k ] , MTU[ k ] , rho[ k ] ); 
   v_Block[ k ] = MCFb;
   }
  }

 // call the base class method to have it done in the sub-Block, if any
 //Block::generate_abstract_variables();

 AR |= HasVar;
 */
 }

/*--------------------------------------------------------------------------*/

void MultiFlowDCRBlock::generate_abstract_constraints( Configuration * stcc )
{

 // do it in the DCR/BKB respectively
 for( auto blck : v_Block )
  blck->generate_abstract_constraints();

 if( ! ( AR & KnapsackRelaxation ) ) {
  
  // initialize the vectors of coefficients, and reset count[]
  std::vector< LinearFunction::v_coeff_pair > coeffs( get_NArcs() );

  for( Index j = 0 ; j < get_NArcs() ; ++j ) {
   coeffs[ j ].resize( NComm );
   }

  for( Index k = 0 ; k < get_NComm() ; ++k )
   for( Index j = 0 ; j < get_NArcs() ; ++j )
    coeffs[ j ][ k ] = std::make_pair(
      static_cast< SingleFlowDCRBlock * >( v_Block[ k ] )->i2p_r( j ) , double( 1 ) );

  // generate the mutual capacity constraints  - - - - - - - - - - - - - - -
  // each constraint is an inequality, i.e., RHS = UTot[ j ]
  MCs.resize( get_NArcs() );
  for( Index j = 0 ; j < get_NArcs() ; ++j ) {
    LinearFunction::v_coeff_pair v_var;
    for( Index k = 0 ; k < get_NComm() ; ++k ) {
      v_var.push_back( coeffs[ j ][ k ] );
    }
    MCs[ j ].set_function( new LinearFunction( std::move( v_var )));
    MCs[ j ].set_rhs( UTot[j] ); 
    MCs[ j ].set_lhs( -Inf< double >() );
    }
  add_static_constraint( MCs , "Mut" );
  }

 AR |= HasMutual;

 }  // end( MultiFlowDCRBlock::generate_abstract_constraints() )

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

void MultiFlowDCRBlock::print( std::ostream & output , char vlvl ) const
{
 output << "MultiFlowDCRBlock with " << get_NComm() << " commodities, "
	<< get_NNodes() << " nodes and " << get_NArcs() << " arcs"
	<< std::endl;
 }

/*--------------------------------------------------------------------------*/

void MultiFlowDCRBlock::serialize( netCDF::NcGroup & group ) const
{ 
 // call the method of Block- - - - - - - - - - - - - - - - - - - - - - - - -

 Block::serialize( group );

 // now the MultiFlowDCRBlock data- - - - - - - - - - - - - - - - - - - - - - - - - -
 netCDF::NcDim nn = group.addDim( "NNodes" , get_NNodes() );
 netCDF::NcDim na = group.addDim( "NArcs" , get_NArcs() );
 netCDF::NcDim nc = group.addDim( "NComm" , get_NComm() );
 netCDF::NcDim ncnst = group.addDim( "NCnst" , NCnst);

 ( group.addVar( "SN" , netCDF::NcUint64() , na ) ).putVar( Startn.data() );

 ( group.addVar( "EN" , netCDF::NcUint64() , na ) ).putVar( Endn.data() );

 ( group.addVar( "Utot" , netCDF::NcDouble() , na ) ).putVar( UTot.data() );

 if( F.size() == NArcs )
  ( group.addVar( "F" , netCDF::NcDouble() , na ) ).putVar( F.data() );

 ::serialize( group, "U", netCDF::NcDouble(), U, {nc,na});
              
 ::serialize( group, "B", netCDF::NcDouble(), B, {nc,nn});
              
 ::serialize( group, "C", netCDF::NcDouble(), C, {nc,na});

 }  // end( MultiFlowDCRBlock::serialize )

/*--------------------------------------------------------------------------*/

void MultiFlowDCRBlock::deserialize( const netCDF::NcGroup & group )
{ 
 // erase previous instance, if any- - - - - - - - - - - - - - - - - - - - - -

 if( NNodes || NComm || get_NArcs() )
   MultiFlowDCRBlock();

 // read problem data- - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 auto nn = group.getDim( "NNodes" );
 if( nn.isNull() )
  throw( std::logic_error( "NNodes dimension is required" ) );
 NNodes = nn.getSize();

 auto na = group.getDim( "NArcs" );
 if( na.isNull() )
  throw( std::logic_error( "NArcs dimension is required" ) );
 NArcs = na.getSize();
 
 auto nc = group.getDim( "NComm" );
 if( nc.isNull() )
  throw( std::logic_error( "NComm dimension is required" ) );
 NComm = nc.getSize();
 
 Index NCnst = NArcs;
 auto ncnst = group.getDim( "NCnst" );
 if( nc.isNull() )
  throw( std::logic_error( "NCnst dimension is required" ) );
 NCnst = ncnst.getSize();
 
 auto sn = group.getVar( "SN" );
 if( sn.isNull() )
  throw( std::logic_error( "Starting Nodes not found" ) );

 Startn.resize( NArcs );
 sn.getVar( Startn.data() );

 auto en = group.getVar( "EN" );
 if( en.isNull() )
  throw( std::logic_error( "Ending Nodes not found" ) );

 Endn.resize( NArcs );
 en.getVar( Endn.data() );
 
 auto ut = group.getVar( "Utot" );
 if( ut.isNull() )
  throw( std::logic_error( "Total capacities not found" ) );

 UTot.resize( NArcs );
 ut.getVar( UTot.data() );

 auto fc = group.getVar( "F" );
 if( ! fc.isNull() ){
  F.resize( NArcs );
  fc.getVar( F.data() );
  }

 U.resize( NComm );
 for( int i = 0 ; i< NComm ; i++ )
  U[ i ].resize( NArcs );
 
 B.resize( NComm );
 for( int i = 0 ; i< NComm ; i++ )
  B[ i ].resize(NNodes); 
 
 C.resize( NComm );
 for( int i = 0 ; i< NComm ; i++)
  C[ i ].resize( NArcs ); 
 
 ::deserialize( group , "U" , U );
 ::deserialize( group , "B" , B );
 ::deserialize( group , "C" , C ); 

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
 {
  const auto sup = FCs.data() + FCs.num_elements();
  for( auto it = FCs.data() ; it != sup ; ++it )
   it->clear();
  }
 {
  const auto sup = SLCs.data() + SLCs.num_elements();
  for( auto it = SLCs.data() ; it != sup ; ++it )
   it->clear();
  }

 MCs.clear();
 FCs.resize( boost::extents[ 0 ][ 0 ] );
 SLCs.resize( boost::extents[ 0 ][ 0 ] );

 for( auto bk : v_Block )
  delete bk;

 v_Block.clear();

 NXtrV = NXtrC = 0;
 IdxBeg.clear();
 CoefIdx.clear();
 CoefVal.clear();
 C.clear();
 U.clear();
 B.clear();
 I.clear();

 UTot.clear();
 F.clear();

 Startn.clear();
 Endn.clear();

 NamesK.clear();
 Active.clear();
 ActiveK.clear();
 PT.clear();

 CIsCpy.clear();
 UIsCpy.clear();
 BIsCpy.clear();

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

 }  // end( guts_of_destructor )

/*--------------------------------------------------------------------------*/
/*---------------------- End File MultiFlowDCRBlock.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
