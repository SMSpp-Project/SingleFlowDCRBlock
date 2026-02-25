/*--------------------------------------------------------------------------*/
/*--------------------- File SingleFlowDCRBlock.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the SingleFlowDCRBlock class.
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
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "SingleFlowDCRBlock.h"
#include <iomanip>

/*--------------------------------------------------------------------------*/
/*--------------------------------- MACROS ---------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef NDEBUG
 #define CHECK_DS 0
 /* Perform long and costly checks on the data structures representing the
  * abstract and the physical representations agree. */
#else
 #define CHECK_DS 0
 // never change this
#endif

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*-------------------------------- TYPES -----------------------------------*/
/*--------------------------------------------------------------------------*/

using Index = Block::Index;
using c_Index = Block::c_Index;

using Range = Block::Range;
using c_Range = Block::c_Range;

using Subset = Block::Subset;
using c_Subset = Block::c_Subset;

using FNumber = SingleFlowDCRBlock::FNumber;

/*--------------------------------------------------------------------------*/
/*-------------------------------- CONSTANTS -------------------------------*/
/*--------------------------------------------------------------------------*/

static constexpr auto dNAN = std::numeric_limits< double >::quiet_NaN();
static const auto cuts_formulation = false;

static constexpr unsigned char FormMsk = 3;
// mask for the first three bits, i.e., the formulation

static constexpr unsigned char PCuts = 1;
/// the P/C formulation is used

static constexpr unsigned char SOCP = 2;
/// the SOCP formulation is used

/*--------------------------------------------------------------------------*/
/*-------------------------------- FUNCTIONS -------------------------------*/
/*--------------------------------------------------------------------------*/

// in the DIMACS format, comment lines start with 'c'

static std::istream & eatDMXcomments( std::istream& is )
{
 for(;;) {
  is >> std::ws;  // skip whitespaces
  if( is.peek() == is.widen( 'c' ) )
   // a comment: skip the rest of line and move to next
   is.ignore( std::numeric_limits< std::streamsize >::max() , is.widen( '\n' ) );
  else
   break;
  }

 return( is );
 }

/*--------------------------------------------------------------------------*/

static void print_UB( std::ostream & os , FNumber ub )
{
 if( ub == Inf< FNumber >() )
  os << "+Inf";
 else
  os << ub;
 }

/*--------------------------------------------------------------------------*/

static FNumber read_UB( std::istream & iStrm )
{
 iStrm >> eatcomments;
 int c = iStrm.peek();
 if( ! iStrm )
  throw( std::invalid_argument( "error reading the input stream" ) );
  
 if( ( c != 'I' ) && ( c != 'i' ) ) {
  FNumber res;
  iStrm >> res;
  if( ! iStrm )
   throw( std::invalid_argument( "error reading the input stream" ) );
  return( res );
  }

 do { c = iStrm.get(); c = iStrm.peek();
      if( ! iStrm )
       throw( std::invalid_argument( "error reading the input stream" ) );

  } while( ( c != iStrm.widen( ' ' ) ) &&
	   ( c != iStrm.widen( '\n' ) ) &&
	   ( c != iStrm.widen( '\t' ) ) );

 return( Inf< FNumber >() );
 }

/*--------------------------------------------------------------------------*/
// returns the number of elements where two vectors differ

template< typename T >
static Index countdiff( T beg , T end , T cmp )
{
 Index ndiff = 0;
 for( ; beg != end ; )
  if( *(beg++) != *(cmp++) )
   ndiff++;

 return( ndiff );
 }

/*--------------------------------------------------------------------------*/
// returns true if two vectors differ, one of them being given as a base
// vector and a subset of indices

template< typename T >
static bool is_equal( std::vector< T > & vec , c_Subset & nms ,
		      typename std::vector< T >::const_iterator cmp ,
		      Index n_max )
{
 for( auto nm : nms ) {
  if( nm >= n_max )
   throw( std::invalid_argument( "invalid name in nms" ) );
  if( vec[ nm ] != *(cmp++) )
   return( false );
  }

 return( true );
 }

/*--------------------------------------------------------------------------*/
// returns the number of elements where two vectors differ, one of them
// being given as a base vector and a subset of indices

template< typename T >
static Index countdiff( std::vector< T > & vec , c_Subset & nms ,
			typename std::vector< T >::const_iterator cmp ,
			Index n_max )
{
 Index ndiff = 0;
 for( auto nm : nms ) {
  if( nm >= n_max )
   throw( std::invalid_argument( "invalid name in nms" ) );
  if( vec[ nm ] != *(cmp++) )
   ndiff++;
  }

 return( ndiff );
 }

/*--------------------------------------------------------------------------*/
// copys one vector to a given subset of another

template< typename T >
static void copyidx( std::vector< T > & vec , c_Subset & nms ,
		     typename std::vector< T >::const_iterator cpy )
{
 for( auto nm : nms )
  vec[ nm ] = *(cpy++);
 }

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register SingleFlowDCRBlock to the Block factory

SMSpp_insert_in_factory_cpp_1( SingleFlowDCRBlock );

// register DCRSolution to the Solution factory

SMSpp_insert_in_factory_cpp_0( DCRSolution );

/*--------------------------------------------------------------------------*/
/*--------------------------- METHODS OF SingleFlowDCRBlock --------------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void SingleFlowDCRBlock::load( Index n , Index m , c_Subset & pEn , c_Subset & pSn ,
		     c_Vec_FNumber & pU , c_Vec_CNumber & pC ,
         c_Vec_FNumber & pNodeDelays , c_Vec_FNumber & pLinkDelays , 
         c_FNumber pFlowBursts , c_FNumber pFlowDeadlines , c_FNumber pMTU , c_FNumber prho )
{
 // sanity checks - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( pSn.size() < m )
  throw( std::invalid_argument( "pSn too small" ) );

 if( pEn.size() < m )
  throw( std::invalid_argument( "pEn too small" ) );

 if( ( pC.size() > 0 ) && ( pC.size() < m ) )
  throw( std::invalid_argument( "pC nonempty but too small" ) );

 if( ( pU.size() > 0 ) && ( pU.size() < m ) )
  throw( std::invalid_argument( "pU nonempty but too small" ) );

/*
 if( ( pB.size() > 0 ) && ( pB.size() < n ) )
  throw( std::invalid_argument( "pB nonempty but too small" ) );
*/

 if( ( pNodeDelays.size() > 0 ) && ( pNodeDelays.size() < n ) )
  throw( std::invalid_argument( "pNodeDelays nonempty but too small" ) );

 if( ( pLinkDelays.size() > 0 ) && ( pLinkDelays.size() < m ) )
  throw( std::invalid_argument( "pLinkDelays nonempty but too small" ) );

 // erase previous instance, if any- - - - - - - - - - - - - - - - - - - - - -

 if( MaxNNodes || get_MaxNArcs() )
  guts_of_destructor();
		   
 // copy over problem data - - - - - - - - - - - - - - - - - - - - - - - - - -
 
 NNodes = n;
 NArcs = m;
 MaxNNodes = NNodes;
 c_Index MaxNArcs = NArcs;

 SN.resize( MaxNArcs );
 if( std::any_of( pSn.begin() , pSn.begin() + m ,
		  [ n ]( c_Index sn ) { return( ( sn < 1 ) || ( sn > n ) ); }
		  ) )
  throw( std::invalid_argument( "wrong starting node" ) );
 std::copy( pSn.begin() , pSn.begin() + m , SN.begin() );

 EN.resize( MaxNArcs );
 if( std::any_of( pEn.begin() , pEn.begin() + m ,
		  [ n ]( c_Index en ) { return( ( en < 1 ) || ( en > n ) ); }
		  ) )
  throw( std::invalid_argument( "wrong ending node" ) );
 std::copy( pEn.begin() , pEn.begin() + m , EN.begin() );

 if( pC.empty() )
  C.assign( MaxNArcs , 0 );
 else {
  C.resize( MaxNArcs );
  std::copy( pC.begin() , pC.begin() + m , C.begin() );
  }

 if( std::any_of( pU.begin() , pU.begin() + m ,
		  []( c_FNumber ui ) { return( ui < Inf< FNumber >() ); } ) ) {
  U.resize( MaxNArcs );
  std::copy( pU.begin() , pU.begin() + m , U.begin() );
  }
 else
  U.clear();

 if( std::any_of( pNodeDelays.begin() , pNodeDelays.begin() + n ,
		  []( c_FNumber NodeDelaysi ) { return( NodeDelaysi < Inf< FNumber >() ); } ) ) {
  NodeDelays.resize( MaxNNodes );
  std::copy( pNodeDelays.begin() , pNodeDelays.begin() + n , NodeDelays.begin() );
  }
 else
  NodeDelays.clear();

 if( std::any_of( pLinkDelays.begin() , pLinkDelays.begin() + m ,
		  []( c_FNumber LinkDelaysi ) { return( LinkDelaysi < Inf< FNumber >() ); } ) ) {
  LinkDelays.resize( MaxNArcs );
  std::copy( pLinkDelays.begin() , pLinkDelays.begin() + m , LinkDelays.begin() );
  }
 else
  LinkDelays.clear();

 FlowBursts = pFlowBursts;
 FlowDeadlines = pFlowDeadlines;
 MTU = pMTU;
 rho = prho;

 // allocate flow variables - - - - - - - - - - - - - - - - - - - - - - - - -

 generate_abstract_variables();

 // the arcs whose cost is infinite have to be closed
 // in addition the cost has to be set to 0 - - - - - - - - - - - - - - - - -

 for( Index j = 0 ; j < C.size() ; ++j )
  if( C[ j ] >= Inf< CNumber >() ) {
   close_arc( j , eNoMod , eNoMod );
   C[ j ] = 0;
   }

 // throw Modification- - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // note: this is a NBModification, the "nuclear option"

 if( anyone_there() )
  add_Modification( std::make_shared< NBModification >( this ) );

 }  // end( SingleFlowDCRBlock::load( memory ) )

/*--------------------------------------------------------------------------*/

void SingleFlowDCRBlock::load( std::istream & input , char frmt )
{
 // erase previous instance, if any- - - - - - - - - - - - - - - - - - - - - -

 if( MaxNNodes || get_MaxNArcs() )
  guts_of_destructor();

 // read first non-comment line - - - - - - - - - - - - - - - - - - - - - - -

 char c;
 if( ! ( input >> eatDMXcomments >> c ) )
  throw( std::invalid_argument( "error reading the input stream" ) );

 if( c != 'p' )
  throw( std::invalid_argument( "format error in the input stream" ) );

 input >> eatDMXcomments;
 input.ignore( 3 , ' ' );  // skip "min"

 if( ! ( input >> eatDMXcomments >> NNodes ) )
  throw( std::invalid_argument( "LoadDMX: error reading number of nodes" ) );

 Index tm;
 if( ! ( input >> eatDMXcomments >> NArcs ) )
  throw( std::invalid_argument( "LoadDMX: error reading number of arcs" ) );

 // allocate memory - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 SN.resize( NArcs );
 EN.resize( NArcs );
 C.assign( NArcs , 0 );
 U.assign( NArcs , Inf< FNumber >() );
 B.assign( NNodes , 0 );

 // read problem data - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Index i = 0;  // arc counter
 int jdeficit = 0; // deficit counter
 for(;;) {
  if( ! ( input >> eatDMXcomments >> c ) )  // read next descriptor
   break;                                   // if none, end

  switch( c ) {
   case( 'n' ):  // description of a node
    Index j;
    if( ! ( input >> j ) )
     throw( std::invalid_argument( "error reading node name" ) );

    if( ( j < 1 ) || ( j > NNodes ) )
     throw( std::invalid_argument( "invalid node name" ) );
    
    FNumber Dfctj;
    if( ! ( input >> Dfctj ) )
     throw( std::invalid_argument( "error reading deficit" ) );

    B[ j - 1 ] = -Dfctj;

    if( jdeficit > 1 )
      throw( std::invalid_argument( "too many deficits" ) );
    
    jdeficit++;
    break;

   case( 'a' ):  // description of an arc
    if( i == NArcs )
     throw( std::invalid_argument( "too many arc descriptors" ) );

    if( ! ( input >> SN[ i ] ) )
     throw( std::invalid_argument( "error reading start node" ) );

    if( ( SN[ i ] < 1 ) || ( SN[ i ] > NNodes ) )
     throw( std::invalid_argument( "invalid start node" ) );

    if( ! ( input >> EN[ i ] ) )
     throw( std::invalid_argument( "error reading end node" ) );

    if( ( EN[ i ] < 1 ) || ( EN[ i ] > NNodes ) )
     throw( std::invalid_argument( "LoadDMX: invalid end node" ) );

    if( SN[ i ] == EN[ i ] )
     throw( std::invalid_argument( "self-loops not permitted" ) );

    FNumber LB;
    if( ! ( input >> LB ) )
     throw( std::invalid_argument( "error reading lower bound" ) );

    U[ i ] = read_UB( input );

    if( ! ( input >> C[ i ] ) )
     throw( std::invalid_argument( "error reading arc cost" ) );

    if( U[ i ] < LB )
     throw( std::invalid_argument( "lower bound > upper bound" ) );

    i++;
    break; 

   default:  // invalid code- - - - - - - - - - - - - - - - - - - - - - - - -
    throw( std::invalid_argument( "invalid DMX code" ) );

   }  // end( switch( c ) )
  }  // end( for( ever ) )

 if( i < NArcs )
  throw( std::invalid_argument( "too few arc descriptors" ) );

 f_cond_lower = dNAN;  // reset conditional bounds

 // simplify out the deta structures- - - - - - - - - - - - - - - - - - - - -

 if( std::all_of( B.begin() , B.end() ,
		  []( c_FNumber bi ) { return( bi == 0 ); } ) )
  B.clear();

 if( std::all_of( U.begin() , U.end() ,
		  []( c_FNumber ui ) { return( ui == Inf< FNumber >() ); } ) )
  U.clear();
 
 // allocate flow variables - - - - - - - - - - - - - - - - - - - - - - - - -

 generate_abstract_variables();

 // the arcs whose cost is infinite have to be closed
 // in addition the cost has to be set to 0 - - - - - - - - - - - - - - - - -

 for( Index j = 0 ; j < C.size() ; ++j )
  if( C[ j ] >= Inf< double >() ) {
   close_arc( j , eNoMod , eNoMod );
   C[ j ] = 0;
   }

 // issue Modification- - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // note: this is a NBModification, the "nuclear option"

 if( anyone_there() )
  add_Modification( std::make_shared< NBModification >( this ) );

 }  // end( SingleFlowDCRBlock::load( istream ) )

/*--------------------------------------------------------------------------*/

void SingleFlowDCRBlock::load_dcr( std::istream & inFile , Index NNodes, Index NArcs )
{
 // read file .dcr (parameters for DCR part of the problem)
 /*
 int l = strlen( FN );
 char *Name = new char[ l + 5 ];  // temporary string containing the constant
 strcpy( Name , FN );             // part of the pathname + space for `.dcr'
 strcpy( Name + l , ".dcr" );

 std::ifstream inFile( Name );
 if( ! inFile.is_open() )
	throw( std::invalid_argument( "invalid cannot open .dcr file" ) );
 */	

 NodeDelays.resize( NNodes );
 LinkDelays.resize( NArcs );

 Index i,j;
 //read node delays
 for(i = 0; i < NNodes; i++){
  if( ! ( inFile >> NodeDelays[ i ] ) )
     throw( std::invalid_argument( "error reading node delays" ) );
 }

 //read link delays
 for(j = 0; j < NArcs; j++){
  if( ! ( inFile >> LinkDelays[ j ] ) )
     throw( std::invalid_argument( "error reading link delays" ) );
 }
 		
 //read flow burst and deadlines
 if( ! ( inFile >> FlowBursts ) )
     throw( std::invalid_argument( "error reading flow bursts" ) );

 if( ! ( inFile >> FlowDeadlines ) )
     throw( std::invalid_argument( "error reading flow deadlines" ) ); 
			
//read MTU;
 if( ! ( inFile >> MTU ) )
     throw( std::invalid_argument( "error reading MTU" ) ); 

//read rho;
 if( ! ( inFile >> rho ) )
     throw( std::invalid_argument( "error reading rho" ) ); 

 //inFile.close();
 //delete [] Name;

} // end SingleFlowDCRBlock::load_dcr( const char *const FN , int NNodes, int NArcs )

/*--------------------------------------------------------------------------*/

void SingleFlowDCRBlock::deserialize( const netCDF::NcGroup & group )
{
 // erase previous instance, if any- - - - - - - - - - - - - - - - - - - - - -

 if( MaxNNodes || get_MaxNArcs() )
  guts_of_destructor();
		   
 // read problem data- - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 netCDF::NcDim nn = group.getDim( "NNodes" );
 if( nn.isNull() )
  throw( std::logic_error( "NNodes dimension is required" ) );
 NNodes = nn.getSize();

 netCDF::NcDim na = group.getDim( "NArcs" );
 if( na.isNull() )
  throw( std::logic_error( "NArcs dimension is required" ) );
 NArcs = na.getSize();

 MaxNNodes = NNodes;
 Index MaxNArcs = NArcs;
 
 netCDF::NcVar sn = group.getVar( "SN" );
 if( sn.isNull() )
  throw( std::logic_error( "Starting Nodes not found" ) );

 SN.resize( MaxNArcs );

 sn.getVar( SN.data() );

 netCDF::NcVar en = group.getVar( "EN" );
 if( en.isNull() )
  throw( std::logic_error( "Ending Nodes not found" ) );

 EN.resize( MaxNArcs );
 en.getVar( EN.data() );

 netCDF::NcVar cst = group.getVar( "C" );
 if( ! cst.isNull() ) {
  C.resize( MaxNArcs );
  cst.getVar( C.data() );
  }
 else
  C.assign( MaxNArcs , 0 );

 netCDF::NcVar cap = group.getVar( "U" );
 if( ! cap.isNull() ) {
  U.resize( MaxNArcs );
  cap.getVar( U.data() );
  if( std::all_of( U.begin() , U.begin() + NArcs ,
		   []( c_FNumber ui ) { return( ui == Inf< FNumber >() ); } ) )
   U.clear();
  }
/*
 netCDF::NcVar dfc = group.getVar( "B" );
 if( ! dfc.isNull() ) {
  B.resize( MaxNNodes );
  std::vector< size_t > countn = { NNodes };
  dfc.getVar( B.data() );
  if( std::all_of( B.begin() , B.begin() + NNodes ,
		   []( c_FNumber bi ) { return( bi == 0 ); } ) )
   B.clear();
  }
*/
 f_cond_lower = dNAN;  // reset conditional bounds

 // allocate flow variables - - - - - - - - - - - - - - - - - - - - - - - - -

 generate_abstract_variables();

 // the arcs whose cost is infinite have to be closed
 // in addition the cost has to be set to 0 - - - - - - - - - - - - - - - - -

 for( Index j = 0 ; j < C.size() ; ++j )
  if( C[ j ] >= Inf< CNumber >() ) {
   close_arc( j , eNoMod , eNoMod );
   C[ j ] = 0;
   }

 // call the method of Block- - - - - - - - - - - - - - - - - - - - - - - - -
 // inside this the NBModification, the "nuclear option",  is issued

 Block::deserialize( group );

 }  // end( SingleFlowDCRBlock::deserialize )

/*--------------------------------------------------------------------------*/

void SingleFlowDCRBlock::generate_abstract_variables( Configuration *stvv )
{
 if( AR & HasVar )  // the variables are there already
  return;           // nothing to do

 if( HasStaticX() ) {
  x.resize( get_NArcs() );
  for( auto & var : x ){
   var.set_type( ColVariable::kBinary );
   //var.set_type( ColVariable::kNonNegative );
  }

  add_static_variable( x );

  r.resize( get_NArcs() );
  for( auto & var : r )
   var.set_type( ColVariable::kNonNegative );

  add_static_variable( r );

  theta.resize( get_NArcs() );
  for( auto & var : theta )
   var.set_type( ColVariable::kNonNegative );

  add_static_variable( theta );
  }

  r_min.set_type( ColVariable::kNonNegative );

  add_static_variable( r_min );

  theta_min.set_type( ColVariable::kNonNegative );

  add_static_variable( theta_min );

 AR |= HasVar;

 }  // end( SingleFlowDCRBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void SingleFlowDCRBlock::generate_dynamic_constraints( Configuration *stcc )
{

 Index FormMsk = 1;  // PCuts formulation
 if( ( ! stcc ) && f_BlockConfig )
  stcc = f_BlockConfig->f_static_variables_Configuration;
 if( auto sci = dynamic_cast< SimpleConfiguration< int > * >( stcc ) )
  FormMsk = sci->f_value;

 if( FormMsk == PCuts ) {
  double tol = 1e-5;  // threshold parameter for P/C separation
  double eps = 1e-4;  // tolerance value to consider a binary variable

  auto extract_parameters = [ & tol , & eps ]( Configuration * c )
    -> bool {
    if( auto tc = dynamic_cast< SimpleConfiguration< double > * >( c ) ) {
      tol = tc->f_value;
      return( true );
    }
    if( auto tc = dynamic_cast<
      SimpleConfiguration< std::pair< double , double > > * >( c ) ) {
      tol = tc->f_value.first;
      eps = tc->f_value.second;
      return( true );
    }
    return( false );
    };

    if( ( ! extract_parameters( stcc ) ) && f_BlockConfig )
    // if the given Configuration is not valid, try the one from the BlockConfig
    extract_parameters( f_BlockConfig->f_dynamic_constraints_Configuration );

  LinearFunction::v_coeff_pair v_var;
  for(Index j = 0; j < get_NArcs(); j++){
    if( x[ j ].get_value() > eps ) {
      if( theta[ j ].get_value() <
          MTU * ( std::pow( x[ j ].get_value() , 2 ) /
            r[j].get_value()) - tol ) {
              std::list< FRowConstraint > cut( 1 ); 
              v_var.push_back( std::make_pair( &theta[ j ] , -1.0 ));
              v_var.push_back( std::make_pair( &x[j] , MTU * 2 * x[ j ].get_value() /  r[ j ].get_value() ));
              v_var.push_back( std::make_pair( &r[j] , - MTU * ( ( std::pow( x[ j ].get_value() , 2 ) / std::pow( r[ j ].get_value() , 2)))));
              LinearFunction* Funct = new LinearFunction( std::move( v_var ));
              cut.front().set_lhs( -Inf< double >() );
              cut.front().set_rhs( ( 0.0 )); 
              cut.front().set_function( Funct );
              add_dynamic_constraints( PC_cuts , cut , eNoBlck );
    }
    }
  }
  add_dynamic_constraint( PC_cuts , "PC_cuts" );  

  LinearFunction::v_coeff_pair v_var_min;
  //if( r_min.get_value() > eps ) {
      if( theta_min.get_value() <
          ( FlowBursts /
            r_min.get_value() ) - tol ) {
              std::list< FRowConstraint > cut_min( 1 ); 
              v_var_min.push_back( std::make_pair( &theta_min , -1.0 ));
              v_var_min.push_back( std::make_pair( &r_min , -( FlowBursts / std::pow( r_min.get_value() , 2 ))));
              LinearFunction* Funct_min = new LinearFunction( std::move( v_var_min ));
              cut_min.front().set_lhs( -Inf< double >() );
              cut_min.front().set_rhs( -2.0 * FlowBursts / r_min.get_value() ); 
              cut_min.front().set_function( Funct_min );
              add_dynamic_constraints( PC_cuts_min , cut_min , eNoBlck );
    }
    //}
    add_dynamic_constraint( PC_cuts_min , "PC_cuts_min" );
  }
 }// end( SingleFlowDCRBlock::generate_dynamic_constraints )

/*--------------------------------------------------------------------------*/

void SingleFlowDCRBlock::generate_abstract_constraints( Configuration *stcc )
{

 if( AR & HasFlw )  // bound constraints there already
  return;           // nothing to do

 if( ! ( AR & HasFlw ) ) {

  // count number of nonzeroes in each constraint, i.e., #FS( i ) + #BS( i )
  Subset count( get_NNodes() );
 
  for( Index i = 0 ; i < get_NArcs() ; ++i ) {
   count[ SN[ i ] - 1 ]++;
   count[ EN[ i ] - 1 ]++;
   }

  // initialize the vectors of coefficients, and reset count[]
  std::vector< LinearFunction::v_coeff_pair > coeffs( get_NNodes() );

  for( Index i = 0 ; i < get_NNodes() ; ++i ) {
   coeffs[ i ].resize( count[ i ] );
   count[ i ] = 0;
   }

  // construct the vector of coefficients, static phase
  Index i = 0;
  for( ; i < get_NArcs() ; ++i ) {
   coeffs[ SN[ i ] - 1 ][ count[ SN[ i ] - 1 ]++ ] =
                                    std::make_pair( &x[ i ] , double( -1 ) );
   coeffs[ EN[ i ] - 1 ][ count[ EN[ i ] - 1 ]++ ] =
                                    std::make_pair( &x[ i ] , double( 1 ) );
   }

  // generate the node-arc incidence matrix - - - - - - - - - - - - - - - - -
  // each constraint is an equality, i.e., LHS = RHS = B[ i ]

  // static part
  if( HasStaticE() ) {
   E.resize( get_NNodes() );

   for( Index i = 0 ; i < get_NNodes() ; ++i ) {
    E[ i ].set_both( B.empty() ? 0 : B[ i ] );
    E[ i ].set_function( new LinearFunction( std::move( coeffs[ i ] ) , 0 ) );
    }

   add_static_constraint( E );
   }
 }

  AR |= HasFlw;

 // generate the DCR constraint- - - - - - - - - - - - - - - - - - - - - -

 LinearFunction::v_coeff_pair v_var;
 for( Index j = 0 ; j < get_NArcs() ; ++j ) {
   v_var.push_back( std::make_pair( &theta[ j ], 1.0 ));
   v_var.push_back( std::make_pair( &x[ j ] , MTU / U[ j ] + LinkDelays[ j ] + NodeDelays[ SN[ j ] - 1 ]));
  }
  v_var.push_back( std::make_pair( &theta_min ,  1.0 ));
  LinearFunction* Funct = new LinearFunction( std::move( v_var ));
  DCR_cnst.set_rhs( FlowDeadlines );
  DCR_cnst.set_lhs( -Inf< double >() ); 
  DCR_cnst.set_function( Funct );
  add_static_constraint( DCR_cnst );

 Index FormMsk = 2;  // PCuts formulation
 if( ( ! stcc ) && f_BlockConfig )
  stcc = f_BlockConfig->f_static_variables_Configuration;
 if( auto sci = dynamic_cast< SimpleConfiguration< int > * >( stcc ) )
  FormMsk = sci->f_value;

 if( FormMsk == SOCP ){

    DQuadFunction::v_coeff_triple v_vars_q;
    DQuadFunction::coeff_triple t1( &theta_min , 0.0 , 0.0 );
    DQuadFunction::coeff_triple t2( &r_min , 0.0 , 0.0 );
    v_vars_q.reserve( 2 );
    v_vars_q.push_back( t1 );
    v_vars_q.push_back( t2 );
    QuadFunction::v_off_diag_term v_nd_vars;
    v_nd_vars.reserve( 1 );
    QuadFunction::off_diag_term term_min( 1 , 0 , 1.0 );
    v_nd_vars.push_back( term_min );
    QuadFunction* quad = new QuadFunction( std::move(v_vars_q), std::move(v_nd_vars) );
    cone_min_cnst.set_rhs( Inf< double >() );
    cone_min_cnst.set_lhs( FlowBursts ); 
    cone_min_cnst.set_function( quad );

    add_static_constraint( cone_min_cnst, "cone_min" );
 
    cone_cnst.resize( get_NArcs() );
    for( Index j = 0 ; j < get_NArcs() ; ++j ) {
      DQuadFunction::v_coeff_triple v_vars_q1;
      DQuadFunction::coeff_triple t1( &theta[j] , 0.0 , 0.0 );
      DQuadFunction::coeff_triple t2( &r[j] , 0.0 , 0.0 );
      DQuadFunction::coeff_triple t3( &x[j] , 0.0, -MTU );
      v_vars_q1.reserve( 3 );
      v_vars_q1.push_back( t1 );
      v_vars_q1.push_back( t2 );
      v_vars_q1.push_back( t3 );
      QuadFunction::v_off_diag_term v_nd_vars1;
      v_nd_vars1.reserve( 1 );
      QuadFunction::off_diag_term term_min( 1 , 0 , 1.0 );
      v_nd_vars1.push_back( term_min );
      QuadFunction* quad1 = new QuadFunction( std::move(v_vars_q1), std::move(v_nd_vars1) );
      cone_cnst[ j ].set_rhs( Inf< double >() );
      cone_cnst[ j ].set_lhs( 0.0 ); 
      cone_cnst[ j ].set_function( quad1 );

    }

    add_static_constraint( cone_cnst, "cone_cnst" );
  }

 int j;
 double U_max = 0.0 ;
 for( Index i = 0 ; i < get_NArcs() ; ++i )
   U_max = std::max( U_max , U[ i ] );

 // generate the indicator constraints- - - - - - - - - - - - - - - - - - - - - -
 
 Indicator_cnst_rmin.resize( NArcs );
 
 for(j = 0; j < NArcs; j++){
  LinearFunction::v_coeff_pair v_var;
  v_var.push_back( std::make_pair( &r_min , 1.0 ));
  v_var.push_back( std::make_pair( &x[ j ] , U_max  ));
  v_var.push_back( std::make_pair( &r[ j ] , -1.0 ));
  LinearFunction* Funct = new LinearFunction( std::move( v_var ));
  Indicator_cnst_rmin[ j ].set_lhs( -Inf< double >() );
  Indicator_cnst_rmin[ j ].set_rhs( U_max ); 
  Indicator_cnst_rmin[ j ].set_function( Funct );
 }

 add_static_constraint( Indicator_cnst_rmin );

 Indicator_cnst_r1.resize( NArcs );
 for(j = 0; j < NArcs; j++){
  LinearFunction::v_coeff_pair v_var;
  v_var.push_back( std::make_pair( &x[ j ] , -U[ j ] ));
  v_var.push_back( std::make_pair( &r[ j ] , 1.0 ));
  LinearFunction* Funct = new LinearFunction( std::move( v_var ));
  Indicator_cnst_r1[ j ].set_lhs( -Inf< double >() );
  Indicator_cnst_r1[ j ].set_rhs( 0.0 ); 
  Indicator_cnst_r1[ j ].set_function( Funct );
 }

 add_static_constraint( Indicator_cnst_r1 );

 Indicator_cnst_r2.resize( NArcs );
 for(j = 0; j < NArcs; j++){
  LinearFunction::v_coeff_pair v_var;
  v_var.push_back( std::make_pair( &x[ j ] , rho ));
  v_var.push_back( std::make_pair( &r[ j ] , -1.0 ));
  LinearFunction* Funct = new LinearFunction( std::move( v_var ));
  Indicator_cnst_r2[ j ].set_lhs( -Inf< double >() );
  Indicator_cnst_r2[ j ].set_rhs( 0.0 ); 
  Indicator_cnst_r2[ j ].set_function( Funct );
 }

 add_static_constraint( Indicator_cnst_r2 );

 // generate the bound constraints- - - - - - - - - - - - - - - - - - - - - -

 if( AR & HasBnd )  // bound constraints there already
  return;           // nothing to do

 AR |= HasBnd;

 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end( SingleFlowDCRBlock::generate_abstract_constraints )

/*--------------------------------------------------------------------------*/

void SingleFlowDCRBlock::generate_objective( Configuration *objc )
{
 if( AR & HasObj )  // the objective is there already
  return;           // cowardly (and silently) return

 // initialize objective function - - - - - - - - - - - - - - - - - - - - - -

 LinearFunction::v_coeff_pair p;

 // construct a "dense" LinearFunction- - - - - - - - - - - - - - - - - - - -

 Index i = 0;
 auto Cit = C.begin();

 // static part
 //if( HasStaticX() )
 for( ; i < get_NArcs() ; ++i ) {
   p[ i ].first = &r[ i ];
   auto ci = *(Cit++);
   p[ i ].second = 1.0; //std::isnan( ci ) ? 0 : ci;
   }

 i = 0;
 for( ; i < get_NArcs() ; ++i ) {
   p.push_back( std::make_pair( &r[ i ], 1.0 ));
   }

 // ensure no Modification is issued: this may happen in case a SingleFlowDCRBlock
 // is re-loaded, so that set_objective( c ) had already been called

 c.set_function( new LinearFunction( std::move( p ) , 0 ) , eNoMod );

 set_objective( & c , eNoMod );

 AR |= HasObj;

 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end( SingleFlowDCRBlock::generate_objective )

/*--------------------------------------------------------------------------*/
/*--------------------- Methods for checking the Block ---------------------*/
/*--------------------------------------------------------------------------*/

bool SingleFlowDCRBlock::flow_feasible( c_FNumber feps , bool useabstract )
{
 if( useabstract && ( AR & HasFlw ) ) {
  // do it using the abstract representation, if possible - - - - - - - - - -

  if( ! RowConstraint::is_feasible( E , feps ) )   // static part
   return( false );
 }

 return( true );

 }  // end( SingleFlowDCRBlock::flow_feasible )

/*--------------------------------------------------------------------------*/

bool SingleFlowDCRBlock::bound_feasible( c_FNumber feps , bool useabstract )
{
 if( useabstract && ( AR & ( HasFlw | HasVar ) ) ) {
  // do it using the abstract representation, if possible - - - - - - - - - -

  // static part
  if( HasStaticX() ) {
   if( UB.empty() ) {
    if( ! ColVariable::is_feasible( x , feps ) )
     return( false );
    }
   else
    if( ! RowConstraint::is_feasible( UB , feps ) )
     return( false );
   }
  }
 else {
  // do it using the physical representation- - - - - - - - - - - - - - - - -
  Index i = 0;

  // static part
  for( ; i < get_NArcs() ; ++i ) {
   c_FNumber Ui = get_U( i );
   c_FNumber xi = x[ i ].get_value();
   if( Ui >= Inf< FNumber >() ) {
    if( xi < - feps )
     return( false );
    }
   else {
    c_FNumber slck = Ui == 0 ? std::abs( xi ) :
                               std::max( - xi , xi - Ui ) / std::abs( Ui );
    if( slck > feps )
     return( false );
    }
   }
  }

 return( true );

 }  // end( SingleFlowDCRBlock::bound_feasible )

/*--------------------------------------------------------------------------*/

bool SingleFlowDCRBlock::is_feasible( bool useabstract , Configuration *fsbc )
{
 // Retrieve the tolerance and the type of violation.
 double tol = 1e-3;
 bool rel_viol = true;

 // Try to extract, from "c", the parameters that determine feasibility.
 // If it succeeds, it sets the values of the parameters and returns
 // true. Otherwise, it returns false.
 auto extract_parameters = [ & tol , & rel_viol ]( Configuration * c )
  -> bool {
  if( auto tc = dynamic_cast< SimpleConfiguration< double > * >( c ) ) {
   tol = tc->f_value;
   return( true );
  }
  if( auto tc = dynamic_cast< SimpleConfiguration< std::pair< double , int > > * >( c ) ) {
   tol = tc->f_value.first;
   rel_viol = tc->f_value.second;
   return( true );
  }
  return( false );
 };

 if( ( ! extract_parameters( fsbc ) ) && f_BlockConfig )
  // if the given Configuration is not valid, try the one from the BlockConfig
  extract_parameters( f_BlockConfig->f_is_feasible_Configuration );

  return(
  // Constraints: notice that the ZOConstraints are not checked, since the
  // corresponding check is made on the ColVariable
  RowConstraint::is_feasible( E , tol , rel_viol )
  && RowConstraint::is_feasible( DCR_cnst , tol , rel_viol )
  && RowConstraint::is_feasible( Indicator_cnst_rmin , tol , rel_viol )
  && RowConstraint::is_feasible( Indicator_cnst_r1 , tol , rel_viol )
  && RowConstraint::is_feasible( Indicator_cnst_r2 , tol , rel_viol )
  && RowConstraint::is_feasible( cone_min_cnst , tol , rel_viol )
  && RowConstraint::is_feasible( cone_cnst , tol , rel_viol )
  );

 }  // end( SingleFlowDCRBlock::is_feasible )

/*--------------------------------------------------------------------------*/

bool SingleFlowDCRBlock::is_feasible_flow( bool useabstract , Configuration *fsbc )
{
 // Retrieve the tolerance and the type of violation.
 double tol = 1e-3;
 bool rel_viol = true;

 // Try to extract, from "c", the parameters that determine feasibility.
 // If it succeeds, it sets the values of the parameters and returns
 // true. Otherwise, it returns false.
 auto extract_parameters = [ & tol , & rel_viol ]( Configuration * c )
  -> bool {
  if( auto tc = dynamic_cast< SimpleConfiguration< double > * >( c ) ) {
   tol = tc->f_value;
   return( true );
  }
  if( auto tc = dynamic_cast< SimpleConfiguration< std::pair< double , int > > * >( c ) ) {
   tol = tc->f_value.first;
   rel_viol = tc->f_value.second;
   return( true );
  }
  return( false );
 };

 if( ( ! extract_parameters( fsbc ) ) && f_BlockConfig )
  // if the given Configuration is not valid, try the one from the BlockConfig
  extract_parameters( f_BlockConfig->f_is_feasible_Configuration );

  return(
  // Constraints: notice that the ZOConstraints are not checked, since the
  // corresponding check is made on the ColVariable
  RowConstraint::is_feasible( E , tol , rel_viol )
  && RowConstraint::is_feasible( DCR_cnst , tol , rel_viol )
  && RowConstraint::is_feasible( Indicator_cnst_rmin , tol , rel_viol )
  && RowConstraint::is_feasible( Indicator_cnst_r1 , tol , rel_viol )
  && RowConstraint::is_feasible( Indicator_cnst_r2 , tol , rel_viol )
  );

 }  // end( SingleFlowDCRBlock::is_feasible )

/*--------------------------------------------------------------------------*/
/*------------------------- Methods for R3 Blocks --------------------------*/
/*--------------------------------------------------------------------------*/

Block * SingleFlowDCRBlock::get_R3_Block( Configuration *r3bc , Block * base  ,
				Block * father )
{
 if( r3bc != nullptr )
  throw( std::invalid_argument( "non-nullptr R3B Configuration" ) );

 SingleFlowDCRBlock *DCRB;
 if( base ) {
  DCRB = dynamic_cast< SingleFlowDCRBlock * >( base );
  if( ! DCRB )
   throw( std::invalid_argument( "base is not a SingleFlowDCRBlock" ) );
  }
 else
  DCRB = new SingleFlowDCRBlock( father );

 DCRB->load( get_NNodes() , get_NArcs() , EN , SN , U , C ,
       NodeDelays , LinkDelays , FlowBursts , FlowDeadlines );
 
 return( DCRB );

 }  // end( SingleFlowDCRBlock::get_R3_Block )

/*--------------------------------------------------------------------------*/

void SingleFlowDCRBlock::map_back_solution( Block *R3B , Configuration *r3bc ,
				               Configuration *solc )
{
 // process Configuration - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 auto DCRB = dynamic_cast< SingleFlowDCRBlock * >( R3B );
 if( ! DCRB )
  throw( std::invalid_argument( "R3B is not a SingleFlowDCRBlock" ) );
 if( r3bc != nullptr )
  throw( std::invalid_argument( "non-nullptr R3B Configuration" ) );

 int wsol = 0;
 auto tsolc = dynamic_cast< SimpleConfiguration< int > * >( solc );

 if( ( ! tsolc ) && f_BlockConfig && f_BlockConfig->f_solution_Configuration )
   tsolc = dynamic_cast< SimpleConfiguration< int > * >(
                                    f_BlockConfig->f_solution_Configuration );
 if( tsolc )
  wsol = tsolc->f_value;

 // if required, map back primal solution - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ( wsol != 2 ) && ( AR & HasVar ) ) {  // ... if any
  if( DCRB->get_NArcs() != get_NArcs() )
   throw( std::invalid_argument( "incompatible static flow size" ) );

  // static part
  if( HasStaticX() )
   for( auto xi = x.begin() , r3bxi = DCRB->x.begin() ; xi != x.end() ;
	++xi , ++r3bxi )
    if( ! xi->is_fixed() )
     xi->set_value( r3bxi->get_value() );
  }
 }  // end( SingleFlowDCRBlock::map_back_solution )

/*--------------------------------------------------------------------------*/

void SingleFlowDCRBlock::map_forward_solution( Block *R3B , Configuration *r3bc ,
				                  Configuration *solc )
{
 // process Configuration - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 auto DCRB = dynamic_cast< SingleFlowDCRBlock * >( R3B );
 if( ! DCRB )
  throw( std::invalid_argument( "R3B is not a SingleFlowDCRBlock" ) );
 if( r3bc != nullptr )
  throw( std::invalid_argument( "non-nullptr R3B Configuration" ) );

 int wsol = 0;
 auto tsolc = dynamic_cast< SimpleConfiguration< int > * >( solc );

 if( ( ! tsolc ) && f_BlockConfig && f_BlockConfig->f_solution_Configuration )
  tsolc = dynamic_cast< SimpleConfiguration< int > * >(
                                    f_BlockConfig->f_solution_Configuration );
 if( tsolc )
  wsol = tsolc->f_value;

 // if required, map forward primal solution- - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ( wsol != 1 ) && ( AR & HasFlw ) ) {  // ... if any
  if( DCRB->get_NArcs() != get_NArcs() )
   throw( std::invalid_argument( "incompatible static flow size" ) );

  // static part
  if( DCRB->HasStaticX() )
   for( auto xi = x.begin() , r3bxi = DCRB->x.begin() ; xi != x.end() ;
	++xi , ++r3bxi )
    if( ! r3bxi->is_fixed() )
     r3bxi->set_value( xi->get_value() );
  }
 }  // end( SingleFlowDCRBlock::map_forward_solution )

/*--------------------------------------------------------------------------*/

bool SingleFlowDCRBlock::map_forward_Modification( Block * R3B , c_p_Mod mod ,
					 Configuration * r3bc ,
					 ModParam issuePMod ,
					 ModParam issueAMod )
{
 if( mod->concerns_Block() )  // an abstract Modification
  return( false );            // none of my business
 
 auto DCRB = dynamic_cast< SingleFlowDCRBlock * >( R3B );
 if( ! DCRB )
  throw( std::invalid_argument( "R3B is not a SingleFlowDCRBlock" ) );
 if( r3bc != nullptr )
  throw( std::invalid_argument( "non-nullptr R3B Configuration" ) );

 auto iPM = issuePMod;
 auto iPA = un_ModBlock( issueAMod );

 /* Use a Lambda to define a "guts" of the method that can be called
    recursively without having to pass "local globals". Note the trick of
    defining the std::function object and "passing" it to the lambda,
    which allows recursive calls. Note the need to explicitly capture
    "this" to use fields/methods of the class. */

 std::function< bool( c_p_Mod ) > guts_of_mfM;
 guts_of_mfM = [ this , & guts_of_mfM , & DCRB , & iPM , & iPA ]( c_p_Mod mod
								  ) {
  // process Modification- - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  /* This requires to patiently sift through the possible Modification types
     to find what this Modification exactly is, and call the appropriate
     method of either DCRB, for a "physical Modification", or of the "abstract
     representation" of DCRB for an "abstract Modification". */

  //!! std::cout << *mod << std::endl;
  
  // GroupModification - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if( auto tmod = dynamic_cast< const GroupModification * >( mod ) ) {
   // open or nest the two channels
   iPM = make_par( par2mod( iPM ) , DCRB->open_channel( par2chnl( iPM ) ) );
   iPA = make_par( par2mod( iPA ) , DCRB->open_channel( par2chnl( iPA ) ) );

   bool ok = true;
   for( const auto & submod : tmod->sub_Modifications() )
    if( ! guts_of_mfM( submod.get() ) )
     ok = false;

   // close or un-nest the channels
   DCRB->close_channel( par2chnl( iPM ) );
   DCRB->close_channel( par2chnl( iPA ) );

   return( ok );
   }

  // SingleFlowDCRBlockRngdMod - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  /* Note: in the following we can assume that C, B and U are nonempty. This
     is because they can be empty only if they are so when the object is
     loaded. But if a Modification has been issued they are no longer empty
     (a Modification changing nothing from the "empty" state is not issued). */

  if( auto tmod = dynamic_cast< const SingleFlowDCRBlockRngdMod * >( mod ) ) {
   switch( tmod->type() ) {
    case( SingleFlowDCRBlockMod::eChgCost ):
     #ifndef NDEBUG
      if( ( tmod->rng().second > get_NArcs() ) ||
	  ( tmod->rng().second > DCRB->get_NArcs() ) )
       throw( std::logic_error(
		     "map_forward_Modification:: incompatible SingleFlowDCRBlock" ) );
     #endif
     if( tmod->rng().second == tmod->rng().first + 1 )
      DCRB->chg_cost( C[ tmod->rng().first ] , tmod->rng().first ,
		      iPM , iPA );
     else
      DCRB->chg_costs( C.begin() + tmod->rng().first , tmod->rng() ,
		       iPM , iPA );
     break;
    case( SingleFlowDCRBlockMod::eChgCaps ):
     #ifndef NDEBUG
      if( ( tmod->rng().second > get_NArcs() ) ||
	  ( tmod->rng().second > DCRB->get_NArcs() ) )
       throw( std::logic_error(
		     "map_forward_Modification:: incompatible SingleFlowDCRBlock" ) );
     #endif
     if( tmod->rng().second == tmod->rng().first + 1 )
      DCRB->chg_ucap( U.empty() ? Inf< FNumber >() : U[ tmod->rng().first ] ,
		      tmod->rng().first , iPM , iPA );
     else
      if( U.empty() ) {
       Vec_FNumber NCap( tmod->rng().second - tmod->rng().first );
       auto NCit = NCap.begin();
       for( Index i = tmod->rng().first ; i < tmod->rng().second ; ++i )
	*(NCit++) = U[ i ];

       DCRB->chg_ucaps( NCap.begin() , tmod->rng() , iPM , iPA );
       }
      else
       DCRB->chg_ucaps( U.begin() + tmod->rng().first , tmod->rng() ,
			iPM , iPA );
     break;

    case( SingleFlowDCRBlockMod::eChgDfct ):
     #ifndef NDEBUG
      if( ( tmod->rng().second > get_NNodes() ) ||
	  ( tmod->rng().second > DCRB->get_NNodes() ) )
       throw( std::logic_error(
		     "map_forward_Modification:: incompatible SingleFlowDCRBlock" ) );
     #endif
     if( tmod->rng().second == tmod->rng().first + 1 )
      DCRB->chg_dfct( B.empty() ? 0 : B[ tmod->rng().first ] ,
		      tmod->rng().first , iPM , iPA );
     else
      if( B.empty() ) {
       Vec_FNumber NDfct( tmod->rng().second - tmod->rng().first );
       auto NDit = NDfct.begin();
       for( Index i = tmod->rng().first ; i < tmod->rng().second ; ++i )
	*(NDit++) = B[ i ];

       DCRB->chg_ucaps( NDfct.begin() , tmod->rng() , iPM , iPA );
       }
      else
       DCRB->chg_dfcts( B.begin() + tmod->rng().first , tmod->rng() ,
			iPM , iPA );
     break;

    case( SingleFlowDCRBlockMod::eOpenArc ):
     #ifndef NDEBUG
      if( ( tmod->rng().second > get_NArcs() ) ||
	  ( tmod->rng().second > DCRB->get_NArcs() ) )
       throw( std::logic_error(
		     "map_forward_Modification:: incompatible SingleFlowDCRBlock" ) );
     #endif
     if( tmod->rng().second == tmod->rng().first + 1 )
      DCRB->open_arc( tmod->rng().first , iPM , iPA );
     else
      DCRB->open_arcs( tmod->rng() , iPM , iPA );
     break;
    case( SingleFlowDCRBlockMod::eCloseArc ):
     #ifndef NDEBUG
      if( ( tmod->rng().second > get_NArcs() ) ||
	  ( tmod->rng().second > DCRB->get_NArcs() ) )
       throw( std::logic_error(
		     "map_forward_Modification:: incompatible SingleFlowDCRBlock" ) );
     #endif
     if( tmod->rng().second == tmod->rng().first + 1 )
      DCRB->close_arc( tmod->rng().first , iPM , iPA );
     else
      DCRB->close_arcs( tmod->rng() , iPM , iPA );
     break;
    default:
     throw( std::invalid_argument( "unknown SingleFlowDCRBlockRngdMod type" ) );
    }
   return( true );
   }

  // SingleFlowDCRBlockSbstMod - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  /* Note that tmod->nms() need be copied, since the chg_*() methods
   * *in principle* "consume" the names vector. This is actually not true
   * if DCRB will *not* issue a physical modification, which one may
   * actually know beforehand, but it has to be done anyway because the
   * SingleFlowDCRBlockSbstMod only provides read-only access to the vector.
   * However, tmod->nms() is guaranteed to be ordered. */

  if( auto tmod = dynamic_cast< const SingleFlowDCRBlockSbstMod * >( mod ) ) {
   switch( tmod->type() ) {
    case( SingleFlowDCRBlockMod::eChgCost ): {
     #ifndef NDEBUG
      if( ( tmod->nms().back() >= get_NArcs() ) ||
	  ( tmod->nms().back() >= DCRB->get_NArcs() ) )
       throw( std::logic_error(
		     "map_forward_Modification:: incompatible SingleFlowDCRBlock" ) );
     #endif
     Vec_CNumber NCost( tmod->nms().size() );
     for( Index i = 0 ; i < NCost.size() ; i++ )
      NCost[ i ] = C[ tmod->nms()[ i ] ];

     DCRB->chg_costs( NCost.begin() , Subset( tmod->nms() ) , true ,
		      iPM , iPA );
     break;
     }
    case( SingleFlowDCRBlockMod::eChgCaps ): {
     #ifndef NDEBUG
      if( ( tmod->nms().back() >= get_NArcs() ) ||
	  ( tmod->nms().back() >= DCRB->get_NArcs() ) )
       throw( std::logic_error(
		     "map_forward_Modification:: incompatible SingleFlowDCRBlock" ) );
     #endif
     Vec_FNumber NCap( tmod->nms().size() , Inf< FNumber >() );
     if( ! U.empty() )
      for( Index i = 0 ; i < NCap.size() ; i++ )
       NCap[ i ] = U[ tmod->nms()[ i ] ];

     DCRB->chg_ucaps( NCap.begin() , Subset( tmod->nms() ) , true ,
		      iPM , iPA );
     break;
     }

    case( SingleFlowDCRBlockMod::eChgDfct ): {
     #ifndef NDEBUG
      if( ( tmod->nms().back() >= get_NNodes() ) ||
	  ( tmod->nms().back() >= DCRB->get_NNodes() ) )
       throw( std::logic_error(
		     "map_forward_Modification:: incompatible SingleFlowDCRBlock" ) );
     #endif
     Vec_FNumber NDfct( tmod->nms().size() , 0 );
     if( ! B.empty() )
      for( Index i = 0 ; i < NDfct.size() ; i++ )
       NDfct[ i ] = B[ tmod->nms()[ i ] ];

     DCRB->chg_dfcts( NDfct.begin() , Subset( tmod->nms() ) , true ,
		      iPM , iPA );
     break;
     }

    case( SingleFlowDCRBlockMod::eOpenArc ):
     #ifndef NDEBUG
      if( ( tmod->nms().back() >= get_NArcs() ) ||
	  ( tmod->nms().back() >= DCRB->get_NArcs() ) )
       throw( std::logic_error(
		     "map_forward_Modification:: incompatible SingleFlowDCRBlock" ) );
     #endif
    DCRB->open_arcs( Subset( tmod->nms() ) , true , iPM , iPA );
     break;
    case( SingleFlowDCRBlockMod::eCloseArc ):
     #ifndef NDEBUG
      if( ( tmod->nms().back() >= get_NArcs() ) ||
	  ( tmod->nms().back() >= DCRB->get_NArcs() ) )
       throw( std::logic_error(
		     "map_forward_Modification:: incompatible SingleFlowDCRBlock" ) );
     #endif
    DCRB->close_arcs( Subset( tmod->nms() ) , true , iPM , iPA );
     break;
    default:
     throw( std::invalid_argument( "unknown SingleFlowDCRBlockSbstMod type" ) );
    }
   return( true );
   }

  // NBModification- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // this is the "nuclear option": the SingleFlowDCRBlock has been re-loaded
  // one should check that the Block is this SingleFlowDCRBlock, but it cannot
  // be otherwise, can it?

  if( auto tmod = dynamic_cast< const NBModification * >( mod ) ) {
   DCRB->load( get_NNodes() , get_NArcs() , EN , SN , U , C ,
         NodeDelays , LinkDelays , FlowBursts , FlowDeadlines );
   return( true );
   }

  return( false );

  };  // end( guts_of_mfM )- - - - - - - - - - - - - - - - - - - - - - - - - -
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // finally, call the "guts of"- - - - - - - - - - - - - - - - - - - - - - - -
 return( guts_of_mfM( mod ) );

 }  // end( SingleFlowDCRBlock::map_forward_Modification )

/*--------------------------------------------------------------------------*/

bool SingleFlowDCRBlock::map_back_Modification( Block *R3B , c_p_Mod mod ,
				      Configuration *r3bc ,
				      ModParam issuePMod ,
				      ModParam issueAMod )
{
 /* Fantastically dirty trick: because the two objects are copies, mapping
  * back a Modification to this from R3B is the same as mapping forward a
  * Modification from R3B to this. */

 auto DCRB = dynamic_cast< SingleFlowDCRBlock * >( R3B );
 if( ! DCRB )
  throw( std::invalid_argument( "R3B is not a SingleFlowDCRBlock" ) );

 return( DCRB->map_forward_Modification( this , mod , r3bc , issuePMod ,
					 issueAMod ) );

 }  // end( SingleFlowDCRBlock::map_back_Modification )

/*--------------------------------------------------------------------------*/
/*----------------------- Methods for handling Solution --------------------*/
/*--------------------------------------------------------------------------*/

Solution * SingleFlowDCRBlock::get_Solution( Configuration * solc , bool emptys )
{
 int wsol = 0;
 if( ( ! solc ) && f_BlockConfig )
  solc = f_BlockConfig->f_solution_Configuration;

 if( auto tsolc = dynamic_cast< SimpleConfiguration< int > * >( solc ) )
  wsol = tsolc->f_value;

 auto *sol = new DCRSolution();

 if( wsol != 2 )
  sol->v_x.resize( get_NArcs() );

 if( wsol != 1 )
  sol->v_r.resize( get_NArcs() );

 if( ! emptys )
  sol->read( this );

 return( sol );

 }  // end( SingleFlowDCRBlock::get_Solution )

/*--------------------------------------------------------------------------*/

void SingleFlowDCRBlock::get_x( Vec_FNumber_it FSol , Range rng ) const
{
 for( ; rng.first < std::min( rng.second , get_NArcs() ) ; )
  *(FSol++) = x[ rng.first++ ].get_value();

 }  // end( SingleFlowDCRBlock::get_x( range ) )

/*--------------------------------------------------------------------------*/

void SingleFlowDCRBlock::get_x( Vec_FNumber_it FSol , c_Subset & nms ) const
{
 if( ! ( AR & HasVar ) )
  throw( std::logic_error( "flow Variable not available" ) );

 auto nmsi = nms.begin();

  while( nmsi != nms.end() ) 
   *(FSol++) = x[ *(nmsi++) ].get_value();

 }  // end( SingleFlowDCRBlock::get_x( subset ) )

/*--------------------------------------------------------------------------*/

void SingleFlowDCRBlock::get_r( Vec_FNumber_it FSol , Range rng ) const
{
 for( ; rng.first < std::min( rng.second , get_NArcs() ) ; )
  *(FSol++) = r[ rng.first++ ].get_value();

 }  // end( SingleFlowDCRBlock::get_r( range ) )

/*--------------------------------------------------------------------------*/

void SingleFlowDCRBlock::get_r( Vec_FNumber_it FSol , c_Subset & nms ) const
{
 if( ! ( AR & HasVar ) )
  throw( std::logic_error( "flow Variable not available" ) );

 auto nmsi = nms.begin();

  while( nmsi != nms.end() ) 
   *(FSol++) = r[ *(nmsi++) ].get_value();

 }  // end( SingleFlowDCRBlock::get_r( subset ) )

/*--------------------------------------------------------------------------*/

void SingleFlowDCRBlock::set_x( c_Vec_FNumber_it fstrt , Range rng )
{
 if( ! ( AR & HasVar ) )  // nowhere to put the value in
  return;                 // cowardly (and silently) return

 if( rng.second > get_NArcs() )
  rng.second = get_NArcs();

 Index i = rng.first;

 if( HasStaticX() )
  for( auto xi = x.begin() + i ;
       i < std::min( rng.second , get_NArcs() ) ; ++i )
   (xi++)->set_value( *(fstrt++) );

 }  // end( SingleFlowDCRBlock::set_x( range ) )

/*--------------------------------------------------------------------------*/

void SingleFlowDCRBlock::set_x( c_Vec_FNumber_it fstrt , c_Subset sbst )
{
 if( ! ( AR & HasVar ) )  // nowhere to put the value in
  return;                 // cowardly (and silently) return

 auto it = sbst.begin();

 if( HasStaticX() )
  while( ( it != sbst.end() ) && ( *it < get_NArcs() ) )
   x[ *(it++) ].set_value( *(fstrt++) );

 return;

 }  // end( SingleFlowDCRBlock::set_x( subset ) )

/*--------------------------------------------------------------------------*/

void SingleFlowDCRBlock::set_r( c_Vec_FNumber_it fstrt , Range rng )
{
 if( ! ( AR & HasVar ) )  // nowhere to put the value in
  return;                 // cowardly (and silently) return

 if( rng.second > get_NArcs() )
  rng.second = get_NArcs();

 Index i = rng.first;

 if( HasStaticX() )
  for( auto xi = r.begin() + i ;
       i < std::min( rng.second , get_NArcs() ) ; ++i )
   (xi++)->set_value( *(fstrt++) );

 return;

 }  // end( SingleFlowDCRBlock::set_r( range ) )

/*--------------------------------------------------------------------------*/

void SingleFlowDCRBlock::set_r( c_Vec_FNumber_it fstrt , c_Subset sbst )
{
 if( ! ( AR & HasVar ) )  // nowhere to put the value in
  return;                 // cowardly (and silently) return

 auto it = sbst.begin();

 if( HasStaticX() )
  while( ( it != sbst.end() ) && ( *it < get_NArcs() ) )
   r[ *(it++) ].set_value( *(fstrt++) );

 return;

 }  // end( SingleFlowDCRBlock::set_x( subset ) )

/*--------------------------------------------------------------------------*/
/*-------------------- Methods for handling Modification -------------------*/
/*--------------------------------------------------------------------------*/

void SingleFlowDCRBlock::add_Modification( sp_Mod mod , ChnlName chnl )
{
 //!! std::cout << *mod << std::endl;

 if( mod->concerns_Block() ) {
  mod->concerns_Block( false );
  guts_of_add_Modification( mod.get() , chnl );
  }

 Block::add_Modification( mod , chnl );
 }

/*--------------------------------------------------------------------------*/
/*------------ METHODS FOR LOADING, PRINTING & SAVING THE SingleFlowDCRBlock*/
/*--------------------------------------------------------------------------*/

void SingleFlowDCRBlock::print( std::ostream  & output , char vlvl ) const
{
 if( vlvl != 'C' ) {  // non-complete version
  // only basic information 
  output << "SingleFlowDCRBlock with: " << NNodes << " nodes and " << SN.size()
	 << " arcs" << std::endl;

  if( ! vlvl ) {     // print the graph

   if( ! B.empty() )
    for( Index i = 0 ; i < get_NNodes() ; ++i )
     if( B[ i ] != 0 )
      output << "B[ " << i + 1 << " ] = " << B[ i ] << std::endl;

   for( Index i = 0 ; i < get_NArcs() ; ++i ) {
    output << "( " << SN[ i ] << " , " << EN[ i ] << " ): C = ";
    if( is_deleted( i ) )
     output << "0";
    else
     output << C[ i ];

    if( ! U.empty() ) {
     output << ", U = ";
     print_UB( output , U[ i ] );
     }

    output << std::endl;
    }
   }
  }
 else  {
  // print header in DIMACS standard format
  output << "p min " << get_NNodes() << " " << get_NArcs() << std::endl;
  output << std::setprecision( 16 );

  // print node descriptors in DIMACS standard format
  if( ! B.empty() )
   for( Index i = 0 ; i < get_NNodes() ; ++i )
    if( B[ i ] != 0 )
     output << "n\t" << i + 1 << "\t" << - B[ i ] << std::endl;
  
  // print arc descriptors in DIMACS standard format
  for( Index i = 0 ; i < get_NArcs() ; ++i ) {
   output << "a\t" << SN[ i ] << "\t" << EN[ i ] << "\t0\t";

   if( ( is_deleted( i ) ) || is_closed( i ) )
    output << "0";
   else {
    if( U.empty() )
     output << "+Inf";
    else
     print_UB( output , U[ i ] );
    }

   output << "\t";
   if( is_deleted( i ) )
    output << "1000000";
   else
    output << C[ i ];
   output << std::endl;
   }
  }
 }  // end( SingleFlowDCRBlock::print )

/*--------------------------------------------------------------------------*/

void SingleFlowDCRBlock::serialize( netCDF::NcGroup & group ) const
{
 // call the method of Block- - - - - - - - - - - - - - - - - - - - - - - - -

 Block::serialize( group );

 // now the SingleFlowDCRBlock data - - - - - - - - - - - - - - - - - - - - - - - - - -

 netCDF::NcDim nn = group.addDim( "NNodes" , get_NNodes() );
 netCDF::NcDim na = group.addDim( "NArcs" , get_NArcs() );

 ( group.addVar( "SN" , netCDF::NcUint64() , na ) ).putVar( SN.data() );

 ( group.addVar( "EN" , netCDF::NcUint64() , na ) ).putVar( EN.data() );

 ( group.addVar( "C" , netCDF::NcDouble() , na ) ).putVar( C.data() );

 if( ! U.empty() )
  ( group.addVar( "U" , netCDF::NcDouble() , na ) ).putVar( U.data() );
/*
 if( ! B.empty() )
  ( group.addVar( "B" , netCDF::NcDouble() , nn ) ).putVar( B.data() );
*/
 }  // end( SingleFlowDCRBlock::serialize )

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/

void SingleFlowDCRBlock::chg_costs( c_Vec_CNumber_it NCost , Range rng ,
			  ModParam issueMod , ModParam issueAMod )
{
 rng.second = std::min( rng.second , get_NArcs() );
 if( rng.second <= rng.first )  // nothing to change
  return;                       // cowardly (and silently) return

 // check to see how many of the initial arcs are either deleted or not
 // really changing the costs
 while( ( std::isnan( C[ rng.first ] ) || ( *NCost == C[ rng.first ] ) )
	&& ( rng.first < rng.second ) ) {
  ++rng.first;
  ++NCost;
  }

 if( rng.second <= rng.first )  // nothing left to change
  return;                       // cowardly (and silently) return

 // check to see how many of the final arcs are either deleted or not
 // really changing the costs
 auto NCEit = NCost + ( rng.second - rng.first );
 while( ( ( std::isnan( C[ rng.second - 1 ] ) ) ||
	  ( *(--NCEit) == C[ rng.second - 1 ] ) )
	&& ( rng.first < rng.second ) )
  --rng.second;

 if( rng.second <= rng.first )  // nothing left to change
  return;                       // cowardly (and silently) return

 if( not_dry_run( issueAMod ) && ( AR & HasObj ) ) {
  // change abstract and physical representation together - - - - - - - - - -
  // in the meantime, if so instructed also issue abstract Modification
  Vec_CNumber NC( rng.second - rng.first );
  auto NCit = NC.begin();

  for( Index i = rng.first ; i < rng.second ; ++i , ++NCost )
   if( std::isnan( C[ i ] ) )  // arc is deleted
    *(NCit++) = 0;             // give it an "harmless" coefficient
   else                        // arc is there    
    *(NCit++) = C[ i ] = *NCost;

  get_lfo()->modify_coefficients( std::move( NC ) , rng ,
				  un_ModBlock( issueAMod ) );
  }
 else
  // only change the physical representation- - - - - - - - - - - - - - - - -
  if( not_dry_run( issueMod ) )
   for( Index i = rng.first ; i < rng.second ; ++i , ++NCost )
    if( ! std::isnan( C[ i ] ) )
     C[ i ] = *NCost;

 f_cond_lower = dNAN;  // reset conditional bounds

 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared< SingleFlowDCRBlockRngdMod >( this ,
					      SingleFlowDCRBlockMod::eChgCost , rng ) ,
			   Observer::par2chnl( issueMod ) );
 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end( SingleFlowDCRBlock::chg_costs( range ) )

/*--------------------------------------------------------------------------*/

void SingleFlowDCRBlock::chg_costs( c_Vec_CNumber_it NCost , Subset && nms ,
			  bool ordered  ,
			  ModParam issueMod , ModParam issueAMod )
{
 if( nms.empty() )  // nothing to change
  return;           // cowardly (and silently) return

 // eliminate from NCost and nms the entries corresponding to either
 // deleted arcs or arcs whose cost actually does not change; meanwhile,
 // if nms is not ordered, order it
 Vec_CNumber NC;
 if( ordered ) {
  NC.resize( nms.size() );
  auto NCit = NC.begin();
  auto nmsit = nms.begin();
  for( auto i : nms ) {
   auto nci = *(NCost++);
   if( ( ! std::isnan( C[ i ] ) ) && ( nci != C[ i ] ) ) {
    *(nmsit++) = i;
    *(NCit++) = nci;
    }
   }
  nms.resize( std::distance( nms.begin() , nmsit ) );
  NC.resize( nms.size() );
  }
 else {
  using TP = std::pair< Index , CNumber >;
  std::vector< TP > pairs;
  pairs.reserve( nms.size() );
  for( auto i : nms ) {
   auto nci = *(NCost++);
   if( ( ! std::isnan( C[ i ] ) ) && ( nci != C[ i ] ) )
    pairs.push_back( std::make_pair( i , nci ) );
   }
  std::sort( pairs.begin() , pairs.end() ,
	     []( auto & a , auto & b ) { return( a.first < b.first ); } );
  NC.resize( pairs.size() );
  auto NCit = NC.begin();
  auto nmsit = nms.begin();
  for( auto & el : pairs ) {
   *(nmsit++) = el.first;
   *(NCit++) = el.second;
   }
  nms.resize( pairs.size() );
  }

 if( nms.empty() )  // nothing left to change
  return;           // cowardly (and silently) return

 if( not_dry_run( issueAMod ) && ( AR & HasObj ) ) {
  // change abstract and physical representation together - - - - - - - - - -
  // in the meantime, if so instructed also issue abstract Modification
  copyidx( C , nms , NC.begin() );

  // note that modify_coefficients owns both vectors, so two copies have
  // to be made
  get_lfo()->modify_coefficients( std::move( NC ) , Subset( nms ) , true ,
				  un_ModBlock( issueAMod ) );
  }
 else
  // only change the physical representation- - - - - - - - - - - - - - - - -
  if( not_dry_run( issueMod ) )
   copyidx( C , nms , NC.begin() );

 f_cond_lower = dNAN;  // reset conditional bounds

 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared< SingleFlowDCRBlockSbstMod >( this ,
                                  SingleFlowDCRBlockMod::eChgCost , std::move( nms ) ) ,
			   Observer::par2chnl( issueMod ) );
 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end( SingleFlowDCRBlock::chg_costs( subset ) )

/*--------------------------------------------------------------------------*/

void SingleFlowDCRBlock::chg_cost( CNumber NCost , Index arc , 
			 ModParam issueMod , ModParam issueAMod )
{
 //if( arc >= get_NArcs() )
  //throw( std::invalid_argument( "invalid arc name" ) );

 if( arc < get_NArcs() ){
 if( std::isnan( C[ arc ] ) || ( C[ arc ] == NCost ) )
  return;  // if the arc is deleted or the cost does not change, all done

 f_cond_lower = dNAN;  // reset conditional bounds

 if( not_dry_run( issueAMod ) && ( AR & HasObj ) ) {
  // change abstract and physical representation together - - - - - - - - - -
  // in the meantime, if so instructed also issue abstract Modification
  C[ arc ] = NCost;

  get_lfo()->modify_coefficient( arc , NCost , un_ModBlock( issueAMod ) );
  }
 else
  // only change the physical representation- - - - - - - - - - - - - - - - -
  if( not_dry_run( issueMod ) )
   C[ arc ] = NCost;

 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared< SingleFlowDCRBlockRngdMod >( this ,
			   SingleFlowDCRBlockMod::eChgCost , Range( arc , arc + 1 ) ) ,
			   Observer::par2chnl( issueMod ) );
 #if CHECK_DS
  CheckAbsVSPhys();
 #endif
  }

 }  // end( SingleFlowDCRBlock::chg_cost )

/*--------------------------------------------------------------------------*/

void SingleFlowDCRBlock::chg_ucaps( c_Vec_FNumber_it NCap , Range rng ,
			  ModParam issueMod , ModParam issueAMod )
{
 rng.second = std::min( rng.second , get_NArcs() );
 if( rng.second <= rng.first )  // nothing to change
  return;                       // cowardly (and silently) return

 if( U.empty() ) {
  if( std::all_of( NCap , NCap + ( rng.second - rng.first ) ,
		   []( c_FNumber cap ) { return( cap >= Inf< FNumber >() ); }
		   ) )
   return;

  U.assign( get_MaxNArcs() , Inf< FNumber >() );
  }

 // check to see how many of the initial arcs are either deleted or not
 // really changing the capacity
 while( ( std::isnan( C[ rng.first ] ) || ( *NCap == U[ rng.first ] ) )
	&& ( rng.first < rng.second ) ) {
  ++rng.first;
  ++NCap;
  }

 if( rng.second <= rng.first )  // nothing left to change
  return;                       // cowardly (and silently) return

 // check to see how many of the final arcs are either deleted or not
 // really changing the capacity
 auto NCEit = NCap + ( rng.second - rng.first );
 while( ( ( std::isnan( C[ rng.second - 1 ] ) ) ||
	  ( *(--NCEit) == U[ rng.second - 1 ] ) )
	&& ( rng.first < rng.second ) )
  --rng.second;

 if( rng.second <= rng.first )  // nothing left to change
  return;                       // cowardly (and silently) return

 f_cond_lower = dNAN;  // reset conditional bounds

 if( not_dry_run( issueAMod ) && ( AR & HasFlw ) ) {
  // change abstract and physical representation together - - - - - - - - - -
  // in the meantime, if so instructed also issue abstract Modification

  if( ! ( AR & HasBnd ) )
   throw( std::logic_error(
		"bound constraints not defined, cannot change capacity" ) );

  not_ModBlock( issueAMod );
  auto ampar = open_if_needed( issueAMod , rng.second - rng.first );

  Index i = rng.first;

  // static part
  for( ; i < std::min( rng.second , get_NArcs() ) ;  ++i , ++NCap )
   if( ( ! std::isnan( C[ i ] ) ) && ( U[ i ] != *NCap ) ) {
    U[ i ] = *NCap;
    UB[ i ].set_rhs( *NCap , ampar );
    }

  close_if_needed( ampar , rng.second - rng.first );
  }
 else
  // only change the physical representation- - - - - - - - - - - - - - - - -
  // note that this also changes the capacity of deleted arcs, but since that
  // is never really used, it does not matter
  if( not_dry_run( issueMod ) )
   std::copy( NCap , NCap + ( rng.second - rng.first ) ,
	      U.begin() + rng.first );

 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared< SingleFlowDCRBlockRngdMod >( this ,
				              SingleFlowDCRBlockMod::eChgCaps , rng ) ,
			   Observer::par2chnl( issueMod ) );
 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end( SingleFlowDCRBlock::chg_ucaps( range ) )

/*--------------------------------------------------------------------------*/

void SingleFlowDCRBlock::chg_ucaps( c_Vec_FNumber_it NCap , Subset && nms ,
			  bool ordered  ,
			  ModParam issueMod , ModParam issueAMod )
{
 if( U.empty() ) {
  if( std::all_of( NCap , NCap + nms.size() ,
		   []( c_FNumber cap ) { return( cap >= Inf< FNumber >() ); }
		   ) )
   return;

  U.assign( get_MaxNArcs() , Inf< FNumber >() );
  }

 // eliminate from NCap and nms the entries corresponding to either
 // deleted arcs or arcs whose capacity actually does not change;
 // meanwhile, if nms is not ordered, order it
 Vec_FNumber NC;
 if( ordered ) {
  NC.resize( nms.size() );
  auto NCit = NC.begin();
  auto nmsit = nms.begin();
  for( auto i : nms ) {
   auto nci = *(NCap++);
   if( ( ! std::isnan( C[ i ] ) ) && ( nci != U[ i ] ) ) {
    *(nmsit++) = i;
    *(NCit++) = nci;
    }
   }
  nms.resize( std::distance( nms.begin() , nmsit ) );
  }
 else {
  using TP = std::pair< Index , FNumber >;
  std::vector< TP > pairs;
  pairs.reserve( nms.size() );
  for( auto i : nms ) {
   auto nci = *(NCap++);
   if( ( ! std::isnan( C[ i ] ) ) && ( nci != U[ i ] ) )
    pairs.push_back( std::make_pair( i , nci ) );
   }
  std::sort( pairs.begin() , pairs.end() ,
	     []( auto & a , auto & b ) { return( a.first < b.first ); } );
  NC.resize( pairs.size() );
  auto NCit = NC.begin();
  auto nmsit = nms.begin();
  for( auto & el : pairs ) {
   *(nmsit++) = el.first;
   *(NCit++) = el.second;
   }
  nms.resize( pairs.size() );
  }

 if( nms.empty() )  // nothing left to change
  return;           // cowardly (and silently) return

 f_cond_lower = dNAN;  // reset conditional bounds

 if( not_dry_run( issueAMod ) && ( AR & HasFlw ) ) {
  // change abstract and physical representation together - - - - - - - - - -
  // in the meantime, if so instructed also issue abstract Modification

  if( ! ( AR & HasBnd ) )
   throw( std::logic_error(
		"bound constraints not defined, cannot change capacity" ) );

  not_ModBlock( issueAMod );
  auto ampar = open_if_needed( issueAMod , nms.size() );

  // static part
  auto nit = nms.begin();
  auto NCit = NC.begin();
  for( ; ( nit != nms.end() ) && ( *nit < get_NArcs() ) ;
       ++NCit , ++nit ) {
   if( ( ! std::isnan( C[ *nit ] ) ) && ( U[ *nit ] != *NCit ) ) {
    U[ *nit ] = *NCit;
    UB[ *nit ].set_rhs( *NCit , ampar );
    }
   }

  close_if_needed( ampar , nms.size() );
  }
 else
  // only change the physical representation- - - - - - - - - - - - - - - - -
  // note that this also changes the capacity of deleted arcs, but since that
  // is never really used, it does not matter
  if( not_dry_run( issueMod ) )
   copyidx( U , nms , NC.begin() );

 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared< SingleFlowDCRBlockSbstMod >( this ,
                                  SingleFlowDCRBlockMod::eChgCaps , std::move( nms ) ) ,
			   Observer::par2chnl( issueMod ) );
 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end( SingleFlowDCRBlock::chg_ucaps( subset ) )

/*--------------------------------------------------------------------------*/

void SingleFlowDCRBlock::chg_ucap( FNumber NCap , Index arc ,
			 ModParam issueMod , ModParam issueAMod )
{
 if( arc >= get_NArcs() )
  throw( std::invalid_argument( "invalid arc name" ) );

 if( U.empty() ) {
  if( NCap < Inf< FNumber >() )
   return;

  U.assign( get_MaxNArcs() , Inf< FNumber >() );
  }

 if( U[ arc ] == NCap )
  return;

 f_cond_lower = dNAN;  // reset conditional bounds

 if( not_dry_run( issueMod ) )
  U[ arc ] = NCap;  // only change the physical representation - - - - - - -

 if( not_dry_run( issueAMod ) && ( AR & HasFlw ) ) {
  // change the abstract representation - - - - - - - - - - - - - - - - - - -
  // in the meantime, if so instructed also issue abstract Modification

  if( ! ( AR & HasBnd ) )
   throw( std::logic_error(
		"bound constraints not defined, cannot change capacity" ) );

   UB[ arc ].set_rhs( NCap , issueAMod );
  }

 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared< SingleFlowDCRBlockRngdMod >( this ,
			   SingleFlowDCRBlockMod::eChgCaps , Range( arc , arc + 1 )  ) ,
			   Observer::par2chnl( issueMod ) );
 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end( SingleFlowDCRBlock::chg_ucap )

/*--------------------------------------------------------------------------*/


void SingleFlowDCRBlock::chg_dfcts( c_Vec_CNumber_it NDfct , Range rng ,
			  ModParam issueMod , ModParam issueAMod )
{
 rng.second = std::min( rng.second , get_NNodes() );
 if( rng.second <= rng.first )  // nothing to change
  return;                 // cowardly (and silently) return

 if( B.empty() ) {
  if( std::all_of( NDfct , NDfct + ( rng.second - rng.first ) ,
		   []( c_FNumber dfct ) { return( dfct == 0 ); } ) )
   return;

  B.assign( get_MaxNNodes() , 0 );
  }

 c_Index ndiff = countdiff( NDfct , NDfct + ( rng.second - rng.first ) ,
			    B.cbegin() + rng.first );
 if( ! ndiff )
  return;

 if( not_dry_run( issueAMod ) && ( AR & HasFlw ) ) {
  // change abstract and physical representation together - - - - - - - - - -
  // in the meantime, if so instructed also issue abstract Modification

  not_ModBlock( issueAMod );
  auto ampar = open_if_needed( issueAMod , ndiff );

  Index i = rng.first;

  // static part
  for( ; i < std::min( rng.second , get_NNodes() ) ;  ++i , ++NDfct )
   if( B[ i ] != *NDfct ) {
    B[ i ] = *NDfct;
    E[ i ].set_both( *NDfct , ampar );
    }

  close_if_needed( ampar , ndiff );
  }
 else
  // only change the physical representation- - - - - - - - - - - - - - - - -
  if( not_dry_run( issueMod ) )
   std::copy( NDfct , NDfct + ( rng.second - rng.first ) ,
	      B.begin() + rng.first );

 // TODO: if some changes are "fake", restrict the range

 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared< SingleFlowDCRBlockRngdMod >( this ,
				              SingleFlowDCRBlockMod::eChgDfct , rng ) ,
			   Observer::par2chnl( issueMod ) );
 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end( SingleFlowDCRBlock::chg_dfcts( range ) )

/*--------------------------------------------------------------------------*/

void SingleFlowDCRBlock::chg_dfcts( c_Vec_CNumber_it NDfct , Subset && nms ,
			  bool ordered ,
			  ModParam issueMod , ModParam issueAMod )
{
 if( B.empty() ) {
  if( std::all_of( NDfct , NDfct + nms.size() ,
		   []( c_FNumber dfct ) { return( dfct == 0 ); } ) )
   return;

  B.assign( get_MaxNNodes() , 0 );
  }

 Index ndiff = countdiff( B , nms , NDfct , get_NNodes() );
 if( ! ndiff )
  return;

 if( not_dry_run( issueAMod ) && ( AR & HasFlw ) ) {
  // change abstract and physical representation together - - - - - - - - - -
  // in the meantime, if so instructed also issue abstract Modification

  not_ModBlock( issueAMod );
  auto ampar = open_if_needed( issueAMod , ndiff );
  }
 else
  // only change the physical representation- - - - - - - - - - - - - - - - -
  if( not_dry_run( issueMod ) )
   copyidx( B , nms , NDfct );

 // TODO: eliminate from nms the "fake" changes

 if( issue_pmod( issueMod ) ) {  // issue "physical Modification" - - - - - -
  // ensure the names are ordered even if they were not so originally
  if( ! ordered )
   std::sort( nms.begin() , nms.end() );

  Block::add_Modification( std::make_shared< SingleFlowDCRBlockSbstMod >( this ,
                                 SingleFlowDCRBlockMod::eChgDfct , std::move( nms ) ) ,
			   Observer::par2chnl( issueMod ) );
  }

 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end( SingleFlowDCRBlock::chg_dfcts( subset ) )

/*--------------------------------------------------------------------------*/

void SingleFlowDCRBlock::chg_dfct( FNumber NDfct , Index nde ,
			 ModParam issueMod , ModParam issueAMod )
{
 if( nde >= get_NNodes() )
  throw( std::invalid_argument( "invalid node name" ) );

 if( B.empty() ) {
  if( ! NDfct )
   return;

  B.assign( get_NNodes() , 0 );
  }

 if( B[ nde ] == NDfct )
  return;

 if( not_dry_run( issueMod ) )
  B[ nde ] = NDfct;  // change the physical representation- - - - - - - - - -

 if( not_dry_run( issueAMod ) && ( AR & HasFlw ) ) {
  // change the abstract representation - - - - - - - - - - - - - - - - - - -
  // in the meantime, if so instructed also issue abstract Modification
  E[ nde ].set_both( NDfct , issueAMod );
  }
 
 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared< SingleFlowDCRBlockRngdMod >( this ,
			   SingleFlowDCRBlockMod::eChgDfct , Range( nde , nde + 1 ) ) ,
			   Observer::par2chnl( issueMod ) );
 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end( SingleFlowDCRBlock::chg_dfct )

 /*--------------------------------------------------------------------------*/

void SingleFlowDCRBlock::chg_st( Index ns , Index nt , ModParam issueMod , 
        ModParam issueAMod )
{
 if( ns >= get_NNodes() || nt >= get_NNodes() )
  throw( std::invalid_argument( "invalid node name" ) );

 for( Index i = 0 ; i < get_NNodes() ; ++i )
  chg_dfct( 0 , i , issueMod, issueAMod );

 chg_dfct( -1 , ns , issueMod, issueAMod );
 chg_dfct( 1 , nt , issueMod, issueAMod );

 }  // end( SingleFlowDCRBlock::chg_st )

/*--------------------------------------------------------------------------*/

void SingleFlowDCRBlock::close_arcs( Range rng ,
			   ModParam issueMod , ModParam issueAMod )
{
 rng.second = std::min( rng.second , get_NArcs() );
 if( rng.second <= rng.first )  // nothing to change
  return;                       // cowardly (and silently) return

 // since the physical and abstract representation are the same, anything
 // that has to do with the abstract representation is skipped in the
 // "dry run" case; but the "physical Modification" is issued anyway

 if( not_dry_run( issueAMod ) ) {
  std::vector< ColVariable * > toclose;
  toclose.reserve( rng.second - rng.first );
  Index i = rng.first;

  // static part
  for( ; i < std::min( rng.second , get_NArcs() ) ; ++i )
   if( ! x[ i ].is_fixed() )
    toclose.push_back( & x[ i ] );

  if( toclose.empty() )
   return;

  // the physical and abstract representation are the same- - - - - - - - - -
  // change both (doh!), and if so instructed also issue abstract Modification

  not_ModBlock( issueAMod );
  auto ampar = open_if_needed( issueAMod , toclose.size() );

  for( auto vi : toclose ) {
   vi->set_value( 0 );
   vi->is_fixed( true , ampar );
   }

  close_if_needed( ampar , toclose.size() );
  }

 f_cond_lower = dNAN;  // reset conditional bounds

 // TODO: if some changes are "fake", restrict the range

 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared< SingleFlowDCRBlockRngdMod >( this ,
				             SingleFlowDCRBlockMod::eCloseArc , rng ) ,
			   Observer::par2chnl( issueMod ) );
 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end( SingleFlowDCRBlock::close_arcs( range ) )

/*--------------------------------------------------------------------------*/

void SingleFlowDCRBlock::close_arcs( Subset && nms , bool ordered  ,
			   ModParam issueMod , ModParam issueAMod )
{
 if( nms.empty() )
  return;

 // ensure the names are ordered even if they were not so originally
 if( ! ordered )
  std::sort( nms.begin() , nms.end() );

 if( nms.back() >= get_NArcs() )
  throw( std::invalid_argument( "invalid arc name" ) );

 // since the physical and abstract representation are the same, anything
 // that has to do with the abstract representation is skipped in the
 // "dry run" case; but the "physical Modification" is issued anyway

 if( not_dry_run( issueAMod ) ) {
  std::vector< ColVariable * > toclose;
  toclose.reserve( nms.size() );

  // static part
  auto nit = nms.begin();
  for( ; ( nit != nms.end() ) && ( *nit < get_NArcs() ) ; ++nit )
   if( ! x[ *nit ].is_fixed() )
    toclose.push_back( & x[ *nit ] );

  if( toclose.empty() )
   return;

  // the physical and abstract representation are the same- - - - - - - - - -
  // change both (doh!), and if so instructed also issue abstract Modification

  not_ModBlock( issueAMod );
  auto ampar = open_if_needed( issueAMod , toclose.size() );

  for( auto vi : toclose ) {
   vi->set_value( 0 );
   vi->is_fixed( true , ampar );
   }

  close_if_needed( ampar , toclose.size() );
  }

 f_cond_lower = dNAN;  // reset conditional bounds

 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared< SingleFlowDCRBlockSbstMod >( this ,
                                 SingleFlowDCRBlockMod::eCloseArc , std::move( nms ) ) ,
			   Observer::par2chnl( issueMod ) );
 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end( SingleFlowDCRBlock::close_arcs( subset ) )

/*--------------------------------------------------------------------------*/

void SingleFlowDCRBlock::close_arc( Index arc ,
			  ModParam issueMod , ModParam issueAMod )
{
 if( arc >= get_NArcs() )
  throw( std::invalid_argument( "invalid arc name" ) );

 // since the physical and abstract representation are the same, anything
 // that has to do with the abstract representation is skipped in the
 // "dry run" case; but the "physical Modification" is issued anyway

 if( not_dry_run( issueAMod ) ) {
  auto xa = i2p_x( arc );

  if( xa->is_fixed() )
   return;

  xa->set_value( 0 );

  // the physical and abstract representation are the same- - - - - - - - - -
  // change both (doh!), and if so instructed also issue abstract Modification

  xa->is_fixed( true , un_ModBlock( issueAMod ) );
  }

 f_cond_lower = dNAN;  // reset conditional bounds

 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared< SingleFlowDCRBlockRngdMod >( this ,
			   SingleFlowDCRBlockMod::eCloseArc , Range( arc , arc + 1 ) ) ,
			   Observer::par2chnl( issueMod ) );
 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end( SingleFlowDCRBlock::close_arc )

/*--------------------------------------------------------------------------*/

void SingleFlowDCRBlock::open_arcs( Range rng ,
			  ModParam issueMod , ModParam issueAMod )
{
 rng.second = std::min( rng.second , get_NArcs() );
 if( rng.second <= rng.first )  // nothing to change
  return;                       // cowardly (and silently) return

 // since the physical and abstract representation are the same, anything
 // that has to do with the abstract representation is skipped in the
 // "dry run" case; but the "physical Modification" is issued anyway

 if( not_dry_run( issueAMod ) ) {
  std::vector< ColVariable * > toopen;
  toopen.reserve( rng.second - rng.first );
  Index i = rng.first;

   // static part
  for( ; i < std::min( rng.second , get_NArcs() ) ; ++i )
   if( x[ i ].is_fixed() && ( ! std::isnan( C[ i ] ) ) )
    toopen.push_back( & x[ i ] );

  if( toopen.empty() )
   return;

  // the physical and abstract representation are the same- - - - - - - - - -
  // change both (doh!), and if so instructed also issue abstract Modification

  not_ModBlock( issueAMod );
  auto ampar = open_if_needed( issueAMod , toopen.size() );

  for( auto vi : toopen )
   vi->is_fixed( false , ampar );

  close_if_needed( ampar , toopen.size() );
  }

 f_cond_lower = dNAN;  // reset conditional bounds

 // TODO: if some changes are "fake", restrict the range

 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared< SingleFlowDCRBlockRngdMod >( this ,
				             SingleFlowDCRBlockMod::eOpenArc , rng ) ,
			   Observer::par2chnl( issueMod ) );
 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end( SingleFlowDCRBlock::open_arcs( range ) )

/*--------------------------------------------------------------------------*/

void SingleFlowDCRBlock::open_arcs( Subset && nms , bool ordered  ,
			  ModParam issueMod , ModParam issueAMod )
{
 if( nms.empty() )
  return;

 // ensure the names are ordered even if they were not so originally
 if( ! ordered )
  std::sort( nms.begin() , nms.end() );

 if( nms.back() >= get_NArcs() )
  throw( std::invalid_argument( "invalid arc name" ) );

 // since the physical and abstract representation are the same, anything
 // that has to do with the abstract representation is skipped in the
 // "dry run" case; but the "physical Modification" is issued anyway

 if( not_dry_run( issueAMod ) ) {
  std::vector< ColVariable * > toopen;
  toopen.reserve( nms.size() );

  // static part
  auto nit = nms.begin();
  for( ; ( nit != nms.end() ) && ( *nit < get_NArcs() ) ; ++nit )
   if( x[ *nit ].is_fixed() && ( ! std::isnan( C[ *nit ] ) ) )
    toopen.push_back( & x[ *nit ] );

  if( toopen.empty() )
   return;

  // the physical and abstract representation are the same- - - - - - - - - -
  // change both (doh!), and if so instructed also issue abstract Modification

  not_ModBlock( issueAMod );
  auto ampar = open_if_needed( issueAMod , toopen.size() );

  for( auto vi : toopen )
   vi->is_fixed( false , ampar );

  close_if_needed( ampar , toopen.size() );
  }

 f_cond_lower = dNAN;  // reset conditional bounds

 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared< SingleFlowDCRBlockSbstMod >( this ,
                                  SingleFlowDCRBlockMod::eOpenArc , std::move( nms ) ) ,
			   Observer::par2chnl( issueMod ) );
 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end( SingleFlowDCRBlock::open_arcs( subset ) )

/*--------------------------------------------------------------------------*/

void SingleFlowDCRBlock::open_arc( Index arc ,
			 ModParam issueMod , ModParam issueAMod )
{
 if( arc >= get_NArcs() )
  throw( std::invalid_argument( "invalid arc name" ) );

 // since the physical and abstract representation are the same, anything
 // that has to do with the abstract representation is skipped in the
 // "dry run" case; but the "physical Modification" is issued anyway

 if( not_dry_run( issueAMod ) ) {
  auto xa = i2p_x( arc );

  if( ( ! xa->is_fixed() ) || std::isnan( C[ arc ] ) )
   return;

  // the physical and abstract representation are the same- - - - - - - - - -
  // change both (doh!), and if so instructed also issue abstract Modification

  xa->is_fixed( false , un_ModBlock( issueAMod ) );
  }

 f_cond_lower = dNAN;  // reset conditional bounds

 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared< SingleFlowDCRBlockRngdMod >( this ,
			   SingleFlowDCRBlockMod::eOpenArc , Range( arc , arc + 1 ) ) ,
			   Observer::par2chnl( issueMod ) );
 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end( SingleFlowDCRBlock::open_arc )

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

void SingleFlowDCRBlock::guts_of_destructor( void )
{
 // clear() all Constraint to ensure that they do not bother to un-register
 // themselves from Variable that are going to be deleted anyway

 // clear the bound constraints
 Constraint::clear( UB );   // static

 // clear the flow conservation constraints
 Constraint::clear( E );   // static

 DCR_cnst.clear();
 Constraint::clear( Indicator_cnst_rmin ); 
 Constraint::clear( Indicator_cnst_r1 );
 Constraint::clear( Indicator_cnst_r2 );
 Constraint::clear( PC_cuts );
 Constraint::clear( PC_cuts_min );
 
 cone_min_cnst.clear();
 Constraint::clear( cone_cnst );

 c.clear();  // clear the Objective

 // delete all Variable
 x.clear();   // static
 r.clear();
 theta.clear();

 U.clear();
 C.clear();
 B.clear();
 NodeDelays.clear();

 // explicitly reset all Constraint and Variable
 // this is done for the case where this method is called prior to re-loading
 // a new instance: if not, the new representation would be added to the
 // (no longer current) one
 reset_static_constraints();
 reset_static_variables();
 reset_dynamic_constraints();
 reset_dynamic_variables();
 reset_objective();

 AR = 0;

 }  // end( SingleFlowDCRBlock::guts_of_destructor )

/*--------------------------------------------------------------------------*/

void SingleFlowDCRBlock::guts_of_add_Modification( p_Mod mod , ChnlName chnl )
{
 // process abstract Modification - - - - - - - - - - - - - - - - - - - - - -
 /* This requires to patiently sift through the possible Modification types
  * to find what this Modification exactly is and appropriately mirror the
  * changes to the "abstract representation" to the "physical one".
  *
  * Note that since SingleFlowDCRBlock is a "leaf" Block (has no sub-Block), this
  * method does not have to deal with GroupModification since these are
  * produced by Block::add_Modification(), but this method is called
  * *before* that one is.
  *
  * As an important consequence,
  *
  *   THE STATE OF THE DATA STRUCTURE IN SingleFlowDCRBlock WHEN THIS METHOD IS
  *   EXECUTED IS PRECISELY THE ONE IN WHICH THE Modification WAS ISSUED:
  *   NO COMPLICATED OPERATIONS (Variable AND/OR Constraint BEING
  *   ADDED/REMOVED ...) CAN HAVE BEEN PERFORMED IN THE MEANTIME
  *
  * This assumption drastically simplifies some of the logic here.*/

   if( const auto tmod = dynamic_cast< C05FunctionModVarsAddd * >( mod ) ) {
    if( ! ( AR & HasObj ) )
      throw( std::invalid_argument( "Modification to non-constructed Objective"
            ) );

    auto lfo = static_cast< LinearFunction * const >( tmod->function() );

    if( static_cast< LinearFunction * const >( c.get_function() ) != lfo ){
      throw( std::invalid_argument( "Modification to non-Objective" ) );
      return;
    }

    if( AR & HasObj){    
    // add the new coefficient in the objective
      const auto mp = Observer::make_par( eNoBlck , chnl );
      LinearFunction::v_coeff_pair p( get_NArcs() );
      for( Index i = 0 ; i < get_NArcs() ; ++i ) {
        ///!
      }
    }

  return;
  }

 // C05FunctionModLinRngd - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( const auto tmod = dynamic_cast< C05FunctionModLinRngd * >( mod ) ) {
  if( ! ( AR & HasObj ) )
   throw( std::invalid_argument( "Modification to non-constructed Objective"
				 ) );

  auto lfo = static_cast< LinearFunction * const >( tmod->function() );
  if( static_cast< LinearFunction * const >( c.get_function() ) != lfo )
   throw( std::invalid_argument( "Modification to non-Objective" ) );

  // note: in the following we can assume that the Range in tmod is
  //       precisely the one we have to use since no Variable can have
  //       been added or deleted, which saves *a lot* of trouble

  if( tmod->range().second == tmod->range().first + 1 )
   // changing one cost only
   chg_cost( lfo->get_coefficient( tmod->range().first ) ,
	     tmod->range().first , make_par( eNoBlck , chnl ) , eDryRun );
  else {                            // changing many costs at once
   Vec_CNumber NC( tmod->range().second - tmod->range().first );
   auto NCit = NC.begin();
   for( Index i = tmod->range().first ; i < tmod->range().second ; ){
    *(NCit++) = lfo->get_coefficient( i++ );
   }

   chg_costs( NC.begin() , tmod->range() ,
	      make_par( eNoBlck , chnl ) , eDryRun );
   }

  return;
  }

 // C05FunctionModLinSbst - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( const auto tmod = dynamic_cast< C05FunctionModLinSbst * >( mod ) ) {
  if( ! ( AR & HasObj ) )
   throw( std::invalid_argument( "Modification to non-constructed Objective"
				 ) );

  auto lfo = static_cast< LinearFunction * const >( tmod->function() );
  if( static_cast< LinearFunction * const >( c.get_function() ) != lfo )
   throw( std::invalid_argument( "Modification to non-Objective" ) );

  // note: in the following we can assume that the Subset in tmod is
  //       precisely the one we have to use since no Variable can have
  //       been added or deleted, which saves *a lot* of trouble
  // note: chg_costs() owns subset, so a copy has to be made

  Vec_CNumber NC( tmod->subset().size() );
  auto NCit = NC.begin();
  for( auto i : tmod->subset() )
   *(NCit++) = lfo->get_coefficient( i++ );

  chg_costs( NC.begin() , Subset( tmod->subset() ) , true ,
	     make_par( eNoBlck , chnl ) , eDryRun );
  return;
  }

 // RowConstraintMod- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( const auto tmod = dynamic_cast< RowConstraintMod * >( mod ) ) {
  if( ! ( AR & HasFlw ) )
   throw( std::invalid_argument(
			    "Modification to non-constructed Constraint" ) );

  if( tmod->type() == RowConstraintMod::eChgRHS ) {
   auto cp = dynamic_cast< LB0Constraint * const >( tmod->constraint() );
   if( ! cp )
    throw( std::invalid_argument( "invalid Modification to Constraint" ) );

   chg_ucap( cp->get_rhs() , p2i_ub( cp ) ,
	     make_par( eNoBlck , chnl ) , eDryRun );
   return;
   }

  if( tmod->type() == RowConstraintMod::eChgBTS ) {
   auto cp = static_cast< FRowConstraint * const >( tmod->constraint() );
   if( ! cp )
    throw( std::invalid_argument( "invalid Modification to Constraint" ) );

   chg_dfct( cp->get_rhs() , p2i_e( cp ) ,
	     make_par( eNoBlck , chnl ) , eDryRun );
   return;
   }

  if( tmod->type() == RowConstraintMod::eChgLHS ) {
   std::cout << "LHS";
   return;
   }

  if( tmod->type() == RowConstraintMod::eChgBTS ) {
   std::cout << "BTS";  
   return;
   }

  if( tmod->type() == RowConstraintMod::eRowConstModLastParam ) {
   return;
   }

  throw( std::invalid_argument( "illegal Modification to Constraint" ) );
  }

 // VariableMod - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( const auto tmod = dynamic_cast< VariableMod * >( mod ) ) {
  auto xi = dynamic_cast< const ColVariable * >( tmod->variable() );
  if( ! xi )
   throw( std::logic_error( "Modification to wrong type of Variable" ) );
  /*
  if( ( xi->get_type() != ColVariable::kNonNegative ) &&
      ( xi->get_type() != ColVariable::kNatural ) )
   throw( std::logic_error( "changing type of flow Variable not allowed" ) );
  */
   
  auto i = p2i_x( xi );
  if( std::isnan( C[ i ] ) )
   throw( std::logic_error( "[un]fixing deleted arc not allowed" ) );
  if( xi->is_fixed() )
   close_arc( i , make_par( eNoBlck , chnl ) , eDryRun );
  else
   open_arc( i , make_par( eNoBlck , chnl ) , eDryRun );

  return;
  }

 //////////////throw( std::invalid_argument( "unsupported Modification to SingleFlowDCRBlock" ) );

 }  // end( SingleFlowDCRBlock::guts_of_add_Modification )

/*--------------------------------------------------------------------------*/

void SingleFlowDCRBlock::compute_conditional_bounds( void )
{
 f_cond_lower = f_cond_upper = 0;

 auto tC = C.begin();
 auto tU = U.begin();

 for( ; tC < C.end() ; ++tC , ++tU ) {
  if( *tC == 0 )
   continue;

  if( *tC < 0 ) {
   if( *tU == Inf< FNumber >() ) {
    f_cond_lower = -Inf< double >();
    break;
    }
   else
    f_cond_lower += *tC * (*tU);
   }
  else
   if( *tU == Inf< FNumber >() ) {
    f_cond_upper = Inf< double >();
    break;
    }
   else
    f_cond_upper += *tC * (*tU);
   }

 if( f_cond_lower > -Inf< double >() ) {
  for( ; tC < C.end() ; ++tC , ++tU )
   if( *tC < 0 ) {
    if( *tU == Inf< FNumber >() ) {
     f_cond_lower = -Inf< double >();
     break;
     }
    else
     f_cond_lower += *tC * (*tU);
    }
  }

 if( f_cond_upper < Inf< double >() ) {
  for( ; tC < C.end() ; ++tC , ++tU )
   if( *tC > 0 ) {
    if( *tU == Inf< FNumber >() ) {
     f_cond_upper = Inf< double >();
     break;
     }
    else
     f_cond_upper += *tC * (*tU);
    }
  }
 }  // end( SingleFlowDCRBlock::compute_conditional_bounds )

/*--------------------------------------------------------------------------*/

#ifndef NDEBUG

void SingleFlowDCRBlock::CheckAbsVSPhys( void )
{
 // check that the (part that has actually been constructed of the) abstract
 // representation coincides with the physical representation

 // check variables, these are always there - - - - - - - - - - - - - - - - -
 if( x.size() != get_NArcs() )
  std::cerr << "x.size() != NArcs" << std::endl;

 // check flow constraints- - - - - - - - - - - - - - - - - - - - - - - - - -
 if( AR & HasFlw ) { 
  if( E.size() != get_NNodes() )
   std::cerr << "E.size() != NNodes" << std::endl;
  
  Subset NRIncid( get_NNodes() , 0 );

  Index objbnd = 0;          // number of actives between obj and bound
  if( AR & HasBnd )
   ++objbnd;
  if( AR & HasObj )
   ++objbnd;
  Index expnc = 2 + objbnd;  // total number of active stuff per variable

  // static arcs
  Index a = 0;
  for( auto xi = x.begin() ; a < get_NArcs() ; ++a , ++xi ) {
   if( is_deleted( a ) ) {
    std::cerr << "static arc " << a << " deleted" << std::endl;
    continue;
    }
   
   if( xi->get_num_active() != expnc )
    std::cerr << "arc " << a << " active in " << xi->get_num_active()
	      << " constraints" << std::endl;

   Index sna = SN[ a ];
   if( ! sna )
    std::cerr << "SN[ " << a << " ] == 0" << std::endl;
   else
    --sna;
   if( sna >= get_NArcs() )
    std::cerr << "SN[ " << a << " ] == " << sna << " >= |A|" << std::endl;

   if( ( SN[ a ] > 0 ) && ( sna < get_NArcs() ) ) {
    ++NRIncid[ sna ];
    auto snc = get_lfc( i2p_e( sna ) );
    auto sni = snc->is_active( &(*xi) );
    if( sni >= snc->get_num_active_var() )
     std::cerr << "static arc " << a
	       << " absent in flow constraint for SN[ a ] == "
	       << sna << std::endl;
    }
    
   Index ena = EN[ a ];
   if( ! ena )
    std::cerr << "EN[ " << a << " ] == 0" << std::endl;
   else
    --ena;
   if( ena >= get_NArcs() )
    std::cerr << "EN[ " << a << " ] == " << ena << " >= |A|" << std::endl;

   if( ( EN[ a ] > 0 ) && ( ena < get_NArcs() ) ) {
    ++NRIncid[ ena ];
    auto enc = get_lfc( i2p_e( ena ) );
    auto eni = enc->is_active( &(*xi) );
    if( eni >= enc->get_num_active_var() )
     std::cerr << "static arc " << a
	       << " absent in flow constraint for EN[ a ] == "
	       << ena << std::endl;
    }
   }

  // static nodes
  Index n = 0;
  for( auto ni = E.begin() ; n < get_NNodes() ; ++n , ++ni ) {
   auto lni = get_lfc( &(*ni) );
   if( lni->get_num_active_var() != NRIncid[ n ] )
    std::cerr << "active variables in static flow constraint " << n
	      << " == " << lni->get_num_active_var()
	      << " do not match with incident arcs " << NRIncid[ n ]
	      << std::endl;
   }
  }  // end( if( AR & HasFlw ) )

 // check bound constraints - - - - - - - - - - - - - - - - - - - - - - - - -
 if( AR & HasBnd ) {
  // static bounds
  Index a = 0;
  auto UBi = UB.begin();
  for( auto xi = x.begin() ; xi != x.end() ; ++a , ++xi , ++UBi ) {
   if( UBi->is_active( &(*xi) ) >= UBi->get_num_active_var() )
    std::cerr << "static arc " << a << " absent in bound constraint"
	      << std::endl;
   }
  }  // end( if( AR & HasBnd ) )

 // check objective - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( AR & HasObj ) {
  auto lfo = get_lfo();
  if( lfo->get_num_active_var() != get_NArcs() )
   std::cerr << "objective has " << lfo->get_num_active_var()
	     << " variables while |A| = " << get_NArcs() << std::endl;

  Index a = 0;
  for( auto & xi : x ) {
   if( lfo->is_active( & xi ) >= get_NArcs() )
    std::cerr << "static arc " << a << " absent from objective " << std::endl;
   ++a;
   }

  }  // end( if( AR & HasObj ) )
 }  // end( SingleFlowDCRBlock::CheckAbsVSPhys )

#endif

/*--------------------------------------------------------------------------*/
/*-------------------------- METHODS OF DCRSolution ------------------------*/
/*--------------------------------------------------------------------------*/

void DCRSolution::deserialize( const netCDF::NcGroup & group )
{
 netCDF::NcDim na = group.getDim( "NumArcs" );
 if( na.isNull() )
  v_x.clear();
 else {
  netCDF::NcVar fs = group.getVar( "ArcSolution" );
  if( fs.isNull() )
   v_x.clear();
  else {
   v_x.resize( na.getSize() );
   fs.getVar( v_x.data() );
   }
  }

 if( na.isNull() )
  v_r.clear();
 else {
  netCDF::NcVar fs = group.getVar( "FlowSolution" );
  if( fs.isNull() )
   v_r.clear();
  else {
   v_r.resize( na.getSize() );
   fs.getVar( v_r.data() );
   }
  }
 }  // end( DCRSolution::deserialize )

/*--------------------------------------------------------------------------*/

void DCRSolution::read( const Block * block )
{
 auto DCRB = dynamic_cast< const SingleFlowDCRBlock * >( block );
 if( ! DCRB )
  throw( std::invalid_argument( "block is not a SingleFlowDCRBlock" ) );

 // read flows- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( ! v_x.empty() ) {
  if( v_x.size() < DCRB->get_NArcs() )
   v_x.resize( DCRB->get_NArcs() );

  DCRB->get_x( v_x.begin() );
  }

 // read potentials - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( ! v_r.empty() ) {
  if( v_r.size() < DCRB->get_NArcs() )
   v_r.resize( DCRB->get_NArcs() );

  DCRB->get_r( v_r.begin() );
  }
 }  // end( DCRSolution::read )

/*--------------------------------------------------------------------------*/

void DCRSolution::write( Block * block ) 
{
 auto DCRB = dynamic_cast< SingleFlowDCRBlock * >( block );
 if( ! DCRB )
  throw( std::invalid_argument( "block is not a SingleFlowDCRBlock" ) );

 // write flows - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( ! v_x.empty() ) {
  if( v_x.size() < DCRB->get_NArcs() )
   throw( std::invalid_argument( "incompatible flow size" ) );

  DCRB->set_x( v_x.begin() );
  }

 if( ! v_r.empty() ) {
  if( v_r.size() < DCRB->get_NArcs() )
   throw( std::invalid_argument( "incompatible flow size" ) );

  DCRB->set_r( v_r.begin() );
  }
 }  // end( DCRSolution::write )

/*--------------------------------------------------------------------------*/

void DCRSolution::serialize( netCDF::NcGroup & group ) const
{
 // always call the method of the base class first
 Solution::serialize( group );

 std::vector< size_t > startp = { 0 };

 if( ! v_x.empty() ) {
  netCDF::NcDim na = group.addDim( "NumArcs" , v_x.size() );

  std::vector< size_t > countpa = { v_x.size() };

  ( group.addVar( "ArcSolution" , netCDF::NcDouble() , na ) ).putVar(
					      startp , countpa , v_x.data() );
  }

 if( ! v_r.empty() ) {
  netCDF::NcDim na = group.addDim( "NumArcs" , v_x.size() );

  std::vector< size_t > countpa = { v_r.size() };

  ( group.addVar( "FlowSolution" , netCDF::NcDouble() , na ) ).putVar(
					      startp , countpa , v_r.data() );
  }
 
 }  // end( DCRSolution::serialize )

/*--------------------------------------------------------------------------*/

DCRSolution * DCRSolution::scale( double factor ) const
{
 auto * sol = DCRSolution::clone( true );

 if( ! v_x.empty() )
  for( SingleFlowDCRBlock::Index i = 0 ; i < v_x.size() ; ++i )
   sol->v_x[ i ] = v_x[ i ];

 if( ! v_r.empty() )
  for( SingleFlowDCRBlock::Index i = 0 ; i < v_r.size() ; ++i )
   sol->v_r[ i ] = v_r[ i ] * factor;

 return( sol );

 }  // end( DCRSolution::scale )

/*--------------------------------------------------------------------------*/

void DCRSolution::sum( const Solution * solution , double multiplier )
{
 auto DCRS = dynamic_cast< const DCRSolution * >( solution );
 if( ! DCRS )
  throw( std::invalid_argument( "solution is not a DCRSolution" ) );

 if( ! v_x.empty() ) {
  if( v_x.size() != DCRS->v_x.size() )
   throw( std::invalid_argument( "incompatible flow size" ) );

  for( SingleFlowDCRBlock::Index i = 0 ; i < v_x.size() ; ++i )
   v_x[ i ] = std::max( v_x[ i ], DCRS->v_x[ i ] );
  }

 if( ! v_r.empty() ) {
  if( v_r.size() != DCRS->v_r.size()  )
   throw( std::invalid_argument( "incompatible reserve size" ) );

  for( SingleFlowDCRBlock::Index i = 0 ; i < v_r.size() ; ++i )
   v_r[ i ] += DCRS->v_r[ i ] * multiplier;
  }
 }  // end( DCRSolution::sum )

/*--------------------------------------------------------------------------*/

DCRSolution * DCRSolution::clone( bool empty ) const
{
 auto *sol = new DCRSolution();

 if( empty ) {
  if( ! v_x.empty() )
   sol->v_x.resize( v_x.size() );

  if( ! v_r.empty() )
   sol->v_r.resize( v_r.size() );
  }
 else {
  sol->v_x = v_x;
  sol->v_r = v_r;
  }

 return( sol );

 }  // end( DCRSolution::clone )

/*--------------------------------------------------------------------------*/
/*------------------- End File SingleFlowDCRBlock.cpp ----------------------*/
/*--------------------------------------------------------------------------*/
