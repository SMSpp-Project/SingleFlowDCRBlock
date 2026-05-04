/*--------------------------------------------------------------------------*/
/*-------------------- File SingleFlowDCRBlock.h ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the *concrete* class SingleFlowDCRBlock, which implements
 * the Block concept [see Block.h] for the solution of Delay-Constrained
 * Routing problems (DCR) relative to a single flow.
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
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __SingleFlowDCRBlock
 #define __SingleFlowDCRBlock  /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "Block.h"

#include "LinearFunction.h"

#include "QuadFunction.h"

#include "FRealObjective.h"

#include "FRowConstraint.h"

#include "OneVarConstraint.h"

#include "Solution.h"

/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{
 class SingleFlowDCRBlock;     // forward declaration of SingleFlowDCRBlock

 class DCRSolution;  // forward declaration of DCRSolution

/*--------------------------------------------------------------------------*/
/*----------------------- SingleFlowDCRBlock-RELATED TYPES -----------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup SingleFlowDCRBlock_TYPES SingleFlowDCRBlock-related types
 *  @{ */

 using p_SingleFlowDCRBlock = SingleFlowDCRBlock *;  ///< a pointer to SingleFlowDCRBlock

 using Vec_SingleFlowDCRBlock = std::vector< p_SingleFlowDCRBlock>;
 ///< a vector of pointers to SingleFlowDCRBlock

 using Vec_SingleFlowDCRBlock_it = Vec_SingleFlowDCRBlock::iterator;
 ///< iterator for a Vec_SingleFlowDCRBlock

 using c_Vec_SingleFlowDCRBlock = const Vec_SingleFlowDCRBlock;
 ///< a const vector of pointers to SingleFlowDCRBlock

 using c_Vec_SingleFlowDCRBlock_it = c_Vec_SingleFlowDCRBlock::iterator;
 ///< iterator for a c_Vec_SingleFlowDCRBlock

/** @}  end( group( SingleFlowDCRBlock_TYPES ) ) */ 
/*--------------------------------------------------------------------------*/
/*------------------------------- CLASSES ----------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup SingleFlowDCRBlock_CLASSES Classes in SingleFlowDCRBlock.h
 *  @{ */

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS SingleFlowDCRBlock ----------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// implementation of the Block concept for the (linear) Min-Cost Flow problem
/** The SingleFlowDCRBlock class implements the Block concept [see Block.h] 
 * for the Min-Cost Flow DCR problem.
 **/

class SingleFlowDCRBlock : public Block
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public types
 *
 * SingleFlowDCRBlock defines three main public types:
 @{ */

/*--------------------------------------------------------------------------*/

 typedef const double c_double;            ///< a read-only double
 typedef std::vector< double > Vec_double; ///< a vector of double
 typedef const Vec_double c_Vec_double;    ///< a const vector of double
 typedef Vec_double::iterator Vec_double_it;   ///< iterator in Vec_double
 typedef Vec_double::const_iterator c_Vec_double_it;
                                           ///< const iterator in Vec_double

/*--------------------------------------------------------------------------*/

/** @} ---------------------------------------------------------------------*/
/*------------------------------- FRIENDS ----------------------------------*/
/*--------------------------------------------------------------------------*/

 friend DCRSolution;  ///< make DCRSolution friend

/*--------------------------------------------------------------------------*/
/*--------------------- PUBLIC METHODS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 *  @{ */

 /// constructor of SingleFlowDCRBlock, taking a pointer to the father (generic) Block
 /** Constructor of SingleFlowDCRBlock. It accepts a pointer to the father Block, which
  * can be of any type, defaulting to nullptr so that this can also be used as
  * the void constructor. */

 explicit SingleFlowDCRBlock( Block *father = nullptr )
  : Block( father ) , NNodes( 0 ) , NArcs( 0 ) , MaxNNodes( 0 ) , AR( 0 ) ,
    f_cond_lower( - Inf< double >() ) , f_cond_upper( - Inf< double >() ) { }

/*--------------------------------------------------------------------------*/
 /// destructor of SingleFlowDCRBlock: deletes the abstract representation, if any

 virtual ~SingleFlowDCRBlock() { guts_of_destructor(); }

/** @} ---------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */

 /// loads the DCR instance from memory
 /** Loads the DCR instance from memory. The parameters are what you expect:
  *
  * - n    is the current number of nodes of the network
  *
  * - m    is the current number of arcs of the network
  *
  * - pSn  is the vector of the arc starting nodes, which must have size at
  *        least m
  *
  * - pEn  is the vector of the arc ending nodes, which must have size at
  *        least m
  *
  * - pU   is the vector of the arc upper capacities; capacities must be
  *        nonnegative, but can be infinite; it must either have size at
  *        least m or be empty, in the latter case all capacities are taken
  *        to be infinite
  *
  * - pC   is the vector of the arc costs; it must either have size at
  *        least m or be empty, in the latter case all arc costs are taken
  *        to be 0
  *
  * - pB   is the vector of the node deficits; source nodes have negative
  *        deficits and sink nodes have positive deficits; it must either
  *        have size at least n or be empty, in the latter case all deficits
  *        are taken to be 0 (a circulation problem)
  *
  * Like load( std::istream & ), if there is any Solver attached to this
  * SingleFlowDCRBlock then a NBModification (the "nuclear option") is issued. */

 void load( Index n , Index m , c_Subset & pEn , c_Subset & pSn ,
	    c_Vec_double & pU = {} , c_Vec_double & pC = {} ,
      c_Vec_double & pNodeDelays = {} , c_Vec_double & pLinkDelays = {} , 
      c_double FlowBursts  = 0 , c_double FlowDeadlines  = 0 ,
      c_double MTU  = 0 , c_double rho  = 0 );

/*--------------------------------------------------------------------------*/
 /// extends Block::deserialize( netCDF::NcGroup )
 /** Extends Block::deserialize( netCDF::NcGroup ) to the specific format of
  * a SingleFlowDCRBlock. Besides what is managed by the serialize() method of the base
  * Block class, the group should contain the following:
  *
  * - the dimension "NNodes" containing the number of nodes in the graph;
  *
  * - the dimension "NArcs" containing the number of arcs in the graph;
  *
  * - the variable "C", of type double and indexed over the dimension "NArcs";
  *   the i-th entry of the variable is assumed to contain the cost of the
  *   i-th arc in the graph, that whose start node and ending node are
  *   specified by the variable "SN" and "EN" (below);
  *
  * - the variable "U", of type double and indexed over the dimension "NArcs";
  *   the i-th entry of the variable is assumed to contain the upper capacity
  *   of the i-th arc in the graph (the lower capacity being fixed to 0),
  *   that whose start node and ending node are specified by the variable
  *   "SN" and "EN" (below);
  *
  * - the variable "B", of type double and indexed over the dimension
  *   "NNodes"; the i-th entry of the variable is assumed to contain the
  *   deficit of the i-th node in the graph (note that node names here go
  *   from 0 to NNodes.getSize() - 1);
  *
  * - the variable "SN", of type int and indexed over the dimension "NArcs";
  *   the i-th entry of the variable is assumed to contain the starting node
  *   of the i-th arc in the graph (note that node names here go from 1 to
  *   NNodes.getSize(), i.e., they are shifted by 1 w.r.t. to the indices
  *   used in the "B" variable);
  *
  * - the variable "EN", of type int and indexed over the dimension "NArcs";
  *   the i-th entry of the variable is assumed to contain the ending node
  *   of the i-th arc in the graph (note that node names here go from 1 to
  *   NNodes.getSize(), i.e., they are shifted by 1 w.r.t. to the indices
  *   used in the "B" variable).
  *
  * The two dimensions "NNodes" and "NArcs" are mandatory, such as are the
  * two variables "SN" and "EN". The three other variables are optional. If
  * "C" is missing, all arc costs are assumed to be 0. If "U" is missing, all
  * arc capacities are assumed to be infinite. If "B" is missing, all node
  * deficits are assumed to be 0. Finally, all the dimensions "DynNNodes",
  * "DynNArcs", "MaxDynNNodes" and "MaxDynNArcs" are optional: if they are
  * missing they are treated as being 0 (this happening for all four means
  * that the graph is "fully static" and cannot be changed). */

 void deserialize( const netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/
 /// loads the DCR instance from file in DIMACS standard format
 /** Protected method for loading a SingleFlowDCRBlock out of a std::istream (which is
  * what operator>> is dispatched to. The std::istream is assumed to contain
  * the description of a DCR instance in DIMACS standard format, which is 
  * the following. The first line must be
  *
  *      p min <number of nodes> <number of arcs>
  *
  * Then the node definition lines must be found, in the form
  *
  *      n <node number> <node supply>
  *
  * Not all nodes need have a node definition line; these are given zero
  * supply, i.e., they are transhipment nodes (supplies are the inverse of
  * deficits, i.e., a node with positive supply is a source node). Finally,
  * the arc definition lines must be found, in the form
  *
  *    a <start node> <end node> <lower bound> <upper bound> <flow cost>
  *
  * There must be exactly <number of arcs> arc definition lines in the file.
  *
  * Note that the file format accepted by load() is more general than the
  * DIMACS standard format, in that node and arc definitions can be mixed in
  * any order, while the DIMACS file requires all node information to appear
  * before all arc information. Also, capacities of arcs can be set to
  * +Inf< double >() by putting "INF", "Inf" or "inf" in the file (actually,
  * any string starting with "I" or "i" where these would be expected).
  *
  * Note that the graph as provided by this method is considered to be
  * "fully static".
  *
  * Since there is only one supported input format, \p frmt is ignored.
  *
  * Like load( memory ), if there is any Solver attached to this SingleFlowDCRBlock
  * then a NBModification (the "nuclear option") is issued. */

 void load( std::istream &input , char frmt = 0 ) override;

 void load_dcr( std::istream &input , Index NNodes, Index NArcs );

/*--------------------------------------------------------------------------*/
 /// generate the abstract variables of the DCR
 /** Method that generates the abstract Variable of the DCR. These are:
  */

 void generate_abstract_variables( Configuration *stvv = nullptr ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// generate the static constraint of the DCR
 /** Method that generates the abstract constraint of the DCR. These are:
  */
 
 void generate_abstract_constraints( Configuration *stcc = nullptr ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// generate the objective of the DCR
 /** Method that generates the objective of the DCR. Although this would seem
  * to be an exceedingly simple object, there is still a nontrivial decision
  * to be made about it, i.e., whether it is represented as a "sparse"
  * LinearFunction or a "dense" one. This is governed by objc: if
  */

 void generate_objective( Configuration *objc = nullptr ) override;

 void generate_dynamic_constraints( Configuration *stcc = nullptr ) override;

/** @} ---------------------------------------------------------------------*/
/*-------------- Methods for reading the data of the SingleFlowDCRBlock --------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for reading the data of the SingleFlowDCRBlock
 *  @{ */

 /// getting the current sense of the Objective, which is minimization

 [[nodiscard]] int get_objective_sense( void ) const override {
  return( Objective::eMin );
  }

/*--------------------------------------------------------------------------*/
 /// get the number of nodes

 [[nodiscard]] Index get_NNodes( void ) const { return( NNodes ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the number of arcs

 [[nodiscard]] Index get_NArcs( void ) const { return( NArcs ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the maximum number of nodes

 [[nodiscard]] Index get_MaxNNodes( void ) const { return( MaxNNodes ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the maximum number of arcs

 [[nodiscard]] Index get_MaxNArcs( void ) const { return( SN.size() ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the MTU

 [[nodiscard]] double get_MTU( void ) const { return( MTU ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the NodeDelays

 [[nodiscard]] Vec_double get_NodeDelays( void ) const { return( NodeDelays ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the NodeDelays

 [[nodiscard]] Vec_double get_LinkDelays( void ) const { return( LinkDelays ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the FlowBurst

 [[nodiscard]] double get_FlowBurst( void ) const { return( FlowBursts ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the FlowDeadline

 [[nodiscard]] double get_FlowDeadline( void ) const { return( FlowDeadlines ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the rho

 [[nodiscard]] double get_rho( void ) const { return( rho ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns true if there are static nodes (= possibly flow constraints)

 [[nodiscard]] bool HasStaticE( void ) const {
  return( get_NNodes() );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns true if there are static arcs (= flow variables if constructed)

 [[nodiscard]] bool HasStaticX( void ) const { return( get_NArcs() ); }

/*--------------------------------------------------------------------------*/
 /// given a pointer to a flow Variable, returns the index of the arc
 /** Given a pointer to a flow Variable (formally a Variable *, but
  * immediately static_cast-ed to a ColVariable * right inside), returns the
  * index of the corresponding arc. Throws exception if the pointer is not to
  * a [Col]Variable of the SingleFlowDCRBlock. */

 [[nodiscard]] Index p2i_x( const Variable * var ) const {
  auto i = p2i_x_s( var );
  if( ( i >= 0 ) && ( i < int( get_NArcs() ) ) )
   return( i );

  throw( std::invalid_argument( "invalid arc pointer" ) );
  return( 0 );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// given an arc, returns the pointer to the corresponding flow variable
 /** Given the index of an arc, returns the pointer to the corresponding flow
  * variable (a ColVariable *). This ASSUMES THE Variable ARE CONSTRUCTED IN
  * THE FIRST PLACE, SEGFAULTS ARE BOUND TO HAPPEN OTHERWISE. */

 [[nodiscard]] ColVariable * i2p_x( Index i ) const {
  if( i >= get_NArcs() )
   throw( std::invalid_argument( "invalid arc name" ) );

  return( const_cast< ColVariable * >( & x[ i ] ) );

 }

/*--------------------------------------------------------------------------*/

 [[nodiscard]] Index p2i_r( const Variable * var ) const {
  auto i = p2i_r_s( var );
  if( ( i >= 0 ) && ( i < int( get_NArcs() ) ) )
   return( i );

  throw( std::invalid_argument( "invalid arc pointer" ) );
  return( 0 );
  }

 [[nodiscard]] ColVariable * i2p_r( Index i ) const {
  if( i < get_NArcs() )
   return( const_cast< ColVariable * >( & r[ i ] ) );

  throw( std::invalid_argument( "invalid arc pointer" ) );
  return( 0 ); 
  }

/*--------------------------------------------------------------------------*/
 /// given a pointer to a UB Constraint, returns the index of the arc
 /** Given a pointer to a UB Constraint (formally a Constraint *, but
  * immediately static_cast-ed to a LB0Constraint * right inside), returns the
  * index of the corresponding arc. Throws exception if the pointer is not to
  * a [LB0]Constraint of the SingleFlowDCRBlock. */

 [[nodiscard]] Index p2i_ub( const Constraint * cns ) const {
  auto i = p2i_ub_s( cns );
  if( ( i >= 0 ) && ( i < int( get_NArcs() ) ) )
   return( i );

  throw( std::invalid_argument( "invalid ub constraint pointer" ) );
  return( 0 );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// given an arc, returns the pointer to the corresponding UB Constraint
 /** Given the index of an arc, returns the pointer to the corresponding UB
  * Constraint (a LB0Constraint *). This ASSUMES THE Constraint ARE
  * CONSTRUCTED IN THE FIRST PLACE, SEGFAULTS ARE BOUND TO HAPPEN OTHERWISE.
  */

 [[nodiscard]] LB0Constraint * i2p_ub( Index i ) const {
  if( i >= get_NArcs() )
   throw( std::invalid_argument( "invalid arc name" ) );

  return( const_cast< LB0Constraint * >( & UB[ i ] ) );
  }

/*--------------------------------------------------------------------------*/
 /// given a pointer to a flow Constraint, returns the index of the arc
 /** Given a pointer to a flow Constraint (formally a Constraint *, but
  * immediately static_cast-ed to a FRowConstraint * right inside), returns
  * the index of the corresponding arc. Throws exception if the pointer is
  * not to a [FRow]Constraint of the SingleFlowDCRBlock. */

 [[nodiscard]] Index p2i_e( const Constraint * cns ) const {
  auto i = p2i_e_s( cns );
  if( ( i >= 0 ) && ( i < int( get_NNodes() ) ) )
   return( i );

  throw( std::invalid_argument( "invalid flow constraint pointer" ) );
  return( 0 );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns true if the arc is closed; deleted arcs are not closed

 [[nodiscard]] bool is_closed( Index arc ) const {
  return( ( ! is_deleted( arc ) ) && i2p_x( arc )->is_fixed() );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns true if the arc is deleted

 [[nodiscard]] bool is_deleted( Index arc ) const {
  return( ( ! C.empty() ) && std::isnan( C[ arc ] ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// given a node, returns the pointer to the corresponding UB Constraint
 /** Given the index of n node, returns the pointer to the corresponding flow
  * Constraint (a FRowConstraint *). This ASSUMES THE Constraint ARE
  * CONSTRUCTED IN THE FIRST PLACE, SEGFAULTS ARE BOUND TO HAPPEN OTHERWISE.
  */

 [[nodiscard]] FRowConstraint * i2p_e( Index i ) const {
  if( i >= get_NNodes() )
   throw( std::invalid_argument( "invalid arc name" ) );

  return( const_cast< FRowConstraint * >( &E[ i ] ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the vector of starting nodes

 [[nodiscard]] c_Subset & get_SN( void ) const { return( SN ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the starting node of arc i (0 <= i < get_NArcs())

 [[nodiscard]] Index get_SN( Index i ) const { return( SN[ i ] ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the vector of ending nodes

 [[nodiscard]] c_Subset & get_EN( void ) const { return( EN ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the ending node of arc i (0 <= i < get_NArcs())

 [[nodiscard]] Index get_EN( Index i ) const { return( EN[ i ] ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the vector of arc costs
 /** Returns a const reference to the vector of arc costs of size
  * get_MaxNArcs(). Note that the cost of a deleted arc is NaN. */

 [[nodiscard]] c_Vec_double & get_C( void ) const { return( C ); }

 /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the cost of arc i (0 <= i < get_NArcs()), NaN if deleted

 [[nodiscard]] double get_C( c_Index i ) const { return( C[ i ] ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the vector of arc upper bounds
 /** Returns a const reference to the vector of arc upper bounds. Note that
  * the returned vector can either be of size get_MaxNArcs() or of size 0, in
  * which case all arc upper bounds are assumed to be +Inf. */

 [[nodiscard]] c_Vec_double & get_U( void ) const { return( U ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the upper bound of arc i (0 <= i < get_NArcs())

 [[nodiscard]] double get_U( Index i ) const {
  return( U.empty() ? Inf< double >() : U[ i ] );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the vector of node deficits
 /** Returns a const reference to the vector of node deficits. Note that the
  * returned vector can either be of size get_MaxNNodes() or of size 0, in
  * which case all node deficits are assumed to be 0. Also, note that the
  * position i (0 <= i < get_NNodes()) in this vector correspond to the node
  * whose name is i + 1 as returned from get_SN() and get_EN(). */

 [[nodiscard]] c_Vec_double & get_B( void ) const { return( B ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the upper deficit of node i (0 <= i < get_NNodes())
 /** Returns the deficit of node i. Note that "node names" here go from 0 to
  * get_NNodes() - 1, despite the fact that get_SN() and get_EN() report node
  * "names" between 1 and get_NNodes(). */
 
 [[nodiscard]] double get_B( Index i ) const {
  return( B.empty() ? 0 : B[ i ] );
  }
  
/** @} ---------------------------------------------------------------------*/
/*--------------------- Methods for checking the Block ---------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for checking the Block
 *  @{ */

 /// returns true if the current solution is (approximately) flow feasible
 /** Returns true if the solution encoded in the current value of the flow
  * (x) Variable of the SingleFlowDCRBlock is approximately feasible w.r.t. the flow
  * conservation constraints only. This clearly requires the Variable of the
  * SingleFlowDCRBlock to have been defined, i.e., that generate_abstract_variables()
  * has been called prior to this method. The parameter feps is the relative
  * accuracy defining "approximately". The parameter "useabstract" has the
  * same meaning as in is_feasible() and is_optimal(). */

 bool flow_feasible( double feps , bool useabstract = false );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns true if the current solution is (approximately) bound feasible
 /** Returns true if the solution encoded in the current value of the flow
  * (x) Variable of the SingleFlowDCRBlock is approximately feasible w.r.t. the bound
  * constraints only. This clearly requires the Variable of the SingleFlowDCRBlock to
  * have been defined, i.e., that generate_abstract_variables() has been
  * called prior to this method. The parameter feps is the relative accuracy
  * defining "approximately". The parameter "useabstract" has the same
  * meaning as in is_feasible() and is_optimal(). */

 bool bound_feasible( double feps , bool useabstract = false );

/*--------------------------------------------------------------------------*/
 /// returns true if the current solution is approximately feasible
 /** Returns true if the solution encoded in the current value of the flow
  * (x) Variable of the SingleFlowDCRBlock is approximately feasible. This clearly
  * requires the Variable of the SingleFlowDCRBlock to have been defined, i.e., that
  * generate_abstract_variables() has been called prior to this method.
  *
  * The parameter for deciding what "approximately feasible" exactly means is
  * a single double value, representing the *relative* tolerance for
  * satisfaction of both flow conservation constraint and flow upper/lower
  * bounds. This value is to be found as:
  *
  * - if fsbc is not nullptr and it is a SimpleConfiguration< double >, then
  *   it is fsbc->f_value;
  *
  * - otherwise, if f_BlockConfig is not nullptr,
  *   f_BlockConfig->f_is_feasible_Configuration is not nullptr and it
  *   is a SimpleConfiguration< double >, then it is
  *   f_BlockConfig->f_is_feasible_Configuration->f_value;
  *
  * - otherwise, it is 0. */
 
 bool is_feasible( bool useabstract = false , Configuration *fsbc = nullptr )
  override;

 bool is_feasible_flow( bool useabstract = false , Configuration *fsbc = nullptr );

/** @} ---------------------------------------------------------------------*/
/*------------------------- Methods for R3 Blocks --------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for R3 Blocks
 *  @{ */

 /// gets an R3 Block of SingleFlowDCRBlock currently only the copy one
 /** Gets an R3 Block of the SingleFlowDCRBlock. The list of currently supported R3
  * Block is:
  *
  * - r3bc == nullptr: the copy (an SingleFlowDCRBlock identical to this)
  */

 Block * get_R3_Block( Configuration *r3bc = nullptr ,
		       Block * base = nullptr , Block * father = nullptr )
  override;

/*--------------------------------------------------------------------------*/
 /// maps back the solution from a copy SingleFlowDCRBlock to the current one
 /** Maps back the solution from a copy SingleFlowDCRBlock to the current one. The
  * parameter r3bc is useless (has to be nullptr). The parameter solc decides
  * which part of the solution is mapped:
  *
  * - if solc != nullptr and it is a SimpleConfiguration< int >, then it
  *   depends on solc->f_value:
  *
  *   = 1 means "only map the primal solution"
  *
  *   = 2 means "only map the dual solution"
  *
  *   = everything else (e.g., 0) means "map everything";
  *
  * - if solc == nullptr, f_BlockConfig != nullptr,
  *   f_BlockConfig->f_solution_Configuration != nullptr and it
  *   is a SimpleConfiguration< int >, then it depends on its f_value as in
  *   the previous case;
  *
  * - otherwise, everything (both the primal and the dual solution) is
  *   mapped.
  *
  * The same format applies verbatim to the case of primal or dual unbounded
  * rays (negative-cost unbounded cycles and cuts, respectively), although
  * one would expect only one of these to be found (but both may
  * theoretically do).
  *
  * Note that R3B may not contain some or all of the required solution, if
  * the corresponding Variable/Constraint have not been constructed yet:
  * this throws an exception. */ 

 void map_back_solution( Block *R3B , Configuration *r3bc = nullptr ,
			 Configuration *solc = nullptr ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// maps the solution of the current SingleFlowDCRBlock to a copy SingleFlowDCRBlock
 /** Maps the solution of the current SingleFlowDCRBlock to a copy SingleFlowDCRBlock. The
  * parameter r3bc is useless (has to be nullptr). The parameter solc decides
  * which part of the solution is mapped:
  *
  * - if solc != nullptr and it is a SimpleConfiguration< int >, then it
  *   depends on solc->f_value:
  *
  *   = 1 means "only map the primal solution"
  *
  *   = 2 means "only map the dual solution"
  *
  *   = everything else (e.g., 0) means "map everything";
  *
  * - if solc == nullptr, f_BlockConfig != nullptr,
  *   f_BlockConfig->f_is_solution_Configuration != nullptr and it
  *   is a SimpleConfiguration< int >, then it depends on its f_value as in
  *   the previous case;
  *
  * - otherwise, everything (both the primal and the dual solution) is
  *   mapped.
  *
  * The same format applies verbatim to the case of primal or dual unbounded
  * rays (negative-cost unbounded cycles and cuts, respectively), although
  * one would expect only one of these to be found (but both may
  * theoretically do).
  *
  * Note that the current SingleFlowDCRBlock may not contain some or all of the
  * required solution, if the corresponding Variable/Constraint have not
  * been constructed yet: this throws an exception. */ 

 void map_forward_solution( Block *R3B , Configuration *r3bc = nullptr ,
			    Configuration *solc = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /** No specific Configuration is required, hence expected, for SingleFlowDCRBlock.
  *
  * IMPORTANT NOTE: map_forward_Modification() only maps "physical"
  * Modification. The point is that if any part of the "abstract
  * representation" of SingleFlowDCRBlock is changed, the corresponding "abstract"
  * Modification is intercepted in add_Modification() and a "physical"
  * Modification is also issued. Hence, for any change in SingleFlowDCRBlock there
  * will always be both Modification "in flight", and therefore there is
  * no need (and good reasons not) to map both.
  *
  * In particular, the method handles the following Modification:
  *
  * - GroupModification
  *
  * - SingleFlowDCRBlockRngdMod
  *
  * - SingleFlowDCRBlockSbstMod
  *
  * - NBModification
  *
  * Any other Modification is ignored (and false is returned).
  *
  *     IMPORTANT NOTE: SingleFlowDCRBlockRngdMod ALLOW TO ADD/DELETE ARCS IN THE
  *     PROBLEM, WHICH ALSO CHANGES THE "NAMES" OF EXISTING ARCS. SingleFlowDCRBlock
  *     IMPLEMENTS map_forward_Modification() IN A WAY THAT IS ONLY
  *     GUARANTEED TO BE CORRECT IF:
  *
  *     = EITHER THE SET OF ARCS IS NEVER CHANGED;
  *
  *     = OR THE Modification ARE MAPPED IMMEDIATELY AFTER THEY ARE ISSUED.
  *
  * This is because otherwise SingleFlowDCRBlock should have to understand whether the
  * set of arc "names" in the Modification is still correct and do something
  * in case it is not, which is too complex to do at the moment.
  *
  * Note that for GroupModification, true is returned only if all the
  * inner Modification of the GroupModification return true.
  *
  * Note that if the issueAMod param is eModBlck, then it is "downgraded" to
  * eNoBlck: the method directly does "physical" changes, hence there is no
  * reason for it to issue "abstract" Modification with concerns_Block() ==
  * true. */

 bool map_forward_Modification( Block *R3B , c_p_Mod mod ,
				Configuration *r3bc = nullptr ,
				ModParam issuePMod = eNoBlck ,
				ModParam issueAMod = eModBlck ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /** No specific Configuration is required, hence expected, for SingleFlowDCRBlock.
  *
  * The current implementation of map_back_Modification() actually uses
  * map_forward_Modification() in reverse, so see the comments to the latter
  * method. */

 bool map_back_Modification( Block *R3B , c_p_Mod mod ,
			     Configuration *r3bc = nullptr ,
			     ModParam issuePMod = eNoBlck ,
			     ModParam issueAMod = eModBlck ) override;

/** @} ---------------------------------------------------------------------*/
/*----------------------- Methods for handling Solution --------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for handling Solution
 *  @{ */

 /// returns a DCRSolution representing the current solution of this SingleFlowDCRBlock
 /** Returns a DCRSolution representing the current solution status of this
  * SingleFlowDCRBlock. 
  */ 

 Solution * get_Solution( Configuration *solc = nullptr ,
			  bool emptys = true ) override;

 /*--------------------------------------------------------------------------*/
 /// returns the objective value of the current solution

 double get_objective_value( void ) {
  if( ! ( AR & HasObj ) )  // the objective is not there
   return( Inf< RealObjective::OFValue >() );
  c.compute();
  return( c.value() );
  }

/*--------------------------------------------------------------------------*/
 /// gets a contiguous interval of the flow solution
 /** Method to get the flow solution; upon return, the current value of the
  * flow solution for the i-th arc in \p rng is written in *( FSol + i ).
  * Note that if the right extreme of the range is >= get_NArcs() it is
  * ignored. */

 void get_x( Vec_double_it FSol , Range rng = Range( 0 , Inf< Index >() ) )
  const;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// gets the flow solution for an arbitrary subset of arcs
 /** Method to get the flow solution; upon return, the current value of the
  * flow solution for arc nms[ i ] for all 0 <= i <  nms.size() is written
  * in *( FSol + i ). Note that
  *
  *     nms IS ASSUMED TO BE ORDERED BY INCREASING Index */

 void get_x( Vec_double_it FSol , c_Subset & nms ) const;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// gets the flow solution of the given arc

 double get_x( Index arc ) const {
  if( arc >= get_NArcs() )
   throw( std::invalid_argument( "invalid arc name" ) );

  return( x[ arc ].get_value() );

  }

/*--------------------------------------------------------------------------*/
 /// gets a contiguous interval of the reserve solution
 /** Method to get the flow solution; upon return, the current value of the
  * flow solution for the i-th arc in \p rng is written in *( FSol + i ).
  * Note that if the right extreme of the range is >= get_NArcs() it is
  * ignored. */

 void get_r( Vec_double_it FSol , Range rng = Range( 0 , Inf< Index >() ) )
  const;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// gets the reserve solution for an arbitrary subset of arcs
 /** Method to get the flow solution; upon return, the current value of the
  * flow solution for arc nms[ i ] for all 0 <= i <  nms.size() is written
  * in *( FSol + i ). Note that
  *
  *     nms IS ASSUMED TO BE ORDERED BY INCREASING Index */

 void get_r( Vec_double_it FSol , c_Subset & nms ) const;


/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// gets the reserve solution of the given arc

  double get_r( Index arc ) const {
  if( arc >= get_NArcs() )
   throw( std::invalid_argument( "invalid arc name" ) );

  return( r[ arc ].get_value() );
  }

/*--------------------------------------------------------------------------*/

  void set_r( c_Vec_double_it fstrt ,
    Range rng = Range( 0 , Inf< Index >() ) );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 void set_r( c_Vec_double_it fstrt , c_Subset sbst );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 void set_r( Index arc , double FSol ) {
  if( arc >= get_NArcs() )
  throw( std::invalid_argument( "invalid arc name" ) );

  r[ arc ].set_value( FSol );
}

/*--------------------------------------------------------------------------*/
 /// sets a contiguous interval of the flow solution
 /** Method to set the flow solution; the values found in the c_Vec_double
  * starting from fstrt are copied into the value of the flow variable
  * x[ i ] for i in rng, in the same order. */

 void set_x( c_Vec_double_it fstrt ,
	     Range rng = Range( 0 , Inf< Index >() ) );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// sets a generic subset of the flow solution
 /** Method to set the flow solution; the values found in the c_Vec_double
  * starting from fstrt are copied into the value of the flow variable
  * x[ i ] for all i in sbst (that must be ordered in increasing sense), in
  * the same order. */

 void set_x( c_Vec_double_it fstrt , c_Subset sbst );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// sets the flow solution of the given arc

 void set_x( Index arc , double FSol ) {
  if( arc >= get_NArcs() )
   throw( std::invalid_argument( "invalid arc name" ) );
  
  x[ arc ].set_value( FSol );
  }

/** @} ---------------------------------------------------------------------*/
/*-------------------- Methods for handling Modification -------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for handling Modification
 *  @{ */

 /// returns true if there is any Solver "listening to this SingleFlowDCRBlock"
 /** Returns true if there is any Solver "listening to this SingleFlowDCRBlock", or if
  * the SingleFlowDCRBlock has to "listen" anyway because the "abstract" representation
  * is constructed, and therefore "abstract" Modification have to be generated
  * anyway to keep the two representations in sync.
  *
  * No, this should not be needed. In fact, if the "abstract" representation
  * is modified with the default eModBlck value of issueMod, it is issued
  * irrespectively to the value of anyone_there(); see Observer::issue_mod().
  * If the value of issueMod is anything else the  "abstract" representation
  * has been modified already and there is no point in issuing the
  * Modification.
  * Note that that Observer::issue_mod() does not check if the "abstract"
  * representation has been constructed, but this is clearly not
  * necessary, as the Modification we are speaking of are issued while
  * changing the "abstract" representation, if that has not been
  * constructed then it cannot issue Modification

 bool anyone_there( void ) const override {
  return( AR ? true : Block::anyone_there() );
  }
 */
/*--------------------------------------------------------------------------*/
 /// adding a new Modification to the SingleFlowDCRBlock
 /** Method for handling Modification.
  *
  * The version of SingleFlowDCRBlock has to intercept any "abstract Modification" that
  * modifies the "abstract representation" of the SingleFlowDCRBlock, and "translate"
  * them into both changes of the actual data structures and corresponding
  * "physical Modification". These Modification are those for which
  * Modification::concerns_Block() is true. Note, however, that before sending
  * the Modification to the Solver and/or the father Block, the
  * concerns_Block() value is set to false. This is because once it is passed
  * through this method, the "abstract Modification" has "already done its
  * duty" of providing the information to the SingleFlowDCRBlock, and this must not be
  * repeated. In particular, this would be an issue if the Modification would
  * be [map_forward or map_back]-ed, because inside of this method a "physical
  * Modification" doing the same job is surely issued. That Modification would
  * also be [map_forward or map_back]-ed, together with the original "abstract
  * Modification" that would pass again through this method (in the other
  * SingleFlowDCRBlock), which would mean that the "physical Modification" would be
  * issued twice.
  *
  * The following "abstract Modification" are handled:
  *
  * - GroupModification, that are simply unpacked into the individual
  *   sub-[Group]Modification and dealt with individually;
  *
  * - C05FunctionModRngd and C05FunctionModSbst changing coefficients coming
  *   from the (LinearFunction into the FRow)Objective, but *not* from the
  *   (LinearFunction into the FRow)Constraint;
  *
  * - RowConstraintMod changing the RHS of the bound constraints and both
  *   sides at once of the flow conservation ones, but not any other
  *   combination; and note that the RHS of the bound constraints may not
  *   be changeable at all if they have not been constructed, in which
  *   case there cannot be any Modification to handle here;
  *
  * - VariableMod fixing and un-fixing a flow ColVariable; however, note
  *   that *fixing is only permitted if the value() of the ColVariable is
  *   zero*, because that corresponds to closing the arc, exception being
  *   thrown otherwise.
  *
  * Any other Modification reaching the SingleFlowDCRBlock will lead to exception
  * being thrown.
  *
  * Note: any "physical" Modification resulting from processing an "abstract"
  *       one will be sent to the same channel (chnl). */

 void add_Modification( sp_Mod mod , ChnlName chnl = 0 ) override;

/** @} ---------------------------------------------------------------------*/
/*--------------- METHODS FOR PRINTING & SAVING THE SingleFlowDCRBlock ---------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for printing & saving the SingleFlowDCRBlock
 *  @{ */

 /// print the SingleFlowDCRBlock on an ostream with the given verbosity
 /** Protected method to print information about the SingleFlowDCRBlock; with the
  * "complete" level ('C') it outputs the SingleFlowDCRBlock in DIMACS format. */

 void print( std::ostream & output , char vlvl = 0 ) const override;

/*--------------------------------------------------------------------------*/
 /// extends Block::serialize( netCDF::NcGroup )
 /** Extends Block::serialize( netCDF::NcGroup ) to the specific format of a
  * SingleFlowDCRBlock. See SingleFlowDCRBlock::deserialize( netCDF::NcGroup ) 
  * for details of the format of the created netCDF group. */

 void serialize( netCDF::NcGroup & group ) const override;

/** @} ---------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/
/** @name Changing the data of the DCR instance
 *
 * All the methods in this section have two parameters issueMod and issueAMod
 * which control if and how the, respectively, "physical Modification" and
 * "abstract Modification" corresponding to the change have to be issued, and
 * where (to which channel). The format of the parameters is that of
 * Observer::make_par(), except that the value eModBlck is ignored and
 * treated it as if it were eNoBlck [see Observer::issue_pmod()]. This is
 * because it makes no sense to issue an "abstract" Modification with
 * concerns_Block() == true, since the changes in the SingleFlowDCRBlock have surely
 * been done already, and this is just not possible for a "physical"
 * Modification.
 *
 * IMPORTANT NOTE: the current implementation of all these methods issues (at
 * most) *two separate* Modification, a "physical" and an "abstract" one. The
 * latter may be a GroupModification bunching together related abstract
 * Modification, but the two Modification are nonetheless separate. A
 * different approach could be to issue a single GroupModification with inside
 * both the "physical" and the "abstract" one (the latter possibly itself a
 * GroupModification). This may allow a more efficient handling of
 * Modification by ensuring that the two are always received together, but at
 * the cost of a more intricate code that is best avoided for now.
 *
 * Note: the methods accept the eDryRun value for the issueAMod parameter for
 * the "abstract" representation. This allows to re-use them within SingleFlowDCRBlock
 * itself when reacting to abstract Modification, where the  "abstract"
 * representation has been changed already.
 *  @{ */

 /// change the costs of a contiguous interval of arcs
 /** Method to change the costs of a subset of arcs with "contiguous names".
  * That is, *( NCost + i - strt ) becomes the new cost of the i-th arc in
  * \p rng. Note that if the right extreme of the range is >= get_NArcs() it 
  * is ignored.
  *
  * Note that if \p rng contains some closed arc, its cost is also changed.
  * While this has no immediate impact on the problem solved, if the arc is
  * re-opened then the cost set with this method when the arc was closed is
  * in effect.
  *
  * If more than one Modification is actually issued and issueAMod specifies
  * an open channel, then the channel is nested so that the three Modification
  * are grouped into a single GroupModification. Similarly, if instead
  * issueAMod specifies the default channel, then a new channel is opened to
  * group the multiple Modification and immediately closed when the last one
  * is issued. If, instead, the Objective is a "dense" LinearFunction, then
  * at most one LinearFunctionMod for modifying the coefficients is issued.
  * Of course this only applies if issueAMod specifies that abstract
  * Modification have to be issued *and* the abstract Objective has been
  * constructed.
  *
  * Also, if issueMod says so then a "physical" SingleFlowDCRBlockRngdMod is issued. */

 void chg_costs( c_Vec_double_it NCost , Range rng = INFRange ,
		 ModParam issueMod = eNoBlck , ModParam issueAMod = eNoBlck );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// change the costs of an arbitrary subset of arcs
 /** Method to change the costs of an arbitrary subset of arc. That is,
  * *( NCost + i ) becomes the new cost of arc nms[ i ] for all 0 <= i <
  * NCost.size(), (which means that nms.size() == NCost.size()). The
  * parameter ordered tells if the nms vector is ordered for increasing
  * index of the arc. As the && tells, nms is "consumed" by the method,
  * typically being shipped to an appropriate SingleFlowDCRBlockSbstMod object.
  *
  * See chg_costs( range ) for Modification issued (except that, of course,
  * the "physical" one is a SingleFlowDCRBlockSbstMod), and about changes in costs
  * of closed arcs. */

 void chg_costs( c_Vec_double_it NCost ,
		 Subset && nms , bool ordered = false ,
		 ModParam issueMod = eNoBlck , ModParam issueAMod = eNoBlck );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// changes the cost of the given arc
 /** Changes the cost of the given arc.
  *
  * Note that this can issue only one Modification of each type; the
  * "physical" one is a SingleFlowDCRBlockRngdMod with rng = [ arc ). */

 void chg_cost( double NCost , Index arc ,
		ModParam issueMod = eNoBlck , ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// change the capacities of a contiguous interval of arcs
 /** Method to change the capacities of a subset of arcs with "contiguous
  * names". That is, *( NCap + i - strt ) becomes the new capacity of the i-th
  * arc in \p rng. Note that if the right extreme of the range is
  * >= get_NArcs() it is ignored. Note that, according to the Configuration of
  * the static Constraint, the capacity of the arcs cannot be changed: trying
  * to do that will result in an exception being thrown.
  *
  * Note that if \p rng contains some closed arc, its capacity is also
  * changed. While this has no immediate impact on the problem solved, if the
  * arc is re-opened then the capacity set with this method when the arc was
  * closed is in effect.
  *
  * Note that changing the capacities can issue as many Modification as there
  * are arcs in the range, in particular OneVarConstraintMod with type
  * RowConstraintMod::eChgRHS. If more than one Modification is actually
  * issued and issueAMod specifies an open channel, then the channel is
  * nested so that all the Modification are grouped into a single
  * GroupModification. Similarly, if instead issueAMod specifies the default
  * channel, then a new channel is opened to group the multiple Modification
  * and immediately closed when the last one is issued. Of course this only
  * applies if issueAMod specifies that abstract Modification have to be
  * issued, *and* the abstract Constraint have been constructed.
  *
  * Note that, according to the Configuration of the static Constraint, the
  * capacity of the arcs cannot be changed: trying to do that will result in
  * an exception being thrown.
  *
  * Also, if issueMod says so then a "physical" SingleFlowDCRBlockRngdMod is issued. */

 void chg_ucaps( c_Vec_double_it NCap , Range rng = INFRange ,
		 ModParam issueMod = eNoBlck , ModParam issueAMod = eNoBlck );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// change the capacities of an arbitrary subset of arcs
 /** Method to change the capacities of an arbitrary subset of arc. That is,
  * *( NCap + i ) becomes the new capacity of arc nms[ i ] for all 0 <= i <
  * NCap.size() (which means that nms.size() == NCap.size()). The parameter
  * ordered tells if the nms vector is ordered for increasing index of the
  * arc. As the && tells, nms is "consumed" by the method, typically
  * being shipped to an appropriate SingleFlowDCRBlockSbstMod object.
  *
  * Note that, according to the Configuration of the static Constraint, the
  * capacity of the arcs cannot be changed: trying to do that will result in
  * an exception being thrown.
  *
  * See chg_ucaps( range ) for Modification issued (except that, of course,
  * the "physical" one is a SingleFlowDCRBlockSbstMod) and about changing capacities
  * of closed arcs. */

 void chg_ucaps( c_Vec_double_it NCap ,
		 Subset && nms , bool ordered = false ,
		 ModParam issueMod = eNoBlck , ModParam issueAMod = eNoBlck );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// change the capacity of the given arc
 /** Method to change the capacity of a given arc: NCap becomes the new
  * capacity of arc arc. Note that, according to the Configuration of the
  * static Constraint, the capacity of the arcs cannot be changed: trying to
  * do that will result in an exception being thrown.
  *
  * Note that this can issue only one Modification; the "physical" one is a
  * SingleFlowDCRBlockRngdMod with rng = [ arc ). */

 void chg_ucap( double NCap , Index arc ,
		ModParam issueMod = eNoBlck , ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// change the deficits of a contiguous interval of nodes
 /** Method to change the deficits of a subset of nodes with "contiguous
  * names". That is, *( NDfct + i - strt ) becomes the new deficit of the i-th
  * node in \p rng. Note that if the right extreme of the range is
  * >= get_NNodes() it is ignored. Note that "node names" here go from 0 to
  * get_NNodes() - 1, despite the fact that get_SN() and get_EN() report node
  * "names" between 1 and get_NNodes().
  *
  * Note that changing the capacities can issue as many Modification as there
  * are nodes in the range, in particular FRowConstraintMod with type
  * RowConstraintMod::eChgBTS. If more than one Modification is actually
  * issued and issueAMod specifies an open channel, then the channel is
  * nested so that all the Modification are grouped into a single
  * GroupModification. Similarly, if instead issueAMod specifies the default
  * channel, then a new channel is opened to group the multiple Modification
  * and immediately closed when the last one is issued. Of course this only
  * applies if issueAMod specifies that abstract Modification have to be
  * issued *and* the abstract Constraint have been constructed.
  *
  * Also, if issueMod says so then a "physical" SingleFlowDCRBlockRngdMod is issued. */

 void chg_dfcts( c_Vec_double_it NDfct , Range rng = INFRange ,
		 ModParam issueMod = eNoBlck , ModParam issueAMod = eNoBlck );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// change the deficits of an arbitrary subset of nodes
 /** Method to change the deficits of an arbitrary subset of nodes. That is,
  * *( NDfct + i ) becomes the new deficit of node nms[ i ] for all 0 <= i <
  * NDfct.size(), (which means that nms.size() == NDfct.size()). The
  * parameter ordered tells if the nms vector is ordered for increasing index
  * of the node. Note that "node names" here go from 0 to get_NNodes() - 1,
  * despite the fact that get_SN() and get_EN() report node "names" between
  * 1 and get_NNodes(). As the && tells, nms is "consumed" by the method,
  * typically being shipped to an appropriate SingleFlowDCRBlockSbstMod object.
  *
  * See chg_dfcts( range ) for Modification issued (except that, of course,
  * the "physical" one is a SingleFlowDCRBlockSbstMod). */

 void chg_dfcts( c_Vec_double_it NDfct ,
		 Subset && nms , bool ordered = false ,
		 ModParam issueMod = eNoBlck , ModParam issueAMod = eNoBlck );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// changes the deficit of the given node
 /** Method to change the deficit of a given node: NDfct becomes the new
  * deficit of node nde. Note that "node names" here go from 0 to
  * get_NNodes() - 1, despite the fact that get_SN() and get_EN() report node
  * "names" between 1 and get_NNodes().
  *
  * Note that this can issue only one Modification; the "physical" one is a
  * SingleFlowDCRBlockRngdMod with rng = [ arc ). */

 void chg_dfct( double NDfct , Index nde ,
		ModParam issueMod = eNoBlck , ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// closes a contiguous interval of arcs
 /** Method to close a subset of arcs with all "contiguous names" given in
  * \rng; note that if the right extreme of the range is >= get_NArcs() it is
  * ignored. The flow on the arcs is fixed to 0 but the arcs are not removed
  * from the problem, and their capacity and cost are not changed, so that
  * they can be easily re-opened later. When the problem is created, all arcs
  * are open. Closing an already closed arc does nothing.
  *
  * Note that closing multiple arcs can issue as many Modification as there
  * are arcs in the range, in particular VariableMod. If more than one
  * Modification is actually issued and issueAMod specifies an open channel,
  * then the channel is nested so that all the Modification are grouped into
  * a single GroupModification. Similarly, if instead issueAMod specifies the
  * default channel, then a new channel is opened to group the multiple
  * Modification and immediately closed when the last one is issued. Of course
  * this only applies if issueAMod specifies that abstract Modification have
  * to be issued.
  *
  * Also, if issueMod says so then a "physical" SingleFlowDCRBlockRngdMod is issued. */

 void close_arcs( Range rng = INFRange ,
		  ModParam issueMod = eNoBlck ,
		  ModParam issueAMod = eNoBlck );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// closes an arbitrary subset of arcs
 /** Method to close an arbitrary subset of arc, i.e., all those whose names
  * are found in the array nms. The flow on the arcs is fixed to 0 but the
  * arcs are not removed from the problem, and their capacity and cost are
  * not changed, so that they can be easily re-opened later. When the problem
  * is created, all arcs are open. Closing an already closed arc does
  * nothing.
  *
  * The parameter ordered tells if the nms vector is ordered for increasing 
  * index of the arc. As the && tells, nms is "consumed" by the method,
  * typically being shipped to an appropriate SingleFlowDCRBlockSbstMod object.
  *
  * See close_arcs( range ) for Modification issued (except that, of course,
  * the "physical" one is a SingleFlowDCRBlockSbstMod). */

 void close_arcs( Subset && nms , bool ordered = false ,
		  ModParam issueMod = eNoBlck ,
		  ModParam issueAMod = eNoBlck );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// closes the given arc
 /** Method to "close" the given arc: the flow on arc is fixed to 0. The arc
  * is not removed from the problem, and its capacity and cost are not
  * changed, so that it can be easily re-opened later. When the problem is
  * created, all arcs are open. Closing an already closed arc does nothing.
  *
  * Note that this can issue only one Modification; the "physical" one is a
  * SingleFlowDCRBlockRngdMod with rng = [ arc ). */

 void close_arc( Index arc , ModParam issueMod = eNoBlck ,
		             ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
/// re-opens a contiguous interval of arcs
 /** Method to "open" a subset of closed arcs with all "contiguous names"
  * given in \rng; note that if the right extreme of the range is >=
  * get_NArcs() it is ignored. Opening an already open arc (which is what all
  * arcs are when the problem is created) does nothing.
  *
  * Note that opening multiple arcs can issue as many Modification as there
  * are arcs in the range, in particular VariableMod. If more than one
  * Modification is actually issued and issueAMod specifies an open channel,
  * then the channel is nested so that all the Modification are grouped into
  * a single GroupModification. Similarly, if instead issueAMod specifies the
  * default channel, then a new channel is opened to group the multiple
  * Modification and immediately closed when the last one is issued. Of course
  * this only applies if issueAMod specifies that abstract Modification.
  *
  * Also, if issueMod says so then a "physical" SingleFlowDCRBlockRngdMod is issued. */

 void open_arcs( Range rng = INFRange ,
		 ModParam issueMod = eNoBlck ,
		 ModParam issueAMod = eNoBlck );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// re-opens an arbitrary subset of arcs
 /** Method to "open" an arbitrary subset of closed arc, i.e., all those
  * whose names are found in the array nms. Opening an already open arc
  * (which is what all arcs are when the problem is created) does nothing.
  *
  * The parameter ordered tells if the nms vector is ordered for increasing 
  * index of the arc. As the && tells, nms is "consumed" by the method,
  * typically being shipped to an appropriate SingleFlowDCRBlockSbstMod object.
  *
  * Note that closing multiple arcs can issue as many Modification as there
  * are arcs in the range, in particular VariableMod. If more than one
  * Modification is actually issued and issueAMod specifies an open channel,
  * then the channel is nested so that all the Modification are grouped into
  * a single GroupModification. Similarly, if instead issueAMod specifies the
  * default channel, then a new channel is opened to group the multiple
  * Modification and immediately closed when the last one is issued.
  * Of course this only applies if issueAMod specifies that abstract
  * Modification have to be issued *and* the abstract Constraint have been
  * constructed.
  *
  * Also, if issueMod says so then a "physical" SingleFlowDCRBlockRngdMod is issued. */

 void open_arcs( Subset && nms , bool ordered = false ,
		 ModParam issueMod = eNoBlck ,
		 ModParam issueAMod = eNoBlck );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// re-opens the given arc
 /** Method to "open" the given closed arc, i.e., allow the flow on arc to
  * vary. Opening an already open arc (which is what all arcs are when the
  * problem is created) does nothing.
  *
  * Note that this can issue only one Modification; the "physical" one is a
  * SingleFlowDCRBlockRngdMod with rng = [ arc ). */

 void open_arc( Index arc , ModParam issueMod = eNoBlck ,
		            ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/

 void chg_st( Index ns , Index nt , ModParam issueMod = eNoBlck , 
                    ModParam issueAMod = eNoBlck );

/** @} ---------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------------- PROTECTED FIELDS  ----------------------------*/
/*--------------------------------------------------------------------------*/

 Index NNodes;                   ///< the current number of nodes
 Index NArcs;                    ///< the current number of arcs
 Index MaxNNodes;                ///< the maximum number of nodes

 Subset SN;                      ///< vector of arc starting nodes
 Subset EN;                      ///< vector of arc ending nodes

 Vec_double C;                  ///< vector of arc costs
 Vec_double U;                  ///< vector of arc upper capacities
 Vec_double B;                  ///< vector of node deficits

 Vec_double NodeDelays;
 Vec_double LinkDelays;
 double FlowBursts;
 double FlowDeadlines; 
 double MTU; 
 double rho; 

 unsigned char AR;               ///< bit-wise coded: what abstract is there

 static constexpr unsigned char HasVar = 1;
 ///< first bit of AR == 1 if the Variable have been constructed
 static constexpr unsigned char HasObj = 2;
 ///< second bit of AR == 1 if the Objective has been constructed
 static constexpr unsigned char HasFlw = 4;
 ///< third bit of AR == 1 if the Flow Conservation have been constructed
 static constexpr unsigned char HasBnd = 8;
 ///< fourth bit of AR == 1 if the Bound have been constructed

 double f_cond_lower;            ///< conditional lower bound, can be -INF
 double f_cond_upper;            ///< conditional upper bound, can be +INF
 
 std::vector< ColVariable > x;     ///< the static flow variables
 std::vector< ColVariable > r;     ///< the static reserve variables
 std::vector< ColVariable > theta; ///< the static theta variables

 ColVariable r_min; ///< the static reserve_min variables
 ColVariable theta_min; ///< the static theta_min variables
 
 std::vector< FRowConstraint> E;   ///< the static flow conservation constrs.
 std::vector< LB0Constraint > UB;  ///< the static bound constraints on flow
 
 FRowConstraint DCR_cnst; /// the DCR constraint
 std::vector< FRowConstraint > Indicator_cnst_rmin; /// the static indicator constraints on reserve min
 std::vector< FRowConstraint > Indicator_cnst_r1; /// the first static indicator constraints on reserve
 std::vector< FRowConstraint > Indicator_cnst_r2; /// the second static indicator constraints on reserve

 FRowConstraint cone_min_cnst; /// the cone constraint
 std::vector< FRowConstraint > cone_cnst; /// the cone constraint
 
 std::list< FRowConstraint > PC_cuts;  /// the perspective dynamic cuts constraints
 std::list< FRowConstraint > PC_cuts_min;  /// the perspective dynamic cuts constraints

 FRealObjective c;               ///< the (linear) objective function

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/
/// register SingleFlowDCRBlock methods into the method factories
/** Although in general private methods should not be commented, this one is
 * because it does the registration of the following SingleFlowDCRBlock methods:
 *
 * - chg_costs() (both range and subset version)
 *
 * - chg_ucaps() (both range and subset version)
 *
 * - chg_dfcts() (both range and subset version)
 *
 * - close_arcs() (both range and subset version)
 *
 * - open_arcs() (both range and subset version)
 *
 * into the corresponding method factories. */

 static void static_initialization( void )
 {
  /*!!
 * Not all C++ compilers enjoy the template wizardry behind the three-args
 * version of register_method<> with the compact MS_*_*::args(), so we just
 * use the slightly less compact one with the explicit argument and be done
 * with it. !!*/
  // register_method< SingleFlowDCRBlock >( "SingleFlowDCRBlock::chg_costs", &SingleFlowDCRBlock::chg_costs,
  //                              MS_dbl_rngd::args() );
  //
  // register_method< SingleFlowDCRBlock >( "SingleFlowDCRBlock::chg_costs", &SingleFlowDCRBlock::chg_costs,
  //                              MS_dbl_sbst::args() );
  //
  // register_method< SingleFlowDCRBlock >( "SingleFlowDCRBlock::chg_ucaps", &SingleFlowDCRBlock::chg_ucaps,
  //                              MS_dbl_rngd::args() );
  //
  // register_method< SingleFlowDCRBlock >( "SingleFlowDCRBlock::chg_ucaps", &SingleFlowDCRBlock::chg_ucaps,
  //                              MS_dbl_sbst::args() );
  //
  // register_method< SingleFlowDCRBlock >( "SingleFlowDCRBlock::chg_dfcts", &SingleFlowDCRBlock::chg_dfcts,
  //                              MS_dbl_rngd::args() );
  //
  // register_method< SingleFlowDCRBlock >( "SingleFlowDCRBlock::chg_dfcts", &SingleFlowDCRBlock::chg_dfcts,
  //                              MS_dbl_sbst::args() );
  //
  // register_method< SingleFlowDCRBlock >( "SingleFlowDCRBlock::close_arcs",
  //                              &SingleFlowDCRBlock::close_arcs,
  //                              MS_rngd::args() );
  //
  // register_method< SingleFlowDCRBlock >( "SingleFlowDCRBlock::close_arcs",
  //                              &SingleFlowDCRBlock::close_arcs,
  //                              MS_sbst::args() );
  //
  // register_method< SingleFlowDCRBlock >( "SingleFlowDCRBlock::open_arcs", &SingleFlowDCRBlock::open_arcs,
  //                              MS_rngd::args() );
  //
  // register_method< SingleFlowDCRBlock >( "SingleFlowDCRBlock::open_arcs", &SingleFlowDCRBlock::open_arcs,
  //                              MS_sbst::args() );


  register_method< SingleFlowDCRBlock , MF_dbl_it , Range >( "SingleFlowDCRBlock::chg_costs" ,
						   & SingleFlowDCRBlock::chg_costs );

  register_method< SingleFlowDCRBlock , MF_dbl_it , Subset && , bool >(
   "SingleFlowDCRBlock::chg_costs" , & SingleFlowDCRBlock::chg_costs );

  register_method< SingleFlowDCRBlock , MF_dbl_it , Range >( "SingleFlowDCRBlock::chg_ucaps" ,
						   & SingleFlowDCRBlock::chg_ucaps );

  register_method< SingleFlowDCRBlock , MF_dbl_it , Subset &&, bool >(
   "SingleFlowDCRBlock::chg_ucaps" , & SingleFlowDCRBlock::chg_ucaps );
/*
  register_method< SingleFlowDCRBlock , MF_dbl_it , Range >( "SingleFlowDCRBlock::chg_dfcts" ,
						   & SingleFlowDCRBlock::chg_dfcts );

  register_method< SingleFlowDCRBlock , MF_dbl_it , Subset && , bool >(
   "SingleFlowDCRBlock::chg_dfcts" , & SingleFlowDCRBlock::chg_dfcts );
*/
  register_method< SingleFlowDCRBlock , Range >( "SingleFlowDCRBlock::close_arcs" ,
				       & SingleFlowDCRBlock::close_arcs );

  register_method< SingleFlowDCRBlock , Subset && , bool >( "SingleFlowDCRBlock::close_arcs" ,
						  & SingleFlowDCRBlock::close_arcs );

  register_method< SingleFlowDCRBlock , Range >( "SingleFlowDCRBlock::open_arcs" ,
				       & SingleFlowDCRBlock::open_arcs );

  register_method< SingleFlowDCRBlock , Subset && , bool >( "SingleFlowDCRBlock::open_arcs" ,
						  & SingleFlowDCRBlock::open_arcs );

  }  // end( static_initialization )

/*--------------------------------------------------------------------------*/

 int p2i_x_s( const Variable * var ) const {
  return( std::distance( x.data() ,
			 static_cast< const ColVariable * >( var ) ) );
  }

  int p2i_r_s( const Variable * var ) const {
  return( std::distance( r.data() ,
			 static_cast< const ColVariable * >( var ) ) );
  }

 int p2i_ub_s( const Constraint * cns ) const {
  return( std::distance( UB.data() ,
			 static_cast< const LB0Constraint * >( cns ) ) );
  }

 int p2i_e_s( const Constraint * cns ) const {
  return( std::distance( E.data() ,
			 static_cast< const FRowConstraint * >( cns ) ) );
  }

 LinearFunction * get_lfo( void ) {
  #ifdef NDEBUG
   return( static_cast< LinearFunction * >( c.get_function() ) );
  #else
   auto lfo = dynamic_cast< LinearFunction * >( c.get_function() );
   assert( lfo );
   return( lfo );
  #endif
  }

 LinearFunction * get_lfc( FRowConstraint * cnsti )
 {
  #ifdef NDEBUG
   return( static_cast< LinearFunction * >( cnsti->get_function() ) );
  #else
   auto lfc = dynamic_cast< LinearFunction * >( cnsti->get_function() );
   assert( lfc );
   return( lfc );
  #endif
  }

 void guts_of_destructor( void );

 void guts_of_add_Modification( p_Mod mod , ChnlName chnl );

 void compute_conditional_bounds( void );

/*--------------------------------------------------------------------------*/

#ifndef NDEBUG

 void CheckAbsVSPhys( void );
 
#endif

/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE FIELDS ------------------------------*/
/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;  // insert SingleFlowDCRBlock in the Block factory

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

 };  // end( class( SingleFlowDCRBlock ) )

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS SingleFlowDCRBlockMod -----------------------------*/
/*--------------------------------------------------------------------------*/
/// derived class from Modification for modifications to a SingleFlowDCRBlock
/** Derived class from Modification to describe modifications to a SingleFlowDCRBlock.
 *  This is actually "sort of abstract", since it does not say exactly what
 *  is changed, this being demanded to derived classes (which do this in
 *  different ways). Note that it is derived from Modification rather than,
 *  say, BlockMod (which has the same structure) because this is a class of
 *  "physical Modification". This means that a SingleFlowDCRBlockMod refers to changes
 *  in the "physical representation" of the SingleFlowDCRBlock; the corresponding
 *  changes in the "abstract representation" of the SingleFlowDCRBlock are dealt with
 *  by means of "abstract Modification", i.e., derived classes from
 *  AModification (as is BlockMod, which is why SingleFlowDCRBlockMod is not derived
 *  from BlockMod). */

class SingleFlowDCRBlockMod : public Modification
{
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/

 public:

/*---------------------------- PUBLIC TYPES --------------------------------*/
 /// public enum for the types of SingleFlowDCRBlockMod
 
 enum DCRB_mod_type {
  eChgCost = 0 ,   ///< change the arc costs
  eChgCaps     ,   ///< change the arc capacities
  eChgDfct     ,   ///< change the node deficits
  eOpenArc     ,   ///< open arcs
  eCloseArc    ,   ///< close arcs
  };

/*---------------------- CONSTRUCTOR & DESTRUCTOR --------------------------*/

 /// constructor: takes the SingleFlowDCRBlock and the type

 SingleFlowDCRBlockMod( SingleFlowDCRBlock * fblock , int type )
  : f_Block( fblock ) , f_type( type ) {}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

 virtual ~SingleFlowDCRBlockMod() = default;   ///< destructor, does nothing

/*-------------------- PUBLIC METHODS OF THE CLASS ------------------------*/

 /// returns the [DCR]Block to which the SingleFlowDCRBlockMod refers

 Block * get_Block( void ) const override  { return( f_Block ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// accessor to the type of modification

 int type( void ) const { return( f_type ); }

/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/

 protected:

/*-------------------------- PROTECTED METHODS -----------------------------*/
 /// print the SingleFlowDCRBlockMod

 void print( std::ostream &output ) const override {
  output << "SingleFlowDCRBlockMod[" << this << "]: ";
  switch( f_type ) {
   case( eChgCost ):  output << "change costs "; break;
   case( eChgCaps ):  output << "change capacities "; break;
   case( eChgDfct ):  output << "change deficits "; break;
   case( eOpenArc ):  output << "open arcs "; break;
   default:           output << "close arcs "; break;
   }
  }

/*--------------------- PROTECTED FIELDS OF THE CLASS ----------------------*/

 SingleFlowDCRBlock *f_Block;
               ///< pointer to the SingleFlowDCRBlock to which the SingleFlowDCRBlockMod refers

 int f_type;   ///< type of Modification

/*--------------------------------------------------------------------------*/

 };  // end( class( SingleFlowDCRBlockMod ) )

/*--------------------------------------------------------------------------*/
/*------------------------ CLASS SingleFlowDCRBlockRngdMod ---------------------------*/
/*--------------------------------------------------------------------------*/
/// derived from SingleFlowDCRBlockMod for "ranged" modifications
/** Derived class from SingleFlowDCRBlockMod to describe "ranged"
 * modifications to a SingleFlowDCRBlock, i.e., modifications that apply to an interval
 * of either arcs or nodes. */

class SingleFlowDCRBlockRngdMod : public SingleFlowDCRBlockMod
{
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/

 public:

/*---------------------- CONSTRUCTOR & DESTRUCTOR --------------------------*/

 /// constructor: takes the SingleFlowDCRBlock, the type, and the range

 SingleFlowDCRBlockRngdMod( SingleFlowDCRBlock * fblock , int type , Block::Range rng )
  : SingleFlowDCRBlockMod( fblock , type ) , f_rng( rng ) {}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

 virtual ~SingleFlowDCRBlockRngdMod() = default;   ///< destructor, does nothing

/*-------------------- PUBLIC METHODS OF THE CLASS ------------------------*/

 /// accessor to the range

 Block::c_Range & rng( void ) const { return( f_rng ); }
 
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/

 protected:

/*-------------------------- PROTECTED METHODS -----------------------------*/
 /// print the SingleFlowDCRBlockRngdMod

 void print( std::ostream &output ) const override {
  SingleFlowDCRBlockMod::print( output );
  output << "[ " << f_rng.first << ", " << f_rng.second << " )" << std::endl;
  }

/*--------------------- PROTECTED FIELDS OF THE CLASS ----------------------*/

 Block::Range f_rng;     ///< the range

/*--------------------------------------------------------------------------*/

 };  // end( class( SingleFlowDCRBlockRngdMod ) )

/*--------------------------------------------------------------------------*/
/*------------------------ CLASS SingleFlowDCRBlockSbstMod ---------------------------*/
/*--------------------------------------------------------------------------*/
/// derived from SingleFlowDCRBlockMod for "subset" modifications
/** Derived class from Modification to describe "subset" modifications to a
 *  SingleFlowDCRBlock, i.e., modifications that apply to an arbitrary subset of either
 * the arcs or the nodes. */

class SingleFlowDCRBlockSbstMod : public SingleFlowDCRBlockMod
{
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/

 public:


/*---------------------- CONSTRUCTOR & DESTRUCTOR --------------------------*/

 ///< constructor: takes the SingleFlowDCRBlock, the type, and the subset
 /**< Constructor: takes the SingleFlowDCRBlock, the type, and the subset. As the the
  * && tells, nms is "consumed" by the constructor and its resources become
  * property of the SingleFlowDCRBlockSbstMod object.
  *
  *   NOTE THAT nms IS REQUIRED TO BE ORDERED IN INCREASING SENSE
  *
  * although this is not checked by the class. */

 SingleFlowDCRBlockSbstMod( SingleFlowDCRBlock * fblock , int type , Block::Subset && nms )
  : SingleFlowDCRBlockMod( fblock , type ) , f_nms( std::move( nms ) ) {}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

 virtual ~SingleFlowDCRBlockSbstMod() = default;  ///< destructor, does nothing

/*-------------------- PUBLIC METHODS OF THE CLASS ------------------------*/

 /// accessor to the subset

 Block::c_Subset & nms( void ) const { return( f_nms ); }

/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/

 protected:

/*-------------------------- PROTECTED METHODS -----------------------------*/
 /// print the SingleFlowDCRBlockSbstMod

 void print( std::ostream &output ) const override {
  SingleFlowDCRBlockMod::print( output );
  output << "(# " << f_nms.size() << ")" << std::endl;
  }

/*--------------------- PROTECTED FIELDS OF THE CLASS ----------------------*/

 Block::Subset f_nms;   ///< the subset

/*--------------------------------------------------------------------------*/

 };  // end( class( SingleFlowDCRBlockSbstMod ) )

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS DCRSolution -----------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// a solution of a SingleFlowDCRBlock
/** The DCRSolution class, derived from Solution, represents a solution of a
 * SingleFlowDCRBlock
 */

class DCRSolution : public Solution {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*------------------------------- FRIENDS ----------------------------------*/

 friend SingleFlowDCRBlock;  ///< make SingleFlowDCRBlock friend

/*---------------- CONSTRUCTING AND DESTRUCTING DCRSolution ----------------*/

 explicit DCRSolution( void ) { }  /// constructor, it has nothing to do

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 void deserialize( const netCDF::NcGroup & group ) override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 ~DCRSolution() = default;  ///< destructor: it is virtual, and empty

/*------------- METHODS DESCRIBING THE BEHAVIOR OF A DCRSolution -----------*/

 void read( const Block * block ) override final;

 void write( Block * block ) override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// serialize a DCRSolution into a netCDF::NcGroup

 void serialize( netCDF::NcGroup & group ) const override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 DCRSolution * scale( double factor ) const override final;

 void sum( const Solution * solution , double multiplier ) override final;

 DCRSolution * clone( bool empty = false ) const override final;

/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/

 protected:

/*-------------------------- PROTECTED METHODS -----------------------------*/

 void print( std::ostream &output ) const override final {
  output << "DCRSolution [" << this << "]: " << v_r.size() 
    << " flows" << std::endl;
  }

/*---------------------- PRIVATE PART OF THE CLASS -------------------------*/

 private:

/*---------------------------- PRIVATE FIELDS ------------------------------*/

 SingleFlowDCRBlock::Vec_double v_x;   ///< the arc integer variables
 SingleFlowDCRBlock::Vec_double v_r;  ///< the arc flows

/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 };  // end( class( DCRSolution ) )

/** @} end( group( SingleFlowDCRBlock_CLASSES ) ) --------------------------*/
/*--------------------------------------------------------------------------*/

 };  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* SingleFlowDCRBlock.h included */

/*--------------------------------------------------------------------------*/
/*------------------- End File SingleFlowDCRBlock.h ------------------------*/
/*--------------------------------------------------------------------------*/
