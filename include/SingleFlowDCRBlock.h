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

 class MCFSolution;  // forward declaration of MCFSolution

/*--------------------------------------------------------------------------*/
/*----------------------- SingleFlowDCRBlock-RELATED TYPES ---------------------------*/
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
/*-------------------------- CLASS SingleFlowDCRBlock --------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// implementation of the Block concept for the (linear) Min-Cost Flow problem
/** The SingleFlowDCRBlock class implements the Block concept [see Block.h] for the
 * (linear) Min-Cost Flow (MCF) problem.
 *
 * The data of the problem consist of a (directed) graph G = ( N , A ) with
 * n = |N| nodes and m = |A| (directed) arcs. Each node `i' has a deficit
 * b[ i ], i.e., the amount of flow that is produced/consumed by the node:
 * source nodes (which produce flow) have negative deficits and sink nodes
 * (which consume flow) have positive deficits. Each arc `(i, j)' has an
 * upper capacity U[ i , j ] and a linear cost coefficient C[ i , j ]. Flow
 * variables X[ i , j ] represents the amount of flow to be sent on arc
 * (i, j). Parallel arcs, i.e., multiple copies of the same arc `(i, j)' are
 * in general allowed; it is expected that they have different costs (for
 * otherwise they can be merged into a unique arc), but this is not strictly
 * enforced. Multiple copies of some arc (i, j) can be seen as "total" flow
 * cost on that arc being a piecewise-linear convex function.
 *
 * The formulation of the problem is:
 * \f[
 *  \min \sum_{ (i, j) \in A } C[ i , j ] X[ i, j ]
 * \f]
 * \f[
 *  \sum_{ (j, i) \in A } X[ j , i ] -
 *  \sum_{ (i, j) \in A } X[ i , j ] = b[ i ] \quad i \in N      (1)
 * \f]
 * \f[
 *   0 \leq X[ i , j ] \leq U[ i , j ]   \quad (i, j) \in A      (2)
 * \f]
 * The n equations (1) are the flow conservation constraints and the 2m
 * inequalities (2) are the flow nonnegativity and capacity constraints.
 * At least one of the flow conservation constraints is redundant, as the
 * demands must be balanced (\f$\sum_{ i \in N } b[ i ] = 0\f$); indeed,
 * exactly n - ConnectedComponents( G ) flow conservation constraints are
 * redundant, as demands must be balanced in each connected component of G.
 *
 * The dual of the problem is:
 * \f[
 *  \max \sum_{ i \in N } Pi[ i ] b[ i ] -
 *       \sum_{ (i, j) \in A } W[ i , j ] U[ i , j ] -
 * \f]
 * \f[
 *  Pi[ j ] - Pi[ i ] - W[ i , j ] + Z[ i , j ] = C[ i , j ]
 *  \quad (i, j) \in A                                           (3)
 * \f]
 * \f[
 *  W[ i , j ] \geq 0 \quad (i, j) \in A                         (4)
 * \f]
 * \f[
 *  Z[ i , j ] \geq 0 \quad (i, j) \in A                         (5)
 * \f]
 *
 * Pi[] is said the vector of node potentials for the problem, W[] are bound
 * variables and Z[] are slack variables. Given Pi[], the quantities
 * \f[
 *  RC[ i , j ] =  C[ i , j ] + Q[ i , j ] * X[ i , j ] - Pi[ j ] + Pi[ i ]
 * \f]
 * are said the "reduced costs" of arcs, and are basically the dual variables
 * of the box constraints.
 *
 * A primal and dual feasible solution pair is optimal if and only if the
 * complementary slackness conditions
 * \f[
 *  RC[ i , j ] > 0 \Rightarrow X[ i , j ] = 0                   (6)
 * \f]
 * \f[
 *  RC[ i , j ] < 0 \Rightarrow X[ i , j ] = U[ i , j ]          (7)
 * \f]
 * are satisfied for all arcs (i, j) of A.
 *
 * The graph G is allowed to be "partly dynamic". The sets of nodes and arcs
 * that are input at the beginning are assumed not to be changed (save for
 * changing costs, capacities, and deficits, and for arcs to be closed or
 * opened). Then, new arcs and nodes, up to a set maximum, can be dynamically
 * added or deleted. This means that the graph can be fully static (if the
 * maximum number of dynamic arcs and nodes is set to zero) as well as fully
 * dynamic (if the initial graph is empty). Note that deleting the very last
 * arc decreases the number of arcs, while deleting one "in the middle" just
 * leaves "a hole": the arc "is not there" and any newly created arc can
 * "take its name", but the reported total number of arcs do not change.
 *
 * Note that changing costs, capacities and deficits is also allowed via the
 * abstract representation. Similarly, opening and closing arcs can be
 * performed by unfixing and fixing (respectively) the corresponding flow
 * variable. However, all other operations require "complex work" on the
 * abstract representation and therefore cannot be performed via that. One
 * could in principle allow it provided that all the Modification be grouped
 * in a GroupModification allowing to check that all the necessary operations
 * to, say, create and delete one arc have been properly done in the
 * abstract representation, but this is not implemented yet (and it's
 * doubtful it ever will). Thus, adding/removing Variable to flow
 * conservation constraints, or even changing their coefficients, via the
 * abstract representation is not allowed, as is (not) adding/removing
 * dynamic Constraint (be them flow conservation or bound ones). Similarly,
 * deleted arcs "in the middle" correspond to flow variables fixed to 0,
 * which cannot be unfixed via the abstract representation. In all these
 * cases, exceptions will be thrown. */

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
 *
 * - FNumber, the type of flow variables, arc capacities, and node deficits;
 *
 * - CNumber, the type of flow costs, node potentials, and arc reduced costs;
 *
 * - FONumber, the type of objective function value.
 *
 * By re-defining the types in this section, some (but not all) solution
 * algorithms may be able to work with the "smallest" choice of data type 
 * that is capable of properly representing the data of the instances to be
 * solved. This may be relevant due to an important property of MCF problems:
 * *if all arc capacities and node deficits are integer, then there exists an
 * integral optimal primal solution*, and *if all arc costs are integer,
 * then there exists an integral optimal dual solution*. Even more
 * importantly, *many solution algorithms will in fact produce an integral
 * primal/dual solution for free*, because *every primal/dual solution they
 * generate during the solution process is naturally integral*. Therefore,
 * one can use integer data types to represent everything connected with
 * flows and/or costs if the corresponding data is integer in all instances
 * one needs to solve. This directly translates in significant memory savings
 * and/or speed improvements.
 *
 * However, while using a SingleFlowDCRBlock as a part of some larger problem, it may
 * be difficult to fully exploit this property: even if some Solver can
 * exploit it, not all of them may be able to (one example are Interior-Point
 * approaches, which require both flow and cost variables to be continuous),
 * and maybe some other aspects of the overall solution algorithm will require
 * general double data anyway. One should actually have Block template over
 * all these types to be able to fully exploit this property, which may be a
 * future evolution but is not what this implementation does. The current
 * choice is to use the "worst case scenario" where FNumber == CNumber ==
 * OFNumber == double, although the data types are left there and it is
 * therefore in principle possible to change this. Note, however, that the
 * above integrality property only holds for *linear* MCF problems. Should
 * the class be extended, by even allowing arc costs to be convex quadratic
 * (the simplest possible nonlinear extension), then a single arc with a
 * nonzero quadratic cost coefficient implies that optimal flows and
 * potentials may be fractional even if all the data of the problem
 * (comprised quadratic cost coefficients) is integer. Hence, for such a
 * setting FNumber == CNumber == OFNumber == double is actually *mandatory*,
 * for any reasonable algorithm will typically misbehave otherwise.
 @{ */

/*--------------------------------------------------------------------------*/

 typedef double FNumber;                     ///< type of arc flow / deficit
 typedef const FNumber c_FNumber;            ///< a read-only FNumber

 typedef std::vector< FNumber > Vec_FNumber; ///< a vector of FNumber
 typedef const Vec_FNumber c_Vec_FNumber;    ///< a const vector of FNumber

 typedef Vec_FNumber::iterator Vec_FNumber_it;   ///< iterator in Vec_FNumber
 typedef Vec_FNumber::const_iterator c_Vec_FNumber_it;
                                           ///< const iterator in Vec_FNumber

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

 typedef double CNumber;                     ///< type of arc cost / potential
 typedef const CNumber c_CNumber;            ///< a read-only CNumber

 typedef std::vector< CNumber > Vec_CNumber;  ///< a vector of CNumber
 typedef const Vec_CNumber c_Vec_CNumber;     ///< a const vector of CNumber

 typedef Vec_CNumber::iterator Vec_CNumber_it;   ///< iterator in Vec_CNumber
 typedef Vec_CNumber::const_iterator c_Vec_CNumber_it;
                                           ///< const iterator in Vec_CNumber

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

 typedef double FONumber; 
 /**< type of the objective function: has to hold sums of products of
    FNumber(s) by CNumber(s) */

 typedef const FONumber c_FONumber;             ///< a read-only FONumber

 typedef std::vector< FONumber > Vec_FONumber;  ///< a vector of FONumber
 typedef const Vec_FONumber c_Vec_FONumber;     ///< a const vector of FONumber

/** @} ---------------------------------------------------------------------*/
/*------------------------------- FRIENDS ----------------------------------*/
/*--------------------------------------------------------------------------*/

 friend MCFSolution;  ///< make MCFSolution friend

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
  : Block( father ) , NNodes( 0 ) , NArcs( 0 ) , MaxNNodes( 0 ) ,
    NStaticNodes( 0 ) , NStaticArcs( 0 ) , AR( 0 ) ,
    f_cond_lower( - Inf< double >() ) , f_cond_upper( - Inf< double >() ) { }

/*--------------------------------------------------------------------------*/
 /// destructor of SingleFlowDCRBlock: deletes the abstract representation, if any

 virtual ~SingleFlowDCRBlock() { guts_of_destructor(); }

/** @} ---------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */

 /// loads the MCF instance from memory
 /** Loads the MCF instance from memory. The parameters are what you expect:
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
  * - dn   the current number (<= n, default 0) of dynamic nodes: all the
  *        nodes between 0 and n - dn - 1 are static, i.e., they cannot be
  *        deleted (and re-created), whereas all those from n - dn to n - 1
  *        are dynamic, i.e., they can deleted and later on re-created; if
  *        dn > n, then it is intended that dn == n, i.e., all nodes are
  *        dynamic (but still n is taken as the current number of nodes)
  *
  * - dm   the current number (<= m, default 0) of dynamic arcs: all the
  *        arcs between 0 and m - dm - 1 are static, i.e., they cannot be
  *        deleted (and re-created), whereas all those from m - dm to m - 1
  *        are dynamic, i.e., they can deleted and later on re-created; if
  *        dm > m, then it is intended that dm == m, i.e., all arcs are
  *        dynamic (but still m is taken as the current number of arcs)
  *
  * - mdn  the maximum number of dynamic nodes (default 0, if mdn < dn
  *        then the value is ignored and dn is used): data in the SingleFlowDCRBlock
  *        is allocated to accommodate for the fact that further mdn - dn
  *        nodes can later on be dynamically created (and deleted); if
  *        mdn < dn the parameter is ignored and n is taken as the maximum
  *        overall number of nodes
  *
  * - mdm  the maximum number of dynamic arcs (default 0, if mdm < dm
  *        then the value is ignored and dm is used): data in the SingleFlowDCRBlock
  *        is allocated to accommodate for the fact that further mdm - dm
  *        arcs can later on be dynamically created (and deleted); if
  *        mdm < dm the parameter is ignored and m is taken as the maximum
  *        overall number of arcs
  *
  * Like load( std::istream & ), if there is any Solver attached to this
  * SingleFlowDCRBlock then a NBModification (the "nuclear option") is issued. */

 void load( Index n , Index m , c_Subset & pEn , c_Subset & pSn ,
	    c_Vec_FNumber & pU = {} , c_Vec_CNumber & pC = {} ,
	    c_Vec_FNumber & pB = {} ,
	    Index dn = 0 , Index dm = 0 , Index mdn = 0 , Index mdm = 0 );

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
  * - the dimension "DynNNodes" containing the current number of dynamic
  *   nodes: all the nodes between 0 and "NNodes" - "DynNNodes" - 1 are
  *   static, i.e., they cannot be deleted (and re-created), whereas all
  *   those from "NNodes" + "DynNNodes" to "NNodes" - 1 are dynamic, i.e.,
  *   they can deleted and later on re-created; if "DynNNodes" > "NNodes",
  *   then it is intended that "DynNNodes" == "NNodes", i.e., all nodes are
  *   dynamic (but still "NNodes" is taken as the current number of nodes);
  *
  * - the dimension "DynNArcs" containing the current number of dynamic
  *   arcs: all the arcs between 0 and "NArcs" - "DynNArcs" - 1 are
  *   static, i.e., they cannot be deleted (and re-created), whereas all
  *   those from "NArcs" + "DynNArcs" to "NArcs" - 1 are dynamic, i.e.,
  *   they can deleted and later on re-created; if "DynNArcs" > "NArcs",
  *   then it is intended that "DynNArcs" == "NArcs", i.e., all arcs are
  *   dynamic (but still "NArcs" is taken as the current number of arcs);
  *
  * - the dimension "MaxDynNNodes" containing the maximum number of dynamic
  *   nodes: data in the SingleFlowDCRBlock is allocated to accommodate for the fact
  *   that further "MaxDynNNodes" - "DynNNodes" nodes can later on be
  *   dynamically created (and deleted); if "MaxDynNNodes" < "DynNNodes"
  *   the dimension is ignored and "NNodes" is taken as the maximum overall
  *   number of nodes;
  *
  * - the dimension "MaxDynNArcs" containing the maximum number of dynamic
  *   arcs: data in the SingleFlowDCRBlock is allocated to accommodate for the fact
  *   that further "MaxDynNArcs" - "DynNArcs" arcs can later on be
  *   dynamically created (and deleted); if "MaxDynNArcs" < "DynNArcs"
  *   the dimension is ignored and "NArcs" is taken as the maximum overall
  *   number of arcs.
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
 /// loads the MCF instance from file in DIMACS standard format
 /** Protected method for loading a SingleFlowDCRBlock out of a std::istream (which is
  * what operator>> is dispatched to. The std::istream is assumed to contain
  * the description of a MCF instance in DIMACS standard format, which is 
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
  * +Inf< FNumber >() by putting "INF", "Inf" or "inf" in the file (actually,
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

/*--------------------------------------------------------------------------*/
 /// generate the abstract variables of the MCF
 /** Method that generates the abstract Variable of the MCF. These are:
  *
  * - if ms = get_NStaticArcs() > 0, a std::vector< ColVariable > with
  *   exactly ms entries, the entry a = 0, ...,  ms - 1 corresponding to the
  *   flow on arc a, i.e., ( SN[ a ] , EN[ a ] );
  *
  * - if md = get_NArcs() - get_NStaticArcs() > 0, a std::list< ColVariable >
  *   with exactly md entries, the entry h = 0, ...,  md - 1 corresponding to
  *   the flow on arc a = ms + h, i.e., ( SN[ ms + h ] , EN[ ms + h ] ).
  *
  * Note that the dynamic Variable are actually created if get_MaxNArcs() >
  * get_NStaticArcs(), which may mean that the list can be empty when it is
  * created (if get_MaxNArcs() > get_NArcs() = get_NStaticArcs()); this is
  * done because new dynamic arcs can be created any time, and the list of
  * dynamic Variable is there ready for when this happens. */

 void generate_abstract_variables( Configuration *stvv = nullptr ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// generate the static constraint of the MCF
 /** Method that generates the abstract constraint of the MCF. These are:
  *
  * - if ns = get_NStaticNodes() > 0, a std::vector< FRowConstraint > with
  *   exactly ns entries, the entry i = 0, ..., ns - 1 being the flow
  *   conservation equation of the node i;
  *
  * - if nd = get_NNodes() - get_NStaticNodes() > 0, a
  *   std::list< FRowConstraint > with exactly nd entries, the entry
  *   h = 0 , ...,  nd - 1 being the flow conservation equation of the node
  *   i = ns + h;
  *
  * - if ms = get_NStaticArcs() > 0, a std::vector< LB0Constraint > with
  *   exactly ms entries, the entry a = 0, ..., ms - 1 being the bound
  *   constraints of the ColVariable x[ a ] corresponding to the flow on arc
  *   ( SN[ a ] , EN[ a ] );
  *
  * - if md = get_NArcs() - get_NStaticArcs() > 0, a
  *   std::list< LB0Constraint > with exactly md entries, the entry
  *   h = 0, ...,  md - 1  being the bound constraints of the ColVariable
  *   dx[ h ] corresponding to the flow on arc a = ms + h, i.e.,
  *   ( SN[ ms + h ] , EN[ ms + h ] ).
  *
  * Note that the dynamic flow conservation constraints are actually created
  * if get_MaxNNodes() > get_NStaticNodes(), which may mean that the list can
  * be empty when it is created (if get_MaxNNodes() > get_NNodes() =
  * get_NStaticNodes()); this is done because new dynamic nodes can be created
  * any time, and the list of dynamic Constraint is there ready for when this
  * happens.
  *
  * Regarding the bound constraints, these have fixed 0 LHS and a generic RHS,
  * which can be Inf< Fnumber >(). If *all* the RHS are +Infty, it is possible
  * to avoid creating the LB0Constraint entirely and just use the fact that
  * the ColVariable can be defined to be non-negative. The parameter stcc is
  * used to decide if this is done: if
  *
  * - all the RHS are +Infty;
  *
  * - either stcc is not nullptr and it is a SimpleConfiguration< int >;
  *
  * - or f_BlockConfig is not nullptr,
  *   f_BlockConfig->f_static_constraints_Configuration is not nullptr,
  *   and it is a SimpleConfiguration< int >;
  *
  * - the f_value of the SimpleConfiguration< int > is != 0
  *
  * the bound constraints are not implemented. Note that this *makes it
  * impossible to change any RHS*; the lower bound of 0 should not be changed
  * anyway, and there is no upper bound to be changed. This is done
  * consistently for the dynamic and static part, i.e., either both of them
  * are created, or none is. Actually, each part is only created if the
  * corresponding set of static/dynamic Variable exists; yet, note that the
  * dynamic Variable are actually created if get_MaxNArcs() >
  * get_NStaticArcs(), and thus the same is done for their dynamic bound
  * constraints (if they are created at all). This may mean that the list can
  * be empty when it is created (if get_MaxNArcs() > get_NArcs() =
  * get_NStaticArcs()); this is done because new dynamic arcs can be created
  * any time, and the list of dynamic bound Constraint is there ready for when
  * this happens.
  *
  * Note that it is only allowed to change the bounds of the LB0Constraint
  * (if any) and both bounds of the FRowConstraint (at the same time, and to
  * the same value): changing *any* the parts of any of the FRowConstraint,
  * such as the coefficients of the LinearFunction inside, is not allowed:
  * the SingleFlowDCRBlock will throw exception while processing the corresponding
  * "abstract" Modification. */
 
 void generate_abstract_constraints( Configuration *stcc = nullptr ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// generate the objective of the MCF
 /** Method that generates the objective of the MCF. Although this would seem
  * to be an exceedingly simple object, there is still a nontrivial decision
  * to be made about it, i.e., whether it is represented as a "sparse"
  * LinearFunction or a "dense" one. This is governed by objc: if
  *
  * - either objc is not nullptr and it is a SimpleConfiguration< double >;
  *
  * - or f_BlockConfig is not nullptr,
  *   f_BlockConfig->f_objective_Configuration is not nullptr,
  *   and it is a SimpleConfiguration< double >;
  *
  * then the f_value of the SimpleConfiguration< int > is taken as the "sparsity
  * parameter" (sprs) of the objective function; otherwise, sprs == 0. The
  * parameter is used as follows: if the number of nonzero arc cost
  * coefficients (at the time in which the method is called) is >=
  * sprs * get_NArcs(), then the LinearFunction in the objective is created
  * "dense": each flow Variable is "active" in it, even if it has a zero cost
  * coefficient. Otherwise, the LinearFunction in the objective is created
  * "sparse": only flow Variable with nonzero coefficient are "active" in it.
  * A "dense" Objective makes it much easier to change the cost coefficients
  * (see chg_cost[s]()), but it comes at the cost of more memory. Besides,
  * Solver using it and "trusting" the SingleFlowDCRBlock about how many nonzeroes are
  * there in the LinearFunction may be sorely disappointed, which may have
  * adverse effects on efficiency.
  *
  * Yet, a "dense" Objective is the default, as with sprs == 0 the Objective
  * is created "dense" even if all arc cost coefficients are zero.
  *
  * Note that the decision is taken at the moment in which this method is
  * called, and never changed later, even if the number of nonzeroes
  * changes dramatically. Also, note that if all cost coefficients are
  * "naturally" nonzero, then the Objective will be "dense" no matter what the
  * value of sprs is. Although this may seem obvious, this also means that the
  * Objective will remain "dense" even if later on many coefficients become
  * zero.
  *
  * IMPORTANT NOTE: ALLOWING SPARSE Objective MAKES IT INORDINATELY MORE
  * DIFFICULT TO REACT TO ABSTRACT Modification, WHILE ITS ACTUAL IMPACT ON
  * PERFORMANCES IS VERY DUBIOUS. THEREFORE, THE SUPPORT FOR IT IS ONLY
  * HALF-BAKED, AND WHATEVER THERE IS CURRENTLY COMMENTED OUT. DEVELOPMENT
  * OF THIS FEATURE WILL ONLY BE RESUMED IF CLEAR PROOF OF ITS WORTHINESS
  * IS ACHIEVED.
  *
  * The consequence is that, currently, THE ONLY Modification POSSIBLE TO THE
  * Objective ARE CHANGING THE COEFFICIENTS: DELETING Variable (AND,
  * THEREFORE, ADDING THEM) IS NOT ALLOWED, the SingleFlowDCRBlock will throw exception
  * while processing the corresponding "abstract" Modification. */

 void generate_objective( Configuration *objc = nullptr ) override;

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
 /// getting upper bounds on the value of the Objective
 /** An upper bound on the optimal value of the problem is computed as
  * \f[
  *  \sum_{ (i,j) \in A : c_{ij} > 0 } c_{ij} u_{ij}
  * \f]
  * If it is finite (which it may not be), this is a conditionally valid upper
  * bound but not a globally valid one because the problem may be empty (and
  * it being a minimization one this would mean that its optimal value is
  * + infinity).
  *
  * TODO: other bounds could be computed by looking at the total amount of
  *       flow to be moved
  *       \f[
  *         D = \sum_{ i \in N : b_i > 0 } b_i
  *       \f]
  *       and the worst possible cost of a simple path, like "max positive
  *       cost of an arc * ( n - 1 )". At least, D = 0 means that the problem
  *       is surely not empty, and thus the conditionally valid upper bound is
  *       also a globally valid upper bound. */

 [[nodiscard]] double get_valid_upper_bound( bool conditional = false )
  override {
  if( ! conditional )
   return( + Inf< double >() );
   
  if( std::isnan( f_cond_upper ) )
   compute_conditional_bounds();

  return( f_cond_upper );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// getting a global valid lower bound on the value of the Objective
 /** A lower bound on the optimal value of the problem is computed as
  * \f[
  *  \sum_{ (i,j) \in A : c_{ij} < 0 } c_{ij} u_{ij}
  * \f]
  * If it is finite (which it may not be), this is both a conditionally valid
  * lower bound abd a globally valid one, since clearly the problem then
  * cannot be unbounded below (although it can still be empty, but that's an
  * issue for upper bound, this being a minimization problem). */

 [[nodiscard]] double get_valid_lower_bound( bool conditional = false )
  override {
  if( std::isnan( f_cond_lower ) )
   compute_conditional_bounds();

  return( f_cond_lower );
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
 /// get the number of static nodes

 [[nodiscard]] Index get_NStaticNodes( void ) const {
  return( NStaticNodes );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the number of static arcs

 [[nodiscard]] Index get_NStaticArcs( void ) const {
  return( NStaticArcs );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns true if there are static nodes (= possibly flow constraints)

 [[nodiscard]] bool HasStaticE( void ) const {
  return( get_NStaticNodes() );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns true if there are dynamic nodes (= possibly flow constraints)

 [[nodiscard]] bool HasDynamicE( void ) const {
  return( get_NNodes() > get_NStaticNodes() );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns true if there may ever be dynamic nodes (= flow constraints)

 [[nodiscard]] bool MayHaveDynE( void ) const {
  return( get_MaxNNodes() > get_NStaticNodes() );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns true if there are static arcs (= flow variables if constructed)

 [[nodiscard]] bool HasStaticX( void ) const { return( get_NStaticArcs() ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns true if there are dynamic arcs (= flow variables if constructed)

 [[nodiscard]] bool HasDynamicX( void ) const {
  return( get_NArcs() > get_NStaticArcs() );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns true if there may ever be dynamic arcs  (= flow variables)

 [[nodiscard]] bool MayHaveDynX( void ) const {
  return( get_MaxNArcs() > get_NStaticArcs() );
  }

/*--------------------------------------------------------------------------*/
 /// given a pointer to a flow Variable, returns the index of the arc
 /** Given a pointer to a flow Variable (formally a Variable *, but
  * immediately static_cast-ed to a ColVariable * right inside), returns the
  * index of the corresponding arc. Throws exception if the pointer is not to
  * a [Col]Variable of the SingleFlowDCRBlock. */

 [[nodiscard]] Index p2i_x( const Variable * var ) const {
  auto i = p2i_x_s( var );
  if( ( i >= 0 ) && ( i < int( get_NStaticArcs() ) ) )
   return( i );

  i = get_NStaticArcs();
  for( auto dxi = dx.begin() ; dxi != dx.end() ; ++i , ++dxi )
   if( &(*dxi) == static_cast< const ColVariable * >( var ) )
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

  if( i < get_NStaticArcs() )
   return( const_cast< ColVariable * >( & x[ i ] ) );
  else
   if( i - get_NStaticArcs() < get_NArcs() - i )
    return( const_cast< ColVariable * >(
		   &( *std::next( dx.begin() , i - get_NStaticArcs() ) ) ) );
   else
    return( const_cast< ColVariable * >(
		           &( *std::prev( dx.end() , get_NArcs() - i ) ) ) );
  }

/*--------------------------------------------------------------------------*/
 /// given a pointer to a UB Constraint, returns the index of the arc
 /** Given a pointer to a UB Constraint (formally a Constraint *, but
  * immediately static_cast-ed to a LB0Constraint * right inside), returns the
  * index of the corresponding arc. Throws exception if the pointer is not to
  * a [LB0]Constraint of the SingleFlowDCRBlock. */

 [[nodiscard]] Index p2i_ub( const Constraint * cns ) const {
  auto i = p2i_ub_s( cns );
  if( ( i >= 0 ) && ( i < int( get_NStaticArcs() ) ) )
   return( i );

  i = get_NStaticArcs();
  for( auto dubi = dUB.begin() ; dubi != dUB.end() ; ++i , ++dubi )
   if( &(*dubi) == static_cast< const LB0Constraint * >( cns ) )
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

  if( i < get_NStaticArcs() )
   return( const_cast< LB0Constraint * >( & UB[ i ] ) );
  else
   if( i - get_NStaticArcs() < get_NArcs() - i )
    return( const_cast< LB0Constraint * >(
		  &( *std::next( dUB.begin() , i - get_NStaticArcs() ) ) ) );
   else
    return( const_cast< LB0Constraint * >(
		          &( *std::prev( dUB.end() , get_NArcs() - i ) ) ) );
  }

/*--------------------------------------------------------------------------*/
 /// given a pointer to a flow Constraint, returns the index of the arc
 /** Given a pointer to a flow Constraint (formally a Constraint *, but
  * immediately static_cast-ed to a FRowConstraint * right inside), returns
  * the index of the corresponding arc. Throws exception if the pointer is
  * not to a [FRow]Constraint of the SingleFlowDCRBlock. */

 [[nodiscard]] Index p2i_e( const Constraint * cns ) const {
  auto i = p2i_e_s( cns );
  if( ( i >= 0 ) && ( i < int( get_NStaticNodes() ) ) )
   return( i );

  i = get_NStaticNodes();
  for( auto dei = dE.begin() ; dei != dE.end() ; ++i , ++dei )
   if( &(*dei) == static_cast< const FRowConstraint * >( cns ) )
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

  if( i < get_NStaticNodes() )
   return( const_cast< FRowConstraint * >( &E[ i ] ) );
  else
   if( i - get_NStaticNodes() < get_NNodes() - i )
    return( const_cast< FRowConstraint * >(
		  &( *std::next( dE.begin() , i - get_NStaticNodes() ) ) ) );
   else
    return( const_cast< FRowConstraint * >(
		          &( *std::prev( dE.end() , get_NNodes() - i ) ) ) );
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

 [[nodiscard]] c_Vec_CNumber & get_C( void ) const { return( C ); }

 /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the cost of arc i (0 <= i < get_NArcs()), NaN if deleted

 [[nodiscard]] CNumber get_C( c_Index i ) const { return( C[ i ] ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the vector of arc upper bounds
 /** Returns a const reference to the vector of arc upper bounds. Note that
  * the returned vector can either be of size get_MaxNArcs() or of size 0, in
  * which case all arc upper bounds are assumed to be +Inf. */

 [[nodiscard]] c_Vec_FNumber & get_U( void ) const { return( U ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the upper bound of arc i (0 <= i < get_NArcs())

 [[nodiscard]] FNumber get_U( Index i ) const {
  return( U.empty() ? Inf< FNumber >() : U[ i ] );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the vector of node deficits
 /** Returns a const reference to the vector of node deficits. Note that the
  * returned vector can either be of size get_MaxNNodes() or of size 0, in
  * which case all node deficits are assumed to be 0. Also, note that the
  * position i (0 <= i < get_NNodes()) in this vector correspond to the node
  * whose name is i + 1 as returned from get_SN() and get_EN(). */

 [[nodiscard]] c_Vec_FNumber & get_B( void ) const { return( B ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the upper deficit of node i (0 <= i < get_NNodes())
 /** Returns the deficit of node i. Note that "node names" here go from 0 to
  * get_NNodes() - 1, despite the fact that get_SN() and get_EN() report node
  * "names" between 1 and get_NNodes(). */

 [[nodiscard]] FNumber get_B( Index i ) const {
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

 bool flow_feasible( FNumber feps , bool useabstract = false );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns true if the current solution is (approximately) bound feasible
 /** Returns true if the solution encoded in the current value of the flow
  * (x) Variable of the SingleFlowDCRBlock is approximately feasible w.r.t. the bound
  * constraints only. This clearly requires the Variable of the SingleFlowDCRBlock to
  * have been defined, i.e., that generate_abstract_variables() has been
  * called prior to this method. The parameter feps is the relative accuracy
  * defining "approximately". The parameter "useabstract" has the same
  * meaning as in is_feasible() and is_optimal(). */

 bool bound_feasible( FNumber feps , bool useabstract = false );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns true if the current solution is (approximately) dual feasible
 /** Returns true if the dual solution encoded in the current value of the
  * dual multipliers of both the flow conservation and bound constraints is
  * feasible. This clearly requires the Constraint of the SingleFlowDCRBlock to have
  * been generated by calling generate_abstract_constraints() prior to this
  * method. The parameter ceps is the relative accuracy defining
  * "approximately". The parameter "useabstract" has the same meaning as in
  * is_feasible() and is_optimal(). */

 bool dual_feasible( CNumber ceps , bool useabstract = false );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns true if complementary slackness are (approximately) satisfied
 /** Returns true if the (primal) solution encoded in the current value of
  * the flow (x) Variable of the SingleFlowDCRBlock and the dual solution encoded in
  * the current value of the dual multipliers of both the bound constraints
  * approximatively satisfy the Complementary Slackness Conditions. This
  * clearly requires both the Variable and the Constraint of the SingleFlowDCRBlock to
  * have been generated by calling generate_abstract_variables() and
  * generate_abstract_constraints() prior to this method. The parameters ceps
  * and feps are the relative accuracy defining "approximately" respectively
  * for "the reduced cost is zero" and "the flow is at the upper/lower bound".
  * The parameter "useabstract" has the same meaning as in is_feasible() and
  * is_optimal(). */

 bool complementary_slackness( CNumber ceps , FNumber feps ,
			       bool useabstract = false );

/*--------------------------------------------------------------------------*/
 /// returns true if the current solution is approximately feasible
 /** Returns true if the solution encoded in the current value of the flow
  * (x) Variable of the SingleFlowDCRBlock is approximately feasible. This clearly
  * requires the Variable of the SingleFlowDCRBlock to have been defined, i.e., that
  * generate_abstract_variables() has been called prior to this method.
  *
  * The parameter for deciding what "approximately feasible" exactly means is
  * a single FNumber value, representing the *relative* tolerance for
  * satisfaction of both flow conservation constraint and flow upper/lower
  * bounds. This value is to be found as:
  *
  * - if fsbc is not nullptr and it is a SimpleConfiguration< FNumber >, then
  *   it is fsbc->f_value;
  *
  * - otherwise, if f_BlockConfig is not nullptr,
  *   f_BlockConfig->f_is_feasible_Configuration is not nullptr and it
  *   is a SimpleConfiguration< FNumber >, then it is
  *   f_BlockConfig->f_is_feasible_Configuration->f_value;
  *
  * - otherwise, it is 0. */
 
 bool is_feasible( bool useabstract = false , Configuration *fsbc = nullptr )
  override;

/*--------------------------------------------------------------------------*/
 /// returns true if the current solution is (approximately) optimal
 /** Returns true if the solution encoded in the current value of the flow
  * (x) Variable of the SingleFlowDCRBlock is approximately optimal, which means that
  * it is approximately feasible, that the dual solution encoded in the
  * current value of the dual multipliers of both the flow conservation and
  * bound constraints is approximately feasible, and that the two
  * approximately satisfies the Complementary Slackness Conditions. This
  * clearly requires that both the Variable and the Constraint of the
  * SingleFlowDCRBlock to have been defined, i.e., that generate_abstract_variables()
  * and generate_abstract_constraints() have been called prior to this method.
  *
  * This requires two parameters for deciding what "approximately feasible"
  * means, one for the primal (feps) and one for the dual (ceps), like in
  * complementary_slackness(). These are found as follows:
  *
  * - if optc is not nullptr and it is a 
  *   SimpleConfiguration< std::pair< CNumber , FNumber > >, then
  *   ceps = optc->f_value.first and feps = optc->f_value.second;
  *
  * - if optc is not nullptr and it is a SimpleConfiguration< CNumber >, then
  *   ceps = optc->f_value, while feps is taken out of
  *   f_BlockConfig->f_is_feasible_Configuration as in is_feasible();
  *
  * - otherwise, if f_BlockConfig is not nullptr, then feps is taken
  *   out of f_BlockConfig->f_is_feasible_Configuration, while ceps
  *   is taken out of f_BlockConfig->f_is_optimal_Configuration
  *   assuming the latter is a SimpleConfiguration< CNumber >;
  *
  * - otherwise, ceps == feps == 0. */
 
 bool is_optimal( bool useabstract = false  , Configuration *optc = nullptr )
  override;

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

 /// returns a MCFSolution representing the current solution of this SingleFlowDCRBlock
 /** Returns a MCFSolution representing the current solution status of this
  * SingleFlowDCRBlock. What kind of solution is saved depends on the integer value ws,
  * obtained as follows:
  *
  * - if solc != nullptr and it is a SimpleConfiguration< int >, then
  *   ws == solc->f_value:
  *
  * - if solc == nullptr, f_BlockConfig != nullptr,
  *   f_BlockConfig->f_solution_Configuration != nullptr and it
  *   is a SimpleConfiguration< int >, ws is its f_value
  *
  * - otherwise ws is 0.
  *
  * The encoding of ws is:
  *
  *   = 1 means "only save the primal solution"
  *
  *   = 2 means "only save the dual solution"
  *
  *   = everything else (e.g., 0) means "save everything";
  *
  * The same format applies verbatim to the case of primal or dual unbounded
  * rays (negative-cost unbounded cycles and cuts, respectively), although
  * one would expect only one of these to be found (but both may
  * theoretically do).
  *
  * Note that SingleFlowDCRBlock may not contain some or all of the required solution,
  * if the corresponding Variable/Constraint have not been constructed yet:
  * this throws an exception, unless emptys = true, in which case the
  * MCFSolution object is only prepped for getting a solution, but it is not
  * really getting one now.
  *
  * Note that, although the method clearly returns a MCFSolution, formally
  * the return type is Solution *. This is because it is not possible to
  * forward declare MCFSolution as a derived class from Solution, nor to
  * define MCFSolution before SingleFlowDCRBlock because the former uses some type
  * information declared in the latter. */ 

 Solution * get_Solution( Configuration *solc = nullptr ,
			  bool emptys = true ) override;

 /*--------------------------------------------------------------------------*/
 /// returns the objective value of the current solution

 FONumber get_objective_value( void ) {
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

 void get_x( Vec_FNumber_it FSol , Range rng = Range( 0 , Inf< Index >() ) )
  const;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// gets the flow solution for an arbitrary subset of arcs
 /** Method to get the flow solution; upon return, the current value of the
  * flow solution for arc nms[ i ] for all 0 <= i <  nms.size() is written
  * in *( FSol + i ). Note that
  *
  *     nms IS ASSUMED TO BE ORDERED BY INCREASING Index */

 void get_x( Vec_FNumber_it FSol , c_Subset & nms ) const;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// gets the flow solution of the given arc

 FNumber get_x( Index arc ) const {
  if( arc >= get_NArcs() )
   throw( std::invalid_argument( "invalid arc name" ) );

  if( arc < get_NStaticArcs() )
   return( x[ arc ].get_value() );
  else
   return( std::next( dx.begin() , arc - get_NStaticArcs() )->get_value() );
  }

/*--------------------------------------------------------------------------*/
 /// gets a contiguous interval of the potential solution
 /** Method to get the potential solution; upon return, the current value of
  * the potential solution for the i-th node in \p rng is written into
  * *( PSol + i ). Note that if the right extreme of the range is >=
  * get_NNodes() it is ignored. Note that "node names" here go from 0 to
  * get_NNodes() - 1, despite the fact that get_SN() and get_EN() report node
  * "names" between 1 and get_NNodes(). */

 void get_pi( Vec_CNumber_it PSol , Range rng = Range( 0 , Inf< Index >() ) )
  const;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// gets the flow potential for an arbitrary subset of nodes
 /** Method to get the potential solution; upon return, the current value of
  * the potential solution for node nms[ i ] for all 0 <= i < nms.size() is
  * written into *( PSol + i ). Note that "node names" here go from 0 to
  * get_NNodes() - 1, despite the fact that get_SN() and get_EN() report node
  * "names" between 1 and get_NNodes(). Also, note that
  *
  *     nms IS ASSUMED TO BE ORDERED BY INCREASING Index */

 void get_pi( Vec_CNumber_it PSol , c_Subset & nms ) const;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// gets the potential solution of the given node
 /** Method to get the potential solution of the given node; note that "node
  * names" here go from 0 to get_NNodes() - 1, despite the fact that get_SN()
  * and get_EN() report node "names" between 1 and get_NNodes(). */

 CNumber get_pi( Index nde ) const {
  if( nde >= get_NNodes() )
   throw( std::invalid_argument( "invalid node name" ) );

  if( ! ( AR & HasFlw ) )
   throw( std::logic_error( "potentials unavailable if Constraint aren't" ) );

  if( nde < get_NStaticNodes() )
   return( E[ nde ].get_dual() );
  else
   return( std::next( dE.begin() , nde - get_NStaticNodes() )->get_dual() );
  }

/*--------------------------------------------------------------------------*/
 /// gets a contiguous interval of the reduced costs
 /** Method to get the reduced costs; upon return, the current value of the
  * reduced cost for the i-th arc in \p rng is written into *( RC + i ). Note
  * that if the right extreme of the range is >= get_NArcs() it is ignored. */

 void get_rc( Vec_CNumber_it RC , Range rng = Range( 0 , Inf< Index >() ) )
  const;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// gets the reduced costs for an arbitrary subset of arcs
 /** Method to get the reduced costs; upon return, the current value of the
  * reduced costs for arc nms[ i ] for all 0 <= i < nms.size() is written
  * into *( RC + i ). Note that
  *
  *     nms IS ASSUMED TO BE ORDERED BY INCREASING Index */

 void get_rc( Vec_CNumber_it RC , c_Subset & nms ) const;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// gets the reduced costs of the given arc

 CNumber get_rc( Index arc ) const {
  if( E.empty() && dE.empty() )
   throw( std::logic_error( "reduced costs unavailable if Constraint aren't"
			    ) );
  if( arc >= get_NArcs() )
   throw( std::invalid_argument( "invalid arc name" ) );

  if( UB.empty() && dUB.empty() )
   return( get_C( arc ) + get_pi( SN[ arc ] - 1 ) - get_pi( EN[ arc ] - 1 ) );
  else
   if( arc < get_NStaticArcs() )
    return( UB[ arc ].get_dual() );
   else
    return( std::next( dUB.begin() , arc - get_NStaticArcs() )->get_dual() );
  }

/*--------------------------------------------------------------------------*/
 /// sets a contiguous interval of the flow solution
 /** Method to set the flow solution; the values found in the c_Vec_FNumber
  * starting from fstrt are copied into the value of the flow variable
  * x[ i ] for i in rng, in the same order. */

 void set_x( c_Vec_FNumber_it fstrt ,
	     Range rng = Range( 0 , Inf< Index >() ) );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// sets a generic subset of the flow solution
 /** Method to set the flow solution; the values found in the c_Vec_FNumber
  * starting from fstrt are copied into the value of the flow variable
  * x[ i ] for all i in sbst (that must be ordered in increasing sense), in
  * the same order. */

 void set_x( c_Vec_FNumber_it fstrt , c_Subset sbst );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// sets the flow solution of the given arc

 void set_x( Index arc , FNumber FSol ) {
  if( arc >= get_NArcs() )
   throw( std::invalid_argument( "invalid arc name" ) );

  if( arc < get_NStaticArcs() )
   x[ arc ].set_value( FSol );
  else
   std::next( dx.begin() , arc - get_NStaticArcs() )->set_value( FSol );
  }

/*--------------------------------------------------------------------------*/
 /// sets a contiguous interval of the potential solution
 /** Method to set the potential solution; the values found in the
  * c_Vec_CNumber starting from pstrt are copied into the potential of node
  * (dual multiplier of the flow balance constraint) i for i in rng, in the
  * same order. */

 void set_pi( c_Vec_CNumber_it pstrt ,
	      Range rng = Range( 0 , Inf< Index >() ) );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// sets a generic subset of the potential solution
 /** Method to set the potential solution; the values found in the
  * c_Vec_FNumber starting from pstrt are copied into the potential of node
  * (dual multiplier of the flow balance constraint) i for all i in sbst
  * (that must be ordered in increasing sense), in the same order. */

 void set_pi( c_Vec_FNumber_it pstrt , c_Subset sbst );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// sets the potential solution of the given node

 void set_pi( CNumber PSol , Index nde ) {
  if( ! ( AR & HasFlw ) )  // nowhere to put the value
   return;                 // cowardly (and silently) return

  if( nde >= get_NNodes() )
   throw( std::invalid_argument( "invalid node name" ) );

  if( nde < get_NStaticNodes() )
   E[ nde ].set_dual( PSol );
  else
   std::next( dE.begin() , nde - get_NStaticNodes() )->set_dual( PSol );
  }

/*--------------------------------------------------------------------------*/
 /// sets a contiguous interval of the reduced costs
 /** Method to set the reduced costs solution; the values found in the
  * c_Vec_CNumber starting from rcstrt are copied into the reduced cost of
  * arc (dual value of the bound constraint) i for i in rng, in the same
  * order. */

 void set_rc( c_Vec_CNumber_it rcstrt ,
	      Range rng = Range( 0 , Inf< Index >() ) );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// sets a generic subset of the reduced costs
 /** Method to set the reduced costs solution; the values found in the
  * c_Vec_FNumber starting from rcstrt are copied into the reduced cost of
  * arc (dual value of the bound constraint) i for all i in sbst (that must
  * be ordered in increasing sense), in the same order. */

 void set_rc( c_Vec_FNumber_it rcstrt , c_Subset sbst );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// sets the reduced cost of the given arc

 void set_rc( CNumber RC , Index arc ) {
 if( ! ( AR & HasBnd ) )  // nowhere to put the value in
  return;                 // cowardly (and silently) return

 if( arc >= get_NArcs() )
  throw( std::invalid_argument( "invalid arc name" ) );

 if( arc < get_NStaticArcs() )
  UB[ arc ].set_dual( RC );
 else
  std::next( dUB.begin() , arc - get_NStaticArcs() )->set_dual( RC );

 }  // end( SingleFlowDCRBlock::set_rc( one ) )

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
  * SingleFlowDCRBlock. See SingleFlowDCRBlock::deserialize( netCDF::NcGroup ) for details of the
  * format of the created netCDF group. */

 void serialize( netCDF::NcGroup & group ) const override;

/** @} ---------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/
/** @name Changing the data of the MCF instance
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

 void chg_costs( c_Vec_CNumber_it NCost , Range rng = INFRange ,
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

 void chg_costs( c_Vec_CNumber_it NCost ,
		 Subset && nms , bool ordered = false ,
		 ModParam issueMod = eNoBlck , ModParam issueAMod = eNoBlck );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// changes the cost of the given arc
 /** Changes the cost of the given arc.
  *
  * Note that this can issue only one Modification of each type; the
  * "physical" one is a SingleFlowDCRBlockRngdMod with rng = [ arc ). */

 void chg_cost( CNumber NCost , Index arc ,
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

 void chg_ucaps( c_Vec_FNumber_it NCap , Range rng = INFRange ,
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

 void chg_ucaps( c_Vec_FNumber_it NCap ,
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

 void chg_ucap( FNumber NCap , Index arc ,
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

 void chg_dfcts( c_Vec_FNumber_it NDfct , Range rng = INFRange ,
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

 void chg_dfcts( c_Vec_FNumber_it NDfct ,
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

 void chg_dfct( FNumber NDfct , Index nde ,
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
 /// add a new arc
 /** Method to add a new arc, providing its starting and ending nodes, cost
  * and capacity.
  *
  * The method returns the "name" that the new arc has received, which is the
  * index by which the arc has to be addressed in all methods (such as
  * chg_[cost/ucap](), [open/close]_arc(), get_[x/rc]()). The chosen name
  * depends on whether or not there are "deleted" arcs (see remove_arc())
  * with "name" < get_NArcs() - 1. If there is any such arc, then the "name"
  * of the new arc will be the smallest index among these; in this case, the
  * value returned by get_NArcs() does *not* change. Otherwise, if get_NArcs()
  * < get_MaxNArcs() then the new arc gets "name" get_NArcs(), which is
  * returned by the method, and the value returned by get_NArcs() increases by
  * one. Otherwise the arc is not actually added, and the method returns
  * Inf< FNumber >().
  *
  * Successfully adding a new arc causes the issuing of several Modification,
  * unless the issueMod and issueAMod parameters prevent this to happen:
  *
  * - a "physical" SingleFlowDCRBlockRngdMod with type eAddArc;
  *
  * - an "abstract" GroupModification containing up to:
  *
  *   = two C05FunctionModVars corresponding to having added the new Variable
  *     to the two flow conservation constraints of its starting and ending
  *     node;
  *
  *   = one OneVarConstraintMod with type RowConstraintMod::eChgRHS for
  *      modifying the flow bound;
  *
  *   = if the "name" of the arc is == get_NArcs() (before the call):
  *
  *     * a BlockModAdd< ColVariable > corresponding to the addition of a
  *       new dynamic Variable (the flow Variable of the arc);
  *
  *     * possibly, a BlockModAdd< LB0Constraint > corresponding to the
  *       addition of a new dynamic Constraint (the bound Constraint of the
  *       arc, if it is defined);
  *
  *     * one more C05FunctionModVars corresponding to having added the new
  *       Variable to the objective.
  *
  *   = if, instead, the "name" of the arc is < get_NArcs() (before the call):
  *
  *     * one C05FunctionModLin for modifying the cost coefficients;
  *
  *     * one VariableMod making the flow variable "free";
  *       
  *  Of course, all the "abstract" Modification are only issued if the
  *  corresponding part of the "abstract" representation is constructed. */
 
 Index add_arc( Index sn , Index en , CNumber cst = 0 ,
		FNumber cap = Inf< FNumber >() ,
		ModParam issueMod = eNoBlck , ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// removes an existing arc
 /** Method to remove the arc which given name. It must be
  * get_NStaticArcs() <= arc < get_NArcs(), otherwise exception is thrown.
  *
  * The operation is performed differently in the case where arc ==
  * get_NArcs() - 1, i.e., the very last arc is eliminated, or
  * arc < get_NArcs() - 1.
  *
  * In the latter case, the elimination of the arc is "virtual", in the
  * sense that the flow variable is kept, together with all corresponding
  * parts of the "abstract" representation (if constructed). Only, the
  * value of the flow Variable is changed to 0 and the Variable is fixed, as
  * when the arc is closed. Furthermore, its starting and ending nodes (as
  * returned by get_SN() and get_EN()) are set to Inf< Index >(). This means
  * that the value returned by get_NArcs() does *not* change.
  *
  * In the former case, the elimination of the arc is "physical": not only
  * of that arc, but also of and all the "deleted" arcs with smaller name
  * up until the first non-deleted arc (or get_NStaticArcs()). The value
  * of get_NArcs() changes accordingly, and all the corresponding parts of
  * the "abstract" representation (Variable, bound constraints, coefficients
  * in the objective function and the constraints) are removed.
  *
  * Removing an existing arc causes the issuing of several Modification,
  * unless the issueMod and issueAMod parameters prevent this to happen:
  *
  * - a "physical" SingleFlowDCRBlockRngdMod with type eRmvArc;
  *
  * - an "abstract" GroupModification containing up to:
  *
  *   = for each of the removed arcs, two C05FunctionModVars corresponding
  *     to having removed the existing Variable from the two flow
  *     conservation constraints of its starting and ending node;
  *
  *   = if the elimination is "virtual", one VariableMod corresponding to
  *     fixing the flow variable;
  *
  *   = if the elimination is "physical":
  *
  *     * a BlockModRmv< ColVariable > corresponding to the removal of the
  *       existing dynamic Variable (the flow Variable of the arcs);
  *
  *     * possibly, a BlockModRmv< LB0Constraint > corresponding to the
  *       removal of the existing dynamic Constraint (the bound Constraint of
  *       the arcs, if they are defined);
  *
  *     * for each of the removed arcs, one more C05FunctionModVars
  *       corresponding to having removed the existing Variable from the
  *       objective.
  *
  *  Of course, all the "abstract" Modification are only issued if the
  *  corresponding part of the "abstract" representation is constructed. */

 void remove_arc( Index arc , ModParam issueMod = eNoBlck ,
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
 Index NStaticNodes;             ///< the number of static nodes
 Index NStaticArcs;              ///< the number of static arcs

 Subset SN;                      ///< vector of arc starting nodes
 Subset EN;                      ///< vector of arc ending nodes

 Vec_CNumber C;                  ///< vector of arc costs
 Vec_FNumber U;                  ///< vector of arc upper capacities
 Vec_FNumber B;                  ///< vector of node deficits

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
 std::vector< FRowConstraint> E;   ///< the static flow conservation constrs.
 std::vector< LB0Constraint > UB;  ///< the static bound constraints
 
 std::list< ColVariable > dx;      ///< the dynamic flow variables
 std::list< FRowConstraint > dE;   ///< the dynamic flow conservation constrs.
 std::list< LB0Constraint > dUB;   ///< the dynamic bound constraints

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

  register_method< SingleFlowDCRBlock , MF_dbl_it , Range >( "SingleFlowDCRBlock::chg_dfcts" ,
						   & SingleFlowDCRBlock::chg_dfcts );

  register_method< SingleFlowDCRBlock , MF_dbl_it , Subset && , bool >(
   "SingleFlowDCRBlock::chg_dfcts" , & SingleFlowDCRBlock::chg_dfcts );

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
 
 enum MCFB_mod_type {
  eChgCost = 0 ,   ///< change the arc costs
  eChgCaps     ,   ///< change the arc capacities
  eChgDfct     ,   ///< change the node deficits
  eOpenArc     ,   ///< open arcs
  eCloseArc    ,   ///< close arcs
  eAddArc      ,   ///< add arcs
  eRmvArc          ///< remove arcs
  };

/*---------------------- CONSTRUCTOR & DESTRUCTOR --------------------------*/

 /// constructor: takes the SingleFlowDCRBlock and the type

 SingleFlowDCRBlockMod( SingleFlowDCRBlock * fblock , int type )
  : f_Block( fblock ) , f_type( type ) {}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

 virtual ~SingleFlowDCRBlockMod() = default;   ///< destructor, does nothing

/*-------------------- PUBLIC METHODS OF THE CLASS ------------------------*/

 /// returns the [MCF]Block to which the SingleFlowDCRBlockMod refers

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
   case( eCloseArc ): output << "close arcs "; break;
   case( eAddArc ):   output << "add arcs "; break;
   default:           output << "remove arcs ";
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
/*-------------------------- CLASS MCFSolution -----------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// a solution of a SingleFlowDCRBlock
/** The MCFSolution class, derived from Solution, represents a solution of a
 * SingleFlowDCRBlock, i.e.:
 *
 * - an m-vector of FNumber for the arc flow values;
 *
 * - an n-vector of CNumber for the node potentials;
 *
 * where m is the number of arcs and n is the number of nodes in the graph.
 * This means that
 *
 *       THE REDUCED COSTS ARE NOT EXPLICITLY SAVED
 *
 * This is OK for feasible dual solutions, as the dual variables of the bound
 * constraints (a.k.a. Reduced Costs) can be cheaply computed out of the
 * potentials. This may not be appropriate in all cases, as one may want to
 * deal with unfeasible dual solutions; if this will ever be the case, the
 * MCFSolution class will have to be changed accordingly.
 *
 * Note that the vectors are (in principle, both) optional: a MCFSolution
 * may have them empty, according to how it is created by a call to
 * SingleFlowDCRBlock::get_Solution(). If a vector is empty it is never read() or
 * write()-n from/to the SingleFlowDCRBlock. There is no support for changing this
 * during the life of the MCFSolution.
 *
 * It is useful to remark that some special cases of MCF would actually have
 * "special" solutions ("less general" ones in the parlance of Solution). In
 * particular:
 *
 * - if all capacities are Inf< FNumber >() and there is only one source or sink
 *   node, then the MCF problem is in fact a Shortest Path (sub-)Tree one, and
 *   its solutions can be represented by means of a predecessor function;
 *
 * - if all (finite) capacities and node deficits are integer, then there
 *   always exist optimal flow solutions of MCF that are integer;
 *
 * - if all arc costs are integer, then there always exist optimal potential
 *   solutions of MCF that are integer.
 *
 * Thus, SingleFlowDCRBlock would have scope for different kinds of Solution objects.
 * The currently implemented one is the "most general" one, so that the
 * Solution::scale() and Solution::sum() operations are always possible;
 * specialized Solution for the specific cases are left for future
 * development. */

class MCFSolution : public Solution {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*------------------------------- FRIENDS ----------------------------------*/

 friend SingleFlowDCRBlock;  ///< make SingleFlowDCRBlock friend

/*---------------- CONSTRUCTING AND DESTRUCTING MCFSolution ----------------*/

 explicit MCFSolution( void ) { }  /// constructor, it has nothing to do

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 void deserialize( const netCDF::NcGroup & group ) override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 ~MCFSolution() = default;  ///< destructor: it is virtual, and empty

/*------------- METHODS DESCRIBING THE BEHAVIOR OF A MCFSolution -----------*/

 void read( const Block * block ) override final;

 void write( Block * block ) override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// serialize a MCFSolution into a netCDF::NcGroup
 /** Serialize a MCFSolution into a netCDF::NcGroup, with the following
  * format:
  *
  * - The dimension "NumNodes" containing the number of nodes. The dimension
  *   is optional, if it is not specified then the corresponding variable
  *   "Potentials" is not read (the MCFSolution object does not contain any
  *   node potentials).
  *
  * - The dimension "NumArcs" containing the number of arcs. The dimension
  *   is optional, if it is not specified then the corresponding variable
  *   "Potentials" is not read (the MCFSolution object does not contain any
  *   flow solution).
  *
  * - The variable "FlowSolution", of type double and indexed over the
  *   dimension NumArcs. The variable is optional, if it is not specified
  *   then the MCFSolution object does not contain any flow solution.
  *
  * - The variable "Potentials", of type double and indexed over the
  *   dimension NumNodes. The variable is optional, if it is not specified
  *   then the MCFSolution object does not contain any node potentials. */
 
 void serialize( netCDF::NcGroup & group ) const override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 MCFSolution * scale( double factor ) const override final;

 void sum( const Solution * solution , double multiplier ) override final;

 MCFSolution * clone( bool empty = false ) const override final;

/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/

 protected:

/*-------------------------- PROTECTED METHODS -----------------------------*/

 void print( std::ostream &output ) const override final {
  output << "MCFSolution [" << this << "]: " << v_x.size() << " flows and "
	 << v_pi.size() << " potentials" << std::endl;
  }

/*---------------------- PRIVATE PART OF THE CLASS -------------------------*/

 private:

/*---------------------------- PRIVATE FIELDS ------------------------------*/

 SingleFlowDCRBlock::Vec_FNumber v_x;   ///< the arc flows

 SingleFlowDCRBlock::Vec_CNumber v_pi;  ///< the node potentials

/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 };  // end( class( MCFSolution ) )

/** @} end( group( SingleFlowDCRBlock_CLASSES ) ) --------------------------*/
/*--------------------------------------------------------------------------*/

 };  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* SingleFlowDCRBlock.h included */

/*--------------------------------------------------------------------------*/
/*------------------- End File SingleFlowDCRBlock.h ------------------------*/
/*--------------------------------------------------------------------------*/
