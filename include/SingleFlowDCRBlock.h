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
 * \author Enrico Sorbera \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni, Laura Galli, Luca Mencarelli,
 *                      Enrico Sorbera
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __SingleFlowDCRBlock
 #define __SingleFlowDCRBlock
                      /* self-identification: #endif at the end of the file */

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

#include <algorithm>

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

 using p_SingleFlowDCRBlock = SingleFlowDCRBlock *;  ///< a pointer to
                                                     ///< SingleFlowDCRBlock

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
/// implementation of the Block concept for the single-flow DCR problem
/** The SingleFlowDCRBlock class implements the Block concept [see Block.h]
 * for the Delay-Constrained Routing (DCR) problem relative to a *single*
 * flow, as opposed to MultiFlowDCRBlock [see MultiFlowDCRBlock.h] that deals
 * with several flows sharing the same network at once.
 *
 * A DCR problem [see DCR.h] consists of a Multi-Commodity Flow part plus a
 * max-delay constraint part expressed via network calculus. In the single-
 * flow case the "flow part" degenerates into finding a single s-t path (or,
 * more precisely, a 0/1 selection of arcs realizing that path) on a
 * (directed) graph G = ( N , A ) with n = |N| nodes and m = |A| arcs, while
 * the "delay part" requires the end-to-end delay incurred by the flow along
 * the selected arcs, as computed by a specific network-calculus formula
 * (currently the Strictly-Rate-Proportional one, see DCR::DCRDelay), not to
 * exceed a given deadline.
 *
 * For each arc (i, j) of A the following variables are defined:
 *
 * - a binary "routing" variable X[ i , j ] (field x), equal to 1 if and only
 *   if the arc is used by the flow;
 *
 * - a continuous "reserved rate" variable R[ i , j ] >= 0 (field r), the
 *   fraction of the arc bandwidth reserved to the flow, which can only be
 *   positive if the arc is selected;
 *
 * - a continuous "burst delay" variable \f$\Theta\f$[ i , j ] >= 0 (field
 *   theta), which is a convexification of the extra queueing delay incurred
 *   on the arc because of rate-limiting: for a Strictly-Rate-Proportional
 *   shaper this delay term is MTU / R[ i , j ], and forcing X[ i , j ] to
 *   multiply it out of the (nonconvex, hyperbolic) fraction leads to the
 *   convex constraint
 *   \f[
 *    \Theta[ i , j ] \, R[ i , j ] \geq MTU \, X[ i , j ]^2
 *   \f]
 *   (see cone_cnst), a rotated second-order cone in the three variables.
 *
 * In addition, two "aggregate" (not arc-specific) variables model the delay
 * contribution of the flow's own burstiness at the source:
 *
 * - the reserved rate r_min >= 0 (field r_min), the minimum reserved rate
 *   guaranteed to the flow along the whole selected path;
 *
 * - the burst delay theta_min >= 0 (field theta_min), convexifying the term
 *   FlowBursts / r_min via the analogous rotated cone constraint
 *   \f[
 *    \Theta_{min} \, R_{min} \geq FlowBursts
 *   \f]
 *   (see cone_min_cnst).
 *
 * The problem is then, denoting with U[ i , j ] the capacity of arc (i, j),
 * with C[ i , j ] its (per-unit-of-reserved-rate) cost, with LinkDelays and
 * NodeDelays the (fixed) propagation/processing delays of arcs and nodes,
 * and with FlowDeadlines the maximum admissible end-to-end delay:
 * \f[
 *  \min \sum_{ (i, j) \in A } C[ i , j ] R[ i , j ]
 * \f]
 * subject to the flow conservation constraints on the (single) commodity
 * (field E, encoding a 0/1 s-t path via the incidence matrix of G), the
 * delay constraint
 * \f[
 *  \sum_{ (i, j) \in A } \Theta[ i , j ] + X[ i , j ] \left(
 *   \frac{MTU}{U[ i , j ]} + LinkDelays[ i , j ] + NodeDelays[ i ] \right)
 *   + \Theta_{min} \; \leq \; FlowDeadlines
 * \f]
 * (field DCR_cnst), the two rotated-cone constraints above, and the "linking"
 * constraints (fields Indicator_cnst_r1, Indicator_cnst_r2 and
 * Indicator_cnst_rmin) that tie the reserved rate to the routing variables:
 * \f[
 *  0 \leq R[ i , j ] \leq U[ i , j ] X[ i , j ] \quad (i, j) \in A
 * \f]
 * \f[
 *  R[ i , j ] \geq \rho \, X[ i , j ] \quad (i, j) \in A
 * \f]
 * \f[
 *  R_{min} \leq R[ i , j ] \quad \mbox{for all (i, j) with } X[ i , j ] = 1
 * \f]
 * (the latter being implemented with a "big-M" reformulation using the
 * largest arc capacity U_max).
 *
 * Since the rotated-cone constraints above make this problem a Mixed-Integer
 * Second-Order Cone Program (the "SOCP" formulation), SingleFlowDCRBlock also
 * supports an alternative "P/C" (Perspective Cuts) formulation in which the
 * two cones are *not* explicitly constructed; rather, they are outer-
 * approximated by dynamically generated linear cuts (fields PC_cuts and
 * PC_cuts_min, populated by generate_dynamic_constraints()) separating the
 * current (x, r, theta) point whenever it violates the corresponding cone
 * by more than a given tolerance. Which of the two formulations is actually
 * constructed is controlled by a Configuration passed to
 * generate_abstract_constraints() [see there].
 *
 * The abstract representation can be changed (costs, capacities, deficits,
 * opening/closing arcs) via the corresponding methods, which keep the
 * "physical" and "abstract" representations of the SingleFlowDCRBlock in
 * sync and issue the appropriate Modification; more complex changes to the
 * structure of the graph are, at the time being, not supported and result
 * in exceptions being thrown. */

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
 * Unlike MCFBlock, SingleFlowDCRBlock does not need distinct types for flow
 * variables, arc costs and objective function values, since (unlike MCF) the
 * DCR problem has no "integrality" property to exploit: because of the
 * binary routing variables and the rotated-cone (or, in the P/C formulation,
 * cutting-plane) constraints linking them to the continuous reserved-rate
 * and burst-delay variables, the problem is never purely integral even if
 * all the input data are. Hence a single "generic" double type, and the
 * corresponding vector/iterator types, are used throughout for all the
 * numerical data (costs, capacities, deficits, delays, ...) of the
 * SingleFlowDCRBlock.
 @{ */

/*--------------------------------------------------------------------------*/

 typedef const double c_double;            ///< a read-only double

 typedef std::vector< double > Vec_double; ///< a vector of double

 typedef const Vec_double c_Vec_double;    ///< a const vector of double

 typedef Vec_double::iterator Vec_double_it;   ///< iterator in Vec_double

 typedef Vec_double::const_iterator c_Vec_double_it;
                                           ///< const iterator in Vec_double

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

 /// constructor of SingleFlowDCRBlock, taking a pointer to the father
 /// (generic) Block
 /** Constructor of SingleFlowDCRBlock. It accepts a pointer to the
  * father Block, which
  * can be of any type, defaulting to nullptr so that this can also be used as
  * the void constructor. */

 explicit SingleFlowDCRBlock( Block *father = nullptr )
  : Block( father ) , NNodes( 0 ) , NArcs( 0 ) , MaxNNodes( 0 ) , AR( 0 ) ,
    f_cond_lower( - Inf< double >() ) , f_cond_upper( - Inf< double >() ) { }

/*--------------------------------------------------------------------------*/
 /// destructor of SingleFlowDCRBlock: deletes the abstract representation, if
 /// any

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
  * - pNodeDelays  is the vector of the (fixed) processing delays of the
  *        nodes, which must either have size at least n or be empty, in the
  *        latter case all node delays are taken to be 0
  *
  * - pLinkDelays  is the vector of the (fixed) propagation delays of the
  *        arcs, which must either have size at least m or be empty, in the
  *        latter case all link delays are taken to be 0
  *
  * - FlowBursts  is the burst of the (single) flow to be routed, i.e., the
  *        amount of data that can be injected into the network instantly
  *        (used in the burst-delay term theta_min, see the general notes)
  *
  * - FlowDeadlines  is the maximum end-to-end delay that the flow can incur
  *        along the selected path (the RHS of the DCR delay constraint,
  *        see DCR_cnst)
  *
  * - MTU  is the Maximum Transmit Unit of the network, i.e., the size of the
  *        largest packet that can be sent over any arc (used in both the
  *        transmission-delay term MTU / U[ i , j ] and the per-arc burst-
  *        delay term theta[ i , j ], see the general notes)
  *
  * - rho  is the minimum rate that has to be reserved on any arc selected
  *        by the flow (the RHS of the r[ i , j ] >= rho * x[ i , j ]
  *        constraint, see Indicator_cnst_r2)
  *
  * Note that, unlike MCFBlock, SingleFlowDCRBlock does not have an explicit
  * vector of node deficits: the source and the sink of the (single) flow to
  * be routed are instead identified by the sign of the entries of the B[]
  * vector filled in by load( std::istream & ) / deserialize(), or set
  * directly with chg_dfct() / chg_st().
  *
  * Like load( std::istream & ), if there is any Solver attached to this
  * SingleFlowDCRBlock then a NBModification (the "nuclear option") is
  * issued. */

 void load( Index n , Index m , c_Subset & pSn , c_Subset & pEn ,
            c_Vec_double & pU = {} , c_Vec_double & pC = {} ,
            c_Vec_double & pNodeDelays = {} ,
            c_Vec_double & pLinkDelays = {} ,
            c_double FlowBursts = 0 , c_double FlowDeadlines = 0 ,
            c_double MTU = 0 , c_double rho = 0 );

/*--------------------------------------------------------------------------*/
 /// extends Block::deserialize( netCDF::NcGroup )
 /** Extends Block::deserialize( netCDF::NcGroup ) to the specific format of
  * a SingleFlowDCRBlock. Besides what is managed by the serialize() method
  * of the base Block class, the group should contain the following:
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
 /** Protected method for loading a SingleFlowDCRBlock out of a std::istream
  * (which is what operator>> is dispatched to. The std::istream is assumed
  * to contain the description of a DCR instance in DIMACS standard format,
  * which is  the following. The first line must be
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
  * Like load( memory ), if there is any Solver attached to this
  * SingleFlowDCRBlock then a NBModification (the "nuclear option") is
  * issued. */

 void load( std::istream &input , char frmt = 0 ) override;

/*--------------------------------------------------------------------------*/
 /// loads the DCR-specific data (delays, burst, deadline, MTU, rho)
 /** Reads from the given std::istream the data of the SingleFlowDCRBlock
  * that are specific to the DCR problem, i.e., are not part of the "plain"
  * MCF instance read by load( std::istream & ): in order, the NNodes node
  * delays, the NArcs link delays, the flow burst, the flow deadline, the
  * MTU and rho. This is meant to be called right after load( std::istream &
  * ) has read the "MCF part" of the instance (in DIMACS format) from the
  * same stream, to complete it with the DCR-specific parameters (typically
  * stored in a separate ".dcr" file); NNodes and NArcs must therefore match
  * those of the MCF part already read. Note that this method does *not*
  * issue any Modification. */

 void load_dcr( std::istream &input , Index NNodes , Index NArcs );

/*--------------------------------------------------------------------------*/
 /// generate the abstract variables of the DCR
 /** Method that generates the abstract Variable of the DCR. These are:
  *
  * - a std::vector< ColVariable > of size get_NArcs(), the binary routing
  *   variables X[ i , j ] (field x), one per arc, each equal to 1 if and
  *   only if the corresponding arc is used by the flow;
  *
  * - a std::vector< ColVariable > of size get_NArcs(), the continuous
  *   nonnegative reserved-rate variables R[ i , j ] (field r), one per arc;
  *
  * - a std::vector< ColVariable > of size get_NArcs(), the continuous
  *   nonnegative burst-delay variables \f$\Theta\f$[ i , j ] (field theta),
  *   one per arc;
  *
  * - the single continuous nonnegative ColVariable r_min, the minimum
  *   reserved rate guaranteed to the flow along the selected path;
  *
  * - the single continuous nonnegative ColVariable theta_min, the burst-
  *   delay term associated with r_min.
  *
  * See the general notes of the class for the exact meaning of these
  * variables. */

 void generate_abstract_variables( Configuration *stvv = nullptr ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// generate the static constraint of the DCR
 /** Method that generates the abstract Constraint of the DCR. These are:
  *
  * - if get_NNodes() > 0, a std::vector< FRowConstraint > with exactly
  *   get_NNodes() entries, the entry i being the flow conservation equation
  *   of node i (field E), exactly as in MCFBlock;
  *
  * - the single FRowConstraint DCR_cnst encoding the end-to-end delay bound
  *   of the flow, i.e.,
  *   \f[
  *    \sum_{ (i, j) \in A } \Theta[ i , j ] + X[ i , j ] \left(
  *     \frac{MTU}{U[ i , j ]} + LinkDelays[ i , j ] + NodeDelays[ i ]
  *     \right) + \Theta_{min} \; \leq \; FlowDeadlines ;
  *   \f]
  *
  * - depending on the Configuration stcc (see below), either the rotated
  *   second-order cone constraints cone_min_cnst and cone_cnst realizing
  *   the "SOCP" formulation of the burst-delay terms, i.e.,
  *   \f$\Theta_{min} R_{min} \geq FlowBursts\f$ and
  *   \f$\Theta[ i , j ] R[ i , j ] \geq MTU X[ i , j ]^2\f$
  *   (these are always *constructed*, but are only *added* to the abstract
  *   representation if the SOCP formulation is selected; the P/C
  *   formulation constructs an outer linearization of the same constraints
  *   instead, see generate_dynamic_constraints());
  *
  * - the three std::vector< FRowConstraint > Indicator_cnst_rmin,
  *   Indicator_cnst_r1 and Indicator_cnst_r2, of size get_NArcs() each,
  *   linking the reserved-rate variables to the routing ones:
  *   \f[
  *    R_{min} + U_{max} X[ i , j ] - R[ i , j ] \leq U_{max}
  *   \f]
  *   \f[
  *    R[ i , j ] \leq U[ i , j ] X[ i , j ]
  *   \f]
  *   \f[
  *    \rho \, X[ i , j ] \leq R[ i , j ]
  *   \f]
  *   where \f$U_{max}\f$ is the largest arc capacity; the first constraint
  *   is a "big-M" reformulation ensuring that, whenever some selected arc
  *   has X[ i , j ] == 1, its reserved rate R[ i , j ] is at least r_min.
  *
  * The Configuration stcc selects which of the two alternative formulations
  * of the burst-delay terms (SOCP vs P/C) is actually constructed: if
  *
  * - either stcc is not nullptr and it is a SimpleConfiguration< int >;
  *
  * - or f_BlockConfig is not nullptr and
  *   f_BlockConfig->f_static_variables_Configuration is not nullptr and it
  *   is a SimpleConfiguration< int >;
  *
  * then its f_value selects the formulation (1 == "P/C", 2 == "SOCP");
  * otherwise, the SOCP formulation is used by default. */

 void generate_abstract_constraints( Configuration *stcc = nullptr ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// generate the objective of the DCR
 /** Method that generates the (linear) Objective of the DCR, i.e.,
  * \f[
  *  \min \sum_{ (i, j) \in A } C[ i , j ] R[ i , j ] ,
  * \f]
  * the cost-weighted sum of the reserved-rate variables (the routing and
  * burst-delay variables do not appear in the objective). Unlike MCFBlock,
  * the choice between a "sparse" and a "dense" LinearFunction is currently
  * not exposed via a Configuration: the objective is always built as a
  * "dense" one, over all the get_NArcs() reserved-rate variables. */

 void generate_objective( Configuration *objc = nullptr ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// generate the dynamic (perspective) cuts for the "P/C" formulation
 /** This method implements the "P/C" (Perspective Cuts) alternative to the
  * rotated-cone constraints cone_cnst and cone_min_cnst [see
  * generate_abstract_constraints()]: rather than explicitly representing
  * the (convex, but nonlinear) constraints
  * \f[
  *  \Theta[ i , j ] \, R[ i , j ] \geq MTU \, X[ i , j ]^2
  *  \qquad \mbox{and} \qquad
  *  \Theta_{min} \, R_{min} \geq FlowBursts ,
  * \f]
  * this method inspects the *current* value of the (x, r, theta) and
  * (r_min, theta_min) variables and, whenever a cone is violated by more
  * than a tolerance tol, adds to PC_cuts (respectively, PC_cuts_min) the
  * linear cut obtained by linearizing the corresponding hyperbolic term
  * MTU X^2 / R (respectively, FlowBursts / R_min) around the current point,
  * i.e., its first-order Taylor approximation. Since the term is convex in
  * (X, R) for X binary, any such linear cut is valid (does not cut off any
  * feasible integer solution) and can be added on the fly, e.g., while
  * solving the LP relaxation of the P/C formulation with a cutting-plane
  * scheme. Routing variables whose value is below the tolerance eps are
  * skipped, as no cut is meaningful for them (the ratio X^2 / R being
  * ill-defined for R == 0).
  *
  * The tolerance tol and the "is-the-variable-really-nonzero" tolerance eps
  * can be provided via the Configuration stcc (or, failing that, via
  * f_BlockConfig->f_dynamic_constraints_Configuration): either a
  * SimpleConfiguration< double > (only changing tol) or a
  * SimpleConfiguration< std::pair< double , double > > (changing both tol
  * and eps, in this order); if none is provided, the default values
  * tol == 1e-5 and eps == 1e-4 are used.
  *
  * Whether the generated cuts are actually added to the abstract
  * representation (as opposed to being merely available in PC_cuts and
  * PC_cuts_min for later use) depends on the same formulation-selecting
  * Configuration used by generate_abstract_constraints(): the cuts are
  * added if and only if the "P/C" formulation is selected. */

 void generate_dynamic_constraints( Configuration *stcc = nullptr ) override;

/** @} ---------------------------------------------------------------------*/
/*-------- Methods for reading the data of the SingleFlowDCRBlock ----------*/
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
 /// get the Maximum Transmit Unit of the network

 [[nodiscard]] double get_MTU( void ) const { return( MTU ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the vector of the (fixed) node processing delays

 [[nodiscard]] c_Vec_double & get_NodeDelays( void ) const {
  return( NodeDelays );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the vector of the (fixed) arc/link propagation delays

 [[nodiscard]] c_Vec_double & get_LinkDelays( void ) const {
  return( LinkDelays );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the burst of the flow

 [[nodiscard]] double get_FlowBurst( void ) const { return( FlowBursts ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the (end-to-end) deadline of the flow

 [[nodiscard]] double get_FlowDeadline( void ) const {
  return( FlowDeadlines );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get rho, the minimum rate to be reserved on any arc used by the flow

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
 /// given a pointer to a reserved-rate Variable, returns the index of the arc
 /** Given a pointer to a reserved-rate Variable (formally a Variable *, but
  * immediately static_cast-ed to a ColVariable * right inside), returns the
  * index of the corresponding arc. Throws exception if the pointer is not to
  * a [Col]Variable of the SingleFlowDCRBlock. */

 [[nodiscard]] Index p2i_r( const Variable * var ) const {
  auto i = p2i_r_s( var );
  if( ( i >= 0 ) && ( i < int( get_NArcs() ) ) )
   return( i );

  throw( std::invalid_argument( "invalid arc pointer" ) );
  return( 0 );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// given an arc, returns the pointer to the corresponding reserved-rate
 /// Variable
 /** Given the index of an arc, returns the pointer to the corresponding
  * reserved-rate variable R[ i , j ] (a ColVariable *). This ASSUMES THE
  * Variable ARE CONSTRUCTED IN THE FIRST PLACE, SEGFAULTS ARE BOUND TO
  * HAPPEN OTHERWISE. */

 [[nodiscard]] ColVariable * i2p_r( Index i ) const {
  if( i < get_NArcs() )
   return( const_cast< ColVariable * >( & r[ i ] ) );

  throw( std::invalid_argument( "invalid arc pointer" ) );
  return( 0 );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns the pointer to the r_min ColVariable
 /** Returns the pointer to r_min, the ColVariable representing the minimum
  * reserved rate guaranteed to the flow along the selected path [see the
  * general notes of the class]. This ASSUMES THE Variable ARE CONSTRUCTED
  * IN THE FIRST PLACE, SEGFAULTS ARE BOUND TO HAPPEN OTHERWISE. */

 [[nodiscard]] ColVariable * p_rmin() const {
  return( const_cast< ColVariable * >( & r_min ) );
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
 *
 * Feasibility of a solution of the DCR problem is checked in four steps:
 * the routing variables have to describe a path from the source to the sink
 * [see flow_feasible()], the routing and reserved-rate variables have to
 * satisfy their bounds [see bound_feasible()] and the constraints linking
 * the two [see link_feasible()], and the end-to-end delay of the resulting
 * path must not exceed the deadline of the flow [see delay_feasible()].
 * is_feasible() checks the four together on the current solution, while
 * is_sol_feasible() does the same on the solution held by a DCRSolution.
 *
 * Each of the four checks comes in two versions: one taking the solution to
 * be checked from the outside, and one reading it out of the abstract
 * representation of the SingleFlowDCRBlock, either "physically" (out of the
 * value of the Variable) or "abstractly" (out of the corresponding
 * Constraint), according to the useabstract parameter. Note that the
 * "physical" version only reads the routing and reserved-rate variables
 * (x, r), the burst-delay ones (theta, theta_min) being *reconstructed*
 * from these: see delay_feasible( double , c_Vec_double & ,
 * c_Vec_double & ).
 *  @{ */

 /// returns true if the given routing is (approximately) flow feasible
 /** Returns true if the routing described by the value of the x Variable of
  * the SingleFlowDCRBlock approximately satisfies the flow conservation
  * constraints, i.e., it describes a path from the source to the sink of
  * the flow. This clearly requires the Variable of the SingleFlowDCRBlock
  * to have been defined, i.e., that generate_abstract_variables() has been
  * called prior to this method. The parameter feps is the relative accuracy
  * defining "approximately". The parameter "useabstract" has the same
  * meaning as in is_feasible(). */

 bool flow_feasible( double feps , bool useabstract = false );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// like flow_feasible(), but the routing is given from the outside
 /** Like flow_feasible( double , bool ), but the routing to be checked is X
  * rather than the one encoded in the Variable of the SingleFlowDCRBlock,
  * which therefore need not even exist. X must have size at least
  * get_NArcs(); the entries corresponding to deleted arcs are ignored. This
  * is what the "physical" version of flow_feasible() and is_sol_feasible()
  * both boil down to. */

 bool flow_feasible( double feps , c_Vec_double & X );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns true if the current solution is (approximately) bound feasible
 /** Returns true if the current value of the x and r Variable of the
  * SingleFlowDCRBlock approximately satisfies their bounds, i.e.,
  * 0 <= x[ i ] <= 1 and 0 <= r[ i ] <= U[ i ]. This clearly requires the
  * Variable of the SingleFlowDCRBlock to have been defined, i.e., that
  * generate_abstract_variables() has been called prior to this method. The
  * parameter feps is the relative accuracy defining "approximately". The
  * parameter "useabstract" has the same meaning as in is_feasible(). */

 bool bound_feasible( double feps , bool useabstract = false );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// like bound_feasible(), but the solution is given from the outside
 /** Like bound_feasible( double , bool ), but the solution to be checked is
  * the pair ( X , R ) rather than the one encoded in the Variable of the
  * SingleFlowDCRBlock, which therefore need not even exist. Both X and R
  * must have size at least get_NArcs(); the entries corresponding to
  * deleted arcs are ignored, while those corresponding to closed arcs are
  * checked against an upper bound of 0. */

 bool bound_feasible( double feps , c_Vec_double & X , c_Vec_double & R );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns true if routing and reserved rates are (approximately) linked
 /** Returns true if the current value of the x and r Variable of the
  * SingleFlowDCRBlock approximately satisfies the three families of
  * constraints linking the reserved rates to the routing decisions [see
  * generate_abstract_constraints()], i.e., rho x[ i ] <= r[ i ] <=
  * U[ i ] x[ i ] and r_min <= r[ i ] for every arc i used by the flow. The
  * parameter feps is the relative accuracy defining "approximately". The
  * parameter "useabstract" has the same meaning as in is_feasible(). */

 bool link_feasible( double feps , bool useabstract = false );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// like link_feasible(), but the solution is given from the outside
 /** Like link_feasible( double , bool ), but the solution to be checked is
  * the pair ( X , R ) rather than the one encoded in the Variable of the
  * SingleFlowDCRBlock, which therefore need not even exist. Both X and R
  * must have size at least get_NArcs(). Note that r_min is not part of the
  * pair: the value it is checked against is the largest feasible one, i.e.,
  * the smallest reserved rate among the arcs used by the flow, which makes
  * the corresponding constraints satisfied by construction. */

 bool link_feasible( double feps , c_Vec_double & X , c_Vec_double & R );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns true if the end-to-end delay (approximately) meets the deadline
 /** Returns true if the current value of the Variable of the
  * SingleFlowDCRBlock approximately satisfies the end-to-end delay
  * constraint DCR_cnst [see generate_abstract_constraints()]. The parameter
  * feps is the relative accuracy defining "approximately". The parameter
  * "useabstract" has the same meaning as in is_feasible(). */

 bool delay_feasible( double feps , bool useabstract = false );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// like delay_feasible(), but the solution is given from the outside
 /** Like delay_feasible( double , bool ), but the solution to be checked is
  * the pair ( X , R ) rather than the one encoded in the Variable of the
  * SingleFlowDCRBlock, which therefore need not even exist. Both X and R
  * must have size at least get_NArcs().
  *
  * Note that the pair does not comprise the burst-delay variables theta and
  * theta_min, which the delay constraint is written in terms of: these are
  * *reconstructed* out of ( X , R ) as the smallest values the two cone
  * constraints allow, i.e., theta[ i ] = MTU x[ i ] / r[ i ] and
  * theta_min = FlowBursts / r_min with r_min the smallest reserved rate
  * among the arcs used by the flow. Since the delay constraint is
  * monotonically increasing in theta and theta_min, ( X , R ) can be
  * completed into a solution that satisfies it if and only if the
  * reconstructed one does, which is what makes this the right question to
  * ask of a solution that only says which arcs are used and at which rate.
  */

 bool delay_feasible( double feps , c_Vec_double & X , c_Vec_double & R );

/*--------------------------------------------------------------------------*/
 /// returns true if the current solution is approximately feasible
 /** Returns true if the current solution, i.e., the current value of the
  * Variable of the SingleFlowDCRBlock, is approximately feasible: the
  * routing describes a path from the source to the sink of the flow
  * [see flow_feasible()], the routing and the reserved rates satisfy their
  * bounds [see bound_feasible()] and the constraints linking the two [see
  * link_feasible()], and the end-to-end delay of the flow does not exceed
  * its deadline [see delay_feasible()]. This clearly requires the Variable
  * of the SingleFlowDCRBlock to have been defined, i.e., that
  * generate_abstract_variables() has been called prior to this method (and,
  * if useabstract == true, that generate_abstract_constraints() has been
  * called as well).
  *
  * If useabstract == true the check is performed on the Constraint of the
  * abstract representation, which is only possible if they have been
  * constructed; otherwise, the value of the Variable is read and checked
  * against the "physical" data of the SingleFlowDCRBlock. The two are not
  * quite the same check: the abstract one takes the burst-delay variables
  * theta and theta_min at their face value, while the physical one
  * reconstructs them out of the reserved rates [see delay_feasible(
  * double , c_Vec_double & , c_Vec_double & )]. In particular, a solution
  * of the "P/C" formulation, in which the cone constraints are only
  * outer-approximated by finitely many linear cuts, can satisfy all the
  * Constraint that are there and still not be feasible for the DCR problem;
  * the physical check is the one that says so.
  *
  * The parameter for deciding what "approximately feasible" exactly means
  * is a single double value, the *relative* tolerance for the satisfaction
  * of all the constraints. This value is to be found as:
  *
  * - if fsbc is not nullptr and it is a SimpleConfiguration< double >, then
  *   it is fsbc->f_value;
  *
  * - otherwise, if f_BlockConfig is not nullptr,
  *   f_BlockConfig->f_is_feasible_Configuration is not nullptr and it is a
  *   SimpleConfiguration< double >, then it is
  *   f_BlockConfig->f_is_feasible_Configuration->f_value;
  *
  * - otherwise, it is 0. */

 bool is_feasible( bool useabstract = false , Configuration *fsbc = nullptr )
  override;

/*--------------------------------------------------------------------------*/
 /// returns true if the solution in the Solution is approximately feasible
 /** Returns true if the solution held by the DCRSolution sol is
  * approximately feasible, i.e., it satisfies the flow conservation
  * constraints, the bounds, the constraints linking the reserved rates to
  * the routing decisions and the end-to-end delay constraint; the
  * SingleFlowDCRBlock is only read for its data, its Variable are not
  * touched and they need not even exist. sol must be a DCRSolution holding
  * both the routing and the reserved-rate variables, otherwise false is
  * returned (a Solution that is not a DCRSolution is an error, and it
  * throws). The tolerance is found exactly as in is_feasible( bool ,
  * Configuration * ). */

 bool is_sol_feasible( Solution * sol ,
                       Configuration * fsbc = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// returns true if the DCR *instance* (not the current solution) is feasible
 /** Unlike is_feasible() and is_sol_feasible(), which check whether a
  * *solution* satisfies the constraints, this method checks whether the DCR
  * *instance* itself admits any feasible solution at all, i.e., whether
  * there exists some path from the (unique) source to the (unique) sink of
  * the flow (identified via the sign of the node deficits B[], see load())
  * whose end-to-end delay does not exceed FlowDeadlines.
  *
  * This is done by translating the SingleFlowDCRBlock data into a DCR::
  * DCRFlow / DCRLink / DCRNode description [see DCR.h] and handing it to a
  * DCR_SPT solver (which implements two Single-Flow Single-Path DCR
  * heuristics, "ERA-I" and "ERA-H", see DCR_SPT.h), using the Strictly-
  * Rate-Proportional delay formula (DCR::SRP) and the "ERA-I" heuristic
  * (DCRsetHeur( '1' )). The method returns false if and only if the
  * heuristic reports DCR::Infeasible, and true otherwise (including the
  * case in which the heuristic merely fails to find a feasible solution
  * that may still exist, since DCR_SPT is not an exact solver). */

 bool is_feasible_instance( void );

/** @} ---------------------------------------------------------------------*/
/*------------------------- Methods for R3 Blocks --------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for R3 Blocks
 *  @{ */

 /// gets an R3 Block of SingleFlowDCRBlock currently only the copy one
 /** Gets an R3 Block of the SingleFlowDCRBlock. The list of currently
  * supported R3 Block is:
  *
  * - r3bc == nullptr: the copy (an SingleFlowDCRBlock identical to this)
  */

 Block * get_R3_Block( Configuration *r3bc = nullptr ,
                       Block * base = nullptr , Block * father = nullptr )
  override;

/*--------------------------------------------------------------------------*/
 /// maps back the solution from a copy SingleFlowDCRBlock to the current one
 /** Maps back the solution from a copy SingleFlowDCRBlock to the current
  * one. The parameter r3bc is useless (has to be nullptr). The parameter
  * solc decides which part of the solution is mapped:
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
 /// maps the solution of the current SingleFlowDCRBlock to a copy
 /** Maps the solution of the current SingleFlowDCRBlock to a copy
  * SingleFlowDCRBlock. The parameter r3bc is useless (has to be nullptr).
  * The parameter solc decides which part of the solution is mapped:
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
  * Note that the current SingleFlowDCRBlock may not contain some or all of
  * the required solution, if the corresponding Variable/Constraint have not
  * been constructed yet: this throws an exception. */

 void map_forward_solution( Block *R3B , Configuration *r3bc = nullptr ,
                            Configuration *solc = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /** No specific Configuration is required, hence expected,
  *
  * IMPORTANT NOTE: map_forward_Modification() only maps "physical"
  * Modification. The point is that if any part of the "abstract
  * representation" of SingleFlowDCRBlock is changed, the corresponding
  * "abstract"
  * Modification is intercepted in add_Modification() and a "physical"
  * Modification is also issued. Hence, for any change in SingleFlowDCRBlock
  * there
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
  * IMPORTANT NOTE: SingleFlowDCRBlockRngdMod ALLOW TO ADD/DELETE ARCS IN THE
  * PROBLEM, WHICH ALSO CHANGES THE "NAMES" OF EXISTING ARCS.
  * SingleFlowDCRBlock
  *     IMPLEMENTS map_forward_Modification() IN A WAY THAT IS ONLY
  *     GUARANTEED TO BE CORRECT IF:
  *
  *     = EITHER THE SET OF ARCS IS NEVER CHANGED;
  *
  *     = OR THE Modification ARE MAPPED IMMEDIATELY AFTER THEY ARE ISSUED.
  *
  * This is because otherwise SingleFlowDCRBlock should have to understand
  * whether the
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
 /** No specific Configuration is required, hence expected
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

 /// returns as a DCRSolution the current solution of the SingleFlowDCRBlock
 /** Returns a new DCRSolution representing the current solution status of
  * this SingleFlowDCRBlock. The parameter solc decides which part of the
  * solution is retained: if it is not nullptr and it is a
  * SimpleConfiguration< int >, then its f_value == 1 means "only the
  * routing variables x[]" (v_r left empty), == 2 means "only the reserved-
  * rate variables r[]" (v_x left empty), and anything else (including
  * solc == nullptr) means "both". If emptys == true (the default) the
  * returned DCRSolution only has the right "shape" (size of v_x and/or
  * v_r), without actually copying the current values; otherwise, read() is
  * also called to fill it in with the current solution. */

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
  if( ! ( AR & HasVar ) )
   throw( std::logic_error( "flow Variable not available" ) );
  if( arc >= get_NArcs() )
   throw( std::invalid_argument( "invalid arc name" ) );

  return( x[ arc ].get_value() );
  }

/*--------------------------------------------------------------------------*/
 /// gets the current value of r_min, the minimum reserved rate

 double get_rmin( void ) const { return( r_min.get_value() ); }

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
  if( ! ( AR & HasVar ) )
   throw( std::logic_error( "reserved-rate Variable not available" ) );
  if( arc >= get_NArcs() )
   throw( std::invalid_argument( "invalid arc name" ) );

  return( r[ arc ].get_value() );
  }

/*--------------------------------------------------------------------------*/
 /// sets a contiguous interval of the reserve solution
 /** Method to set the reserve solution; the values found in the
  * c_Vec_double starting from fstrt are copied into the value of the
  * reserved-rate variable r[ i ] for i in rng, in the same order. */

  void set_r( c_Vec_double_it fstrt ,
              Range rng = Range( 0 , Inf< Index >() ) );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// sets a generic subset of the reserve solution
 /** Method to set the reserve solution; the values found in the
  * c_Vec_double starting from fstrt are copied into the value of the
  * reserved-rate variable r[ i ] for all i in sbst (that must be ordered in
  * increasing sense), in the same order. */

 void set_r( c_Vec_double_it fstrt , c_Subset sbst );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// sets the reserve solution of the given arc

 void set_r( Index arc , double FSol ) {
  if( ! ( AR & HasVar ) )
   throw( std::logic_error( "reserved-rate Variable not available" ) );
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
  if( ! ( AR & HasVar ) )
   throw( std::logic_error( "flow Variable not available" ) );
  if( arc >= get_NArcs() )
   throw( std::invalid_argument( "invalid arc name" ) );

  x[ arc ].set_value( FSol );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// sets the value of r_min, the minimum reserved rate

 void set_rmin( double FSol ) {
  r_min.set_value( FSol );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// sets the value of theta_min, the burst-delay term associated to r_min

 void set_theta_min( double FSol ) {
  theta_min.set_value( FSol );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// sets the value of theta[ arc ], the burst-delay term of the given arc

 void set_theta( Index arc , double FSol ) {
  if( ! ( AR & HasVar ) )
   throw( std::logic_error( "burst-delay Variable not available" ) );
  if( arc >= get_NArcs() )
   throw( std::invalid_argument( "invalid arc name" ) );

  theta[ arc ].set_value( FSol );
  }

/** @} ---------------------------------------------------------------------*/
/*-------------------- Methods for handling Modification -------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for handling Modification
 *  @{ */

 /// adding a new Modification to the SingleFlowDCRBlock
 /** Method for handling Modification.
  *
  * The version of SingleFlowDCRBlock has to intercept any "abstract
  * Modification" that modifies the "abstract representation" of the
  * SingleFlowDCRBlock, and "translate" them into both changes of the actual
  * data structures and corresponding "physical Modification". These
  * Modification are those for which Modification::concerns_Block() is true.
  * Note, however, that before sending the Modification to the Solver and/or
  * the father Block, the concerns_Block() value is set to false. This is
  * because once it is passed through this method, the "abstract Modification"
  * has "already done its duty" of providing the information to the
  * SingleFlowDCRBlock, and this must not be repeated. In particular, this
  * would be an issue if the Modification would be [map_forward or
  * map_back]-ed, because inside of this method a "physical Modification"
  * doing the same job is surely issued. That Modification would also be
  * [map_forward or map_back]-ed, together with the original "abstract
  * Modification" that would pass again through this method (in the other
  * SingleFlowDCRBlock), which would mean that the "physical Modification"
  * would be issued twice.
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
  * Any other Modification reaching the SingleFlowDCRBlock will lead to
  * exception being thrown.
  *
  * Note: any "physical" Modification resulting from processing an "abstract"
  *       one will be sent to the same channel (chnl). */

 void add_Modification( sp_Mod mod , ChnlName chnl = 0 ) override;

/** @} ---------------------------------------------------------------------*/
/*---------- METHODS FOR PRINTING & SAVING THE SingleFlowDCRBlock ----------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for printing & saving the SingleFlowDCRBlock
 *  @{ */

 /// print the SingleFlowDCRBlock on an ostream with the given verbosity
 /** Protected method to print information about the SingleFlowDCRBlock;
  * with the "complete" level ('C') it outputs the SingleFlowDCRBlock in
  * DIMACS format. */

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
 * concerns_Block() == true, since the changes in the SingleFlowDCRBlock have
 * surely been done already, and this is just not possible for a "physical"
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
 * the "abstract" representation. This allows to re-use them within
 * SingleFlowDCRBlock itself when reacting to abstract Modification, where
 * the  "abstract" representation has been changed already.
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
  * Also, if issueMod says so then a "physical" SingleFlowDCRBlockRngdMod is
  * issued. */

 void chg_costs( c_Vec_double_it NCost , Range rng = INFRange ,
                 ModParam issueMod = eNoBlck , ModParam issueAMod = eNoBlck );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// change the costs of an arbitrary subset of arcs
 /** Method to change the costs of an arbitrary subset of arc. That is,
  * *( NCost + i ) becomes the new cost of arc nms[ i ] for all 0 <= i <
  * NCost.size(), (which means that nms.size() == NCost.size()). The
  * parameter ordered tells if the nms vector is ordered for increasing
  * index of the arc. As the && tells, nms is "consumed" by the method,
  * typically being shipped to an appropriate SingleFlowDCRBlockSbstMod
  * object.
  *
  * See chg_costs( range ) for Modification issued (except that, of course,
  * the "physical" one is a SingleFlowDCRBlockSbstMod), and about changes in
  * costs of closed arcs. */

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
  * Also, if issueMod says so then a "physical" SingleFlowDCRBlockRngdMod is
  * issued. */

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
  * the "physical" one is a SingleFlowDCRBlockSbstMod) and about changing
  * capacities of closed arcs. */

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
  * Also, if issueMod says so then a "physical" SingleFlowDCRBlockRngdMod is
  * issued. */

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
  * typically being shipped to an appropriate SingleFlowDCRBlockSbstMod
  * object.
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
  * Also, if issueMod says so then a "physical" SingleFlowDCRBlockRngdMod is
  * issued. */

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
  * typically being shipped to an appropriate SingleFlowDCRBlockSbstMod
  * object.
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
  * Also, if issueMod says so then a "physical" SingleFlowDCRBlockRngdMod is
  * issued. */

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
  * typically being shipped to an appropriate SingleFlowDCRBlockSbstMod
  * object.
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
  * Also, if issueMod says so then a "physical" SingleFlowDCRBlockRngdMod is
  * issued. */

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
 /// changes the source and the sink of the flow
 /** Sets node ns as the (unique) source and node nt as the (unique) sink of
  * the flow to be routed: all node deficits are first reset to 0 and then
  * the deficit of ns is set to -1 and that of nt to +1 [see chg_dfct()],
  * hence issuing (up to) get_NNodes() + 2 separate Modification of type
  * SingleFlowDCRBlockRngdMod, one per changed node deficit; issueMod and
  * issueAMod have the same meaning as in chg_dfct(). */

 void chg_st( Index ns , Index nt , ModParam issueMod = eNoBlck ,
                    ModParam issueAMod = eNoBlck );


/*--------------------------------------------------------------------------*/
 /// adds a pair of static Constraint fixing r_min to the given value
 /** Adds to the SingleFlowDCRBlock two new static FRowConstraint,
  * "rmin_fixed1" (r_min >= value) and "rmin_fixed2" (r_min <= value), that
  * together fix the ColVariable r_min to the given value. Unlike
  * set_rmin(), which merely changes the *value* currently held by the
  * r_min ColVariable (without altering the feasible region), this method
  * permanently restricts the feasible region of the abstract representation
  * by adding two extra rows; it is meant to be used, e.g., to explore what
  * happens to the problem if r_min is pinned to a specific value. Note
  * that, unlike the other "changing the data" methods of this section, this
  * method neither checks for a Solver being attached nor issues any
  * Modification: the two Constraint are simply appended to the abstract
  * representation via add_static_constraint(). */

 void fix_rmin( double value ){
  auto rmin_fixed1 = new FRowConstraint();
  auto rmin_fixed2 = new FRowConstraint();

  LinearFunction::v_coeff_pair r1;
  r1.push_back( std::make_pair( &r_min , 1.0 ) );
  rmin_fixed1->set_lhs( value );
  rmin_fixed1->set_rhs( Inf<double>() );
  LinearFunction* Funct1 = new LinearFunction( std::move( r1 ));
  rmin_fixed1->set_function( Funct1 );

  LinearFunction::v_coeff_pair r2;
  r2.push_back( std::make_pair( &r_min , 1.0 ) );
  rmin_fixed2->set_lhs( -Inf<double>() );
  rmin_fixed2->set_rhs( value );
  LinearFunction* Funct2 = new LinearFunction( std::move( r2 ));
  rmin_fixed2->set_function( Funct2 );

  add_static_constraint( *rmin_fixed1 , "rmin_fixed1" );
  add_static_constraint( *rmin_fixed2 , "rmin_fixed2" );
  }

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
 Vec_double B;                  ///< vector of node deficits (source/sink)

 Vec_double NodeDelays;    ///< vector of the (fixed) node processing delays
 Vec_double LinkDelays;    ///< vector of the (fixed) arc propagation delays
 double FlowBursts;        ///< the burst of the flow
 double FlowDeadlines;     ///< the end-to-end deadline of the flow
 double MTU;                ///< the Maximum Transmit Unit of the network
 double rho;   ///< minimum rate to be reserved on any arc used by the flow

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

 std::vector< ColVariable > x;     ///< the binary routing variables X[i,j]
 std::vector< ColVariable > r;     ///< the reserved-rate variables R[i,j]
 std::vector< ColVariable > theta; ///< the per-arc burst-delay variables

 ColVariable r_min;      ///< the minimum reserved rate along the path
 ColVariable theta_min;  ///< the burst-delay term associated with r_min

 std::vector< FRowConstraint> E;   ///< the static flow conservation constrs.
 std::vector< LB0Constraint > UB;
 /**< the static bound constraints on the flow; NOTE: unlike in MCFBlock,
  * generate_abstract_constraints() does currently *not* construct/add UB
  * (the arc-capacity link between X[i,j] and R[i,j] being instead enforced
  * via Indicator_cnst_r1), even though it does set the AR & HasBnd bit; the
  * chg_ucap[s]() methods will still try to update UB[] whenever that bit is
  * set, which requires UB to have been sized elsewhere first. */

 FRowConstraint DCR_cnst;  ///< the end-to-end DCR delay constraint

 std::vector< FRowConstraint > Indicator_cnst_rmin;
 ///< "big-M" constraints linking r_min to R[i,j] on arcs used by the flow

 std::vector< FRowConstraint > Indicator_cnst_r1;
 ///< constraints R[i,j] <= U[i,j] X[i,j] (reserved rate needs an open arc)

 std::vector< FRowConstraint > Indicator_cnst_r2;
 ///< constraints R[i,j] >= rho X[i,j] (minimum rate on arcs used by the flow)

 FRowConstraint cone_min_cnst;
 ///< rotated cone theta_min * r_min >= FlowBursts ("SOCP" formulation)

 std::vector< FRowConstraint > cone_cnst;
 ///< rotated cones theta[i,j] * r[i,j] >= MTU * X[i,j]^2 ("SOCP" formulation)

 std::list< FRowConstraint > PC_cuts;
 ///< dynamically separated perspective cuts outer-approximating cone_cnst

 std::list< FRowConstraint > PC_cuts_min;
 ///< dynamically separated perspective cuts outer-approx. cone_min_cnst

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
 * because it does the registration of the following SingleFlowDCRBlock
 * methods:
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

  register_method< SingleFlowDCRBlock , MF_dbl_it , Range >(
                                          "SingleFlowDCRBlock::chg_costs" ,
                                          & SingleFlowDCRBlock::chg_costs );

  register_method< SingleFlowDCRBlock , MF_dbl_it , Subset && , bool >(
   "SingleFlowDCRBlock::chg_costs" , & SingleFlowDCRBlock::chg_costs );

  register_method< SingleFlowDCRBlock , MF_dbl_it , Range >(
                                          "SingleFlowDCRBlock::chg_ucaps" ,
                                          & SingleFlowDCRBlock::chg_ucaps );

  register_method< SingleFlowDCRBlock , MF_dbl_it , Subset && , bool >(
   "SingleFlowDCRBlock::chg_ucaps" , & SingleFlowDCRBlock::chg_ucaps );

  register_method< SingleFlowDCRBlock , Range >(
                                       "SingleFlowDCRBlock::close_arcs" ,
                                       & SingleFlowDCRBlock::close_arcs );

  register_method< SingleFlowDCRBlock , Subset && , bool >(
                                       "SingleFlowDCRBlock::close_arcs" ,
                                       & SingleFlowDCRBlock::close_arcs );

  register_method< SingleFlowDCRBlock , Range >(
                                       "SingleFlowDCRBlock::open_arcs" ,
                                       & SingleFlowDCRBlock::open_arcs );

  register_method< SingleFlowDCRBlock , Subset && , bool >(
                                       "SingleFlowDCRBlock::open_arcs" ,
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
 /// which of the two formulations of the burst-delay terms is in use
 /** Returns PCuts or SOCP according to what stcc says, or, failing that,
  * what the BlockConfig says in the Configuration of the static Constraint
  * (or in that of the static Variable, where the formulation used to be
  * asked for); the SOCP one, which is the exact formulation of the problem,
  * if nobody says anything. */

 Index formulation( Configuration * stcc = nullptr ) const;

/*--------------------------------------------------------------------------*/
 /// extracts out of fsbc, or of the BlockConfig, the is_feasible() tolerance

 void feps_of( Configuration * fsbc , double & feps , bool & rel_viol ) const;

/*--------------------------------------------------------------------------*/
 /// throws if X or R are too small to describe a solution of this Block

 void check_sizes( c_Vec_double & X , c_Vec_double & R ,
                   const char * method ) const;

/*--------------------------------------------------------------------------*/
 /// the largest r_min compatible with the solution ( X , R )
 /** The smallest reserved rate among the arcs that the routing X uses,
  * which is the largest value the Indicator_cnst_rmin constraints allow
  * r_min to take, and therefore the one that makes the burst-delay term
  * FlowBursts / r_min of the delay constraint as small as possible; 0 if
  * the flow uses no arc at all. */

 double min_rate( c_Vec_double & X , c_Vec_double & R ) const;

/*--------------------------------------------------------------------------*/
 /// true if lhs <= rhs is violated by more than the tolerance feps
 /** The tolerance is relative to rhs, save that a rhs smaller than \p scale
  * (one by default) counts as \p scale, so that the check does not become
  * an exact one where the right-hand side happens to be zero: \p scale is
  * the size of the terms of the constraint, which a Solver measures its own
  * tolerance against, e.g., the capacity U of r <= U x. */

 static bool violated( c_double lhs , c_double rhs , c_double feps ,
                       c_double scale = 1 ) {
  return( lhs - rhs > feps * std::max( { double( 1 ) , std::abs( scale ) ,
                                         std::abs( rhs ) } ) );
  }

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
/*--------------------- CLASS SingleFlowDCRBlockMod ------------------------*/
/*--------------------------------------------------------------------------*/
/// derived class from Modification for modifications to a SingleFlowDCRBlock
/** Derived class from Modification to describe modifications to a
 * SingleFlowDCRBlock. This is actually "sort of abstract", since it does not
 * say exactly what is changed, this being demanded to derived classes (which
 * do this in different ways). Note that it is derived from Modification
 * rather than,  say, BlockMod (which has the same structure) because this is
 * a class of "physical Modification". This means that a
 * SingleFlowDCRBlockMod refers to changes in the "physical representation"
 * of the SingleFlowDCRBlock; the corresponding changes in the "abstract
 * representation" of the SingleFlowDCRBlock are dealt with  by means of
 * "abstract Modification", i.e., derived classes from AModification (as is
 * BlockMod, which is why SingleFlowDCRBlockMod is not derived from
 * BlockMod). */

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

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 virtual ~SingleFlowDCRBlockMod() = default;   ///< destructor, does nothing

/*---------------------- PUBLIC METHODS OF THE CLASS -----------------------*/

 /// returns the [DCR]Block to which the SingleFlowDCRBlockMod refers

 Block * get_Block( void ) const override  { return( f_Block ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
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
               ///< pointer to the SingleFlowDCRBlock to which the
               ///< SingleFlowDCRBlockMod refers

 int f_type;   ///< type of Modification

/*--------------------------------------------------------------------------*/

 };  // end( class( SingleFlowDCRBlockMod ) )

/*--------------------------------------------------------------------------*/
/*------------------- CLASS SingleFlowDCRBlockRngdMod ----------------------*/
/*--------------------------------------------------------------------------*/
/// derived from SingleFlowDCRBlockMod for "ranged" modifications
/** Derived class from SingleFlowDCRBlockMod to describe "ranged"
 * modifications to a SingleFlowDCRBlock, i.e., modifications that apply to
 * an interval of either arcs or nodes. */

class SingleFlowDCRBlockRngdMod : public SingleFlowDCRBlockMod
{
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/

 public:

/*---------------------- CONSTRUCTOR & DESTRUCTOR --------------------------*/

 /// constructor: takes the SingleFlowDCRBlock, the type, and the range

 SingleFlowDCRBlockRngdMod( SingleFlowDCRBlock * fblock , int type ,
                            Block::Range rng )
  : SingleFlowDCRBlockMod( fblock , type ) , f_rng( rng ) {}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 virtual ~SingleFlowDCRBlockRngdMod() = default;
 ///< destructor, does nothing

/*---------------------- PUBLIC METHODS OF THE CLASS -----------------------*/

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
/*------------------- CLASS SingleFlowDCRBlockSbstMod ----------------------*/
/*--------------------------------------------------------------------------*/
/// derived from SingleFlowDCRBlockMod for "subset" modifications
/** Derived class from Modification to describe "subset" modifications to a
 *  SingleFlowDCRBlock, i.e., modifications that apply to an arbitrary
 * subset of either the arcs or the nodes. */

class SingleFlowDCRBlockSbstMod : public SingleFlowDCRBlockMod
{
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/

 public:


/*---------------------- CONSTRUCTOR & DESTRUCTOR --------------------------*/

 ///< constructor: takes the SingleFlowDCRBlock, the type, and the subset
 /**< Constructor: takes the SingleFlowDCRBlock, the type, and the subset.
  * As the the && tells, nms is "consumed" by the constructor and its
  * resources become property of the SingleFlowDCRBlockSbstMod object.
  *
  *   NOTE THAT nms IS REQUIRED TO BE ORDERED IN INCREASING SENSE
  *
  * although this is not checked by the class. */

 SingleFlowDCRBlockSbstMod( SingleFlowDCRBlock * fblock , int type ,
                            Block::Subset && nms )
  : SingleFlowDCRBlockMod( fblock , type ) , f_nms( std::move( nms ) ) {}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 virtual ~SingleFlowDCRBlockSbstMod() = default;
 ///< destructor, does nothing

/*---------------------- PUBLIC METHODS OF THE CLASS -----------------------*/

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
/** The DCRSolution class, derived from Solution, represents a (primal)
 * solution of a SingleFlowDCRBlock: the value of the binary routing
 * variables X[i,j] (field v_x) and of the continuous reserved-rate
 * variables R[i,j] (field v_r), one entry per arc. Either vector can be
 * left empty to represent "only part of the solution is saved", exactly as
 * in MCFSolution [see MCFBlock.h]; the aggregate (r_min, theta_min)
 * variables and the burst-delay variables theta[i,j] are currently *not*
 * part of a DCRSolution. */

class DCRSolution : public Solution
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*------------------------------- FRIENDS ----------------------------------*/

 friend SingleFlowDCRBlock;  ///< make SingleFlowDCRBlock friend

/*---------------- CONSTRUCTING AND DESTRUCTING DCRSolution ----------------*/

 explicit DCRSolution( void ) { }  /// constructor, it has nothing to do

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// de-serialize a DCRSolution out of a netCDF::NcGroup
 /** Extends Solution::deserialize( netCDF::NcGroup ) to the specific format
  * of a DCRSolution. The group may contain the dimension "NumArcs" and, if
  * so, the variables "ArcSolution" (the routing variables, into v_x) and/or
  * "FlowSolution" (the reserved-rate variables, into v_r), both indexed
  * over "NumArcs"; either (or both) can be missing, in which case the
  * corresponding vector is left empty. */

 void deserialize( const netCDF::NcGroup & group ) override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 ~DCRSolution() = default;  ///< destructor: it is virtual, and empty

/*------------- METHODS DESCRIBING THE BEHAVIOR OF A DCRSolution -----------*/
 /// reads the solution from the given SingleFlowDCRBlock
 /** Reads into v_x and v_r (whichever of the two is not empty) the current
  * value of, respectively, the routing variables x[] and the reserved-rate
  * variables r[] of the SingleFlowDCRBlock pointed by block, which must
  * actually be a SingleFlowDCRBlock (an exception is thrown otherwise). */

 void read( const Block * block ) override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// writes the solution into the given SingleFlowDCRBlock
 /** Writes the value of (whichever of) v_x and v_r (is not empty) into,
  * respectively, the routing variables x[] and the reserved-rate variables
  * r[] of the SingleFlowDCRBlock pointed by block, which must actually be a
  * SingleFlowDCRBlock (an exception is thrown otherwise); an exception is
  * also thrown if the size of v_x or v_r does not match get_NArcs(). */

 void write( Block * block ) override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// serialize a DCRSolution into a netCDF::NcGroup

 void serialize( netCDF::NcGroup & group ) const override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns a new DCRSolution with v_r scaled by the given factor
 /** Returns a new (empty, i.e., with the same "shape" but no meaningful
  * value in v_x) DCRSolution obtained by copying v_x unchanged and scaling
  * every entry of v_r by factor; note that the routing variables, being
  * binary, are *not* scaled. */

 DCRSolution * scale( double factor ) const override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// adds (a multiple of) another DCRSolution to this one
 /** Combines this DCRSolution with the one pointed by solution (which must
  * actually be a DCRSolution, and have vectors of the same size, otherwise
  * exception is thrown): each entry of v_x is replaced with the maximum of
  * the two corresponding entries (a reasonable way of "summing" 0/1
  * routing decisions), while each entry of v_r is incremented by the
  * corresponding entry of solution scaled by multiplier. */

 void sum( const Solution * solution , double multiplier ) override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// clones this DCRSolution
 /** Returns a new DCRSolution with the same "shape" (which of v_x and v_r
  * are non-empty, and their size) as this one; if empty == true the content
  * of v_x and v_r is not copied (only their size is replicated), otherwise
  * an exact copy of this DCRSolution is returned. */

 DCRSolution * clone( bool empty = false ) const override final;

/*----------- METHODS FOR READING AND WRITING THE SOLUTION -----------------*/
 /// returns the routing variables saved in this DCRSolution
 /** Returns the value of the routing variables x[] saved in this
  * DCRSolution, which is empty if it does not save them [see
  * SingleFlowDCRBlock::get_Solution()]. */

 [[nodiscard]] SingleFlowDCRBlock::c_Vec_double & get_x( void ) const {
  return( v_x );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns the reserved-rate variables saved in this DCRSolution
 /** Returns the value of the reserved-rate variables r[] saved in this
  * DCRSolution, which is empty if it does not save them [see
  * SingleFlowDCRBlock::get_Solution()]. */

 [[nodiscard]] SingleFlowDCRBlock::c_Vec_double & get_r( void ) const {
  return( v_r );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// sets the routing variables saved in this DCRSolution
 /** Sets the value of the routing variables x[] saved in this DCRSolution.
  * This is what a Solver fills the Solution with directly out of its own
  * data structures, rather than writing the solution in the Variable of the
  * SingleFlowDCRBlock and having it read back from there, which requires
  * the Variable to exist at all [see
  * SingleFlowDCRBendersSolver::get_Solution()]. */

 void set_x( SingleFlowDCRBlock::Vec_double && x ) { v_x = std::move( x ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// sets the reserved-rate variables saved in this DCRSolution
 /** The counterpart of set_x() for the reserved rates. */

 void set_r( SingleFlowDCRBlock::Vec_double && r ) { v_r = std::move( r ); }

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

 SingleFlowDCRBlock::Vec_double v_x;  ///< the arc routing (binary) variables

 SingleFlowDCRBlock::Vec_double v_r;  ///< the arc reserved-rate variables

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
