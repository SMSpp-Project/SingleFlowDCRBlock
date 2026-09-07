/*--------------------------------------------------------------------------*/
/*------------------------ File MultiFlowDCRBlock.h ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the *concrete* class MultiFlowDCRBlock, which implements
 * the Block concept [see Block.h] for a Multicommodity Delay-Constrained
 * Routing (DCR) problem, i.e., the multi-flow extension of the single-flow
 * DCR problem implemented by SingleFlowDCRBlock [see SingleFlowDCRBlock.h].
 * Each commodity (flow) is represented by its own SingleFlowDCRBlock
 * sub-Block, and the sub-Blocks are coupled together by "mutual capacity"
 * constraints limiting the total rate that all the commodities can jointly
 * reserve on each arc.
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

#ifndef __MultiFlowDCRBlock
 #define __MultiFlowDCRBlock  /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "Block.h"
#include "SingleFlowDCRBlock.h"
#include "ColVariable.h"
#include "FRowConstraint.h"
#include "Configuration.h"
#include "Objective.h"

/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{
/*--------------------------------------------------------------------------*/
/*-------------------- MultiFlowDCRBlock-RELATED TYPES ---------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public Types
 *
 * "Import" basic types from SingleFlowDCRBlock.
 *
 *  @{ */

 using c_RHSValue = RowConstraint::c_RHSValue;  ///< a read-only RHS value

 using Vec_double = SingleFlowDCRBlock::Vec_double;    ///< a vector of double
 using c_Vec_double = SingleFlowDCRBlock::c_Vec_double;
 ///< a const vector of double

 using MultiVector = std::vector< Vec_double >;
 ///< a "matrix" (vector of vectors of double), one row per commodity

 using MultiSubset = std::vector< Block::Subset >;
 ///< a vector of Subset, one per commodity

 using Vec_Bool = std::vector< bool >;  ///< a vector of bool

/** @}  end( types ) */
/*--------------------------------------------------------------------------*/
/*------------------------------- CLASSES ----------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup MultiFlowDCRBlock_CLASSES Classes in MultiFlowDCRBlock.h
 *  @{ */

/*--------------------------------------------------------------------------*/
/*------------------------ CLASS MultiFlowDCRBlock -------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// implementation of the Block concept for the Multi-Flow DCR problem
/** The MultiFlowDCRBlock class implements the Block concept [see Block.h]
 * for the Multi-Flow Delay-Constrained Routing (DCR) problem, i.e., the
 * extension of the (single-flow) DCR problem implemented by
 * SingleFlowDCRBlock [see SingleFlowDCRBlock.h] to the case where several
 * flows ("commodities") k = 0, ..., NComm - 1 have to be simultaneously
 * routed on the same (directed) graph G = ( N , A ), with n = |N| nodes and
 * m = |A| arcs, while sharing the physical bandwidth of the arcs.
 *
 * Each commodity k is a DCR instance in its own right: it has its own
 * source/sink pair (encoded by its own node deficits), its own arc costs,
 * its own traffic parameters (burst, rate, deadline) and its own worst-case
 * delay bound to be respected (computed via network calculus, exactly as
 * in SingleFlowDCRBlock). Accordingly, MultiFlowDCRBlock does not implement
 * the whole per-commodity model itself: rather, for every commodity it
 * builds (see generate_abstract_variables()) a father-owned
 * SingleFlowDCRBlock sub-Block, which is charged with the per-commodity
 * flow-conservation, capacity and delay constraints and with the reserved
 * rate variables r^k[ i , j ] for each arc ( i , j ) of A. In principle an
 * alternative "knapsack" formulation, based on a BinaryKnapsackBlock per
 * commodity, is also part of the design (see the KnapsackRelaxation bit of
 * AR); however, the current implementation always uses the SingleFlowDCRBlock
 * ("flow relaxation") formulation, whatever the Configuration passed to
 * generate_abstract_variables() may say.
 *
 * What MultiFlowDCRBlock itself is responsible for is the *coupling* among
 * the commodities, in the form of "mutual capacity" constraints
 * \f[
 *  \sum_{ k = 0 }^{ NComm - 1 } r^k[ i , j ] \leq CapTot[ i , j ]
 *  \quad ( i , j ) \in A
 * \f]
 * stating that the sum of the rates reserved by all the commodities on a
 * given arc cannot exceed the (shared) capacity CapTot[] made available on
 * that arc; these are generated (see generate_abstract_constraints()) as a
 * static std::vector< FRowConstraint > MCs with one entry per arc, whose
 * active Variable are the r^k[ i , j ] Variable of the NComm
 * SingleFlowDCRBlock sub-Block. CapTot[] need not coincide with the "raw"
 * per-arc capacity UTot[] read from the instance: in load(), CapTot[] is
 * obtained from UTot[] by an oversubscription/discount factor (see the
 * static FACTOR constant in MultiFlowDCRBlock.cpp), reflecting the fact that
 * not all commodities are expected to simultaneously need their maximum
 * individual capacity.
 *
 * Since the whole per-commodity structure (Variable, flow-conservation and
 * delay Constraint, Objective) is delegated to the NComm SingleFlowDCRBlock
 * sub-Block, is_feasible() as implemented here only checks the mutual
 * capacity coupling Constraint MCs; feasibility of the individual
 * commodities has to be separately assessed on the corresponding
 * sub-Block. */

class MultiFlowDCRBlock : public Block
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/

 //enum MCFType { kMCF , kSPT };

/*--------------------------------------------------------------------------*/
/*--------------------- PUBLIC METHODS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 *  @{ */

 /// constructor of MultiFlowDCRBlock, taking a pointer to the father Block
 /** Constructor of MultiFlowDCRBlock. It accepts a pointer to the father
  * Block, which can be of any type, defaulting to nullptr so that this can
  * also be used as the void constructor. Note that the actual per-commodity
  * SingleFlowDCRBlock sub-Block are *not* created here, but rather by
  * load() / deserialize() (which read the instance data) followed by
  * generate_abstract_variables(). */

 MultiFlowDCRBlock( Block *father = nullptr ) : Block( father ) , AR( 0 ) { }

/*--------------------------------------------------------------------------*/
 /// destructor of MultiFlowDCRBlock
 /** Destructor of MultiFlowDCRBlock: deletes the abstract representation
  * (the mutual capacity Constraint), and destroys all the commodity
  * sub-Block, if any (see guts_of_destructor()). */

 virtual ~MultiFlowDCRBlock() { guts_of_destructor(); }

/** @} ---------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */

 /// loads the instance from the given file in the given format
 /** Loads a MMCF instance using filename as the "base filename". This method
  * supports several formats depending on \p frmt, that is case-insensitive.
  * In particular, for two single-file formats
  *
  * - frmt == 0 (default) or frmt == 'c': PPRN format
  *
  * - frmt == 's': Canad format
  *
  * it behaves just as the Block method (just open an ifstream and
  * dispatch it load( std::istream & ). However, it also supports 5
  * multi-file formats:
  *
  * - 'm': Mnetgen format
  * - 'p': Jones-Lustig PSP (product-specific problem) format
  * - 'o': Jones-Lustig OSP (origin-specific problem) format
  * - 'd': Jones-Lustig OSP (origin-destination problem) format
  * - 'u': same as 'd' but supply information is looked at in file
  *        input + ".od" rather than input + ".sup" as in all the
  *        other cases
  *
  * where input (prefixed as set by set_filename_prefix(), if any) is
  * completed by the appropriate suffixes ".nod", ".arc", ".mut", ".sup"
  * or ".od" to load different parts of the description of the MMCF
  * instance.
  *
  * TODO: properly document all the formats.
  *
  * If there is any Solver attached to this MultiFlowDCRBlock then a NBModification
  * (the "nuclear option") is issued.
  *
  * IMPLEMENTATION NOTE: the above describes the intended, general design.
  * The *current* implementation, however, ignores \p frmt entirely and only
  * supports one fixed multi-file, per-commodity DCR format: it reads the
  * global data (number of commodities/nodes/arcs, MTU) from input + ".nod",
  * the per-commodity source/sink/rate from input + ".sup" and the
  * burst/deadline from input + ".param", and, for each commodity, the arc
  * topology/individual-capacity data from input + ".arc" together with the
  * delay-related data from input + ".dcr"; the latter two are used to build
  * a temporary DIMACS-like description that is fed to a freshly created
  * SingleFlowDCRBlock sub-Block (via SingleFlowDCRBlock::load() and
  * SingleFlowDCRBlock::load_dcr()). The "mutual capacity" CapTot[] of each
  * arc is *not* read from file but computed from the individual arc
  * capacity UTot[] by discounting it by a fixed oversubscription FACTOR
  * (see MultiFlowDCRBlock.cpp). */

 void load( const std::string & input , char frmt = 0 ) override;

/*--------------------------------------------------------------------------*/
 /// load the MultiFlowDCRBlock out of an istream
 /** Load the MultiFlowDCRBlock out of an istream. This is intended to
  * handle the two single-file formats, i.e., Canad and PPRN (as opposed to
  * the multi-file ones only supported by load( const std::string & )).
  *
  * TODO: properly document the formats.
  *
  * If there is any Solver attached to this MultiFlowDCRBlock then a NBModification
  * (the "nuclear option") is issued.
  *
  * IMPLEMENTATION NOTE: this method is currently a stub (its body is
  * empty): no format is actually read and no Modification is issued. */

 void load( std::istream & input , char frmt = 0 ) override;

/*--------------------------------------------------------------------------*/
 /// extends Block::deserialize( netCDF::NcGroup )
 /** Extends Block::deserialize( netCDF::NcGroup ) to the specific format of
  * a MultiFlowDCRBlock. Besides what is managed by the serialize() method
  * of the base Block class, the group is expected to contain:
  *
  * - the dimension "NNodes" containing the number of nodes in the graph;
  *
  * - the dimension "NArcs" containing the number of arcs in the graph;
  *
  * - the dimension "NComm" containing the number of commodities (flows);
  *
  * - the dimension "NCnst" containing the number of arcs having a mutual
  *   capacity constraint (currently always == "NArcs");
  *
  * - the variable "SN", of type int and indexed over "NArcs", containing
  *   the starting node of each arc;
  *
  * - the variable "EN", of type int and indexed over "NArcs", containing
  *   the ending node of each arc;
  *
  * - the variable "Utot", of type double and indexed over "NArcs",
  *   containing the (individual) total capacity of each arc;
  *
  * - the variable "U", of type double and indexed over ( "NComm" , "NArcs" ),
  *   containing the per-commodity individual arc capacities;
  *
  * - the variable "B", of type double and indexed over ( "NComm" , "NNodes" ),
  *   containing the per-commodity node deficits;
  *
  * - the variable "C", of type double and indexed over ( "NComm" , "NArcs" ),
  *   containing the per-commodity arc costs.
  *
  * The "NNodes", "NArcs" and "NComm" dimensions and the "SN", "EN" variables
  * are mandatory. Note that this method only reads the *coupling* data
  * (topology, mutual capacities, and the raw per-commodity U/B/C matrices);
  * it does *not* construct the NComm SingleFlowDCRBlock sub-Block, which is
  * instead the job of generate_abstract_variables(). */

 void deserialize( const netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/
 /// simplifies the problem
 /** Performs various pre-processing of the data, trying to make the instance
  * more easily solvable. The parameters to be given are the following:
  *
  * IncUk , DecUk   => (>= 0) upper bounds on the increase and decrease of the
  *                    mutual capacities: may be Inf<double>() if unknown;
  *
  * IncUjk , DecUjk => (>= 0) same as above for single-commodity capacities;
  *
  * ChgDfct         => (>= 0) upper bound on the maximum change, in absolute
  *                    value, of the node deficits: it must be a finite
  *                    number, since it is used to generate "loose" but
  *                    finite individual capacities for arcs that have none;
  *
  * DecCsts         => (>= 0) upper bound on the decrease of arc Costs: must
  *                    be < Inf<double>().
  *
  * Giving tight bounds (0 is the best, obviously) may cause the preprocessor
  * to find more redundant coupling constraints, to squeeze down individual
  * arc capacities, to remove more unused arcs and in general to do a better
  * preprocessing; for instance, IncUjk == 0 allows PreProcess() to declare
  * un-existent (set the cost to Inf<double>()) any arc with 0 individual
  * capacity.
  *
  * For all k such that, after the pre-processing, the graph has only a source
  * and no (existing) arcs have a "real" capacity, the type of the subproblem
  * is set to kSPT: all other problem types are left unchanged.
  *
  * Important note: in order for PreProcess() to work, it has to be able to
  * guess at least an upper bound on the maximum quantity of each commodity
  * in the graph. In order to do that, *all arcs* with potentially *negative
  * costs* (ChgCsts is used to estimate that) must have a *finite capacity*.
  *
  * PreProcess() will also look for redundancy in the data structures (e.g.
  * identical costs/deficits/individual capacities for some commodities) and
  * eliminate them, thus possibly saving some memory.
  *
  * It can be called *only once*. The ideal would be that it is automatically
  * called after load(), deserialize() ecc. but this would not allow to set
  * the proper parameters, therefore it has to be done independently (if
  * ever).
  *
  * IMPLEMENTATION NOTE: this method is currently declared but not defined;
  * calling it will result in a link-time error. */

 void PreProcess( double IncUk = 0 , double DecUk = 0 ,
		  double IncUjk = 0 , double DecUjk = 0 ,
		  double ChgDfct = 0 , double DecCsts = 0 );

/*--------------------------------------------------------------------------*/
 /// generate the "abstract representation" of the Variable of the Block
 /** This method generates the "abstract representation" of the Variable of
  * the MultiFlowDCRBlock, and in fact it decides which formulation of the MMCF
  * problem is implemented. This is controlled by the parameter stvv. If stvv
  * is not nullptr and it is a SimpleConfiguration< int >, or if
  * f_BlockConfig->f_static_variables_Configuration is not nullptr and it is a
  * SimpleConfiguration< int >, then the f_value (an int) dictates which
  * MMCF formulation as follows:
  *
  * - [1]: the standard knapsack formulation in which get_NArcs()
  *   BinaryKnapsackBlock sub-Block are constructed, one for each commodity,
  *   and the flow constraints are handled in the father MultiFlowDCRBlock;
  *
  * - [0]: the standard flow formulation in which get_NComm() SingleFlowDCRBlock
  *   sub-Block are constructed, one for each commodity, and the
  *   linking constraints are handled in the father MultiFlowDCRBlock;
  *
  * - [other ones possibly to follow].
  *
  *  by default is considered the Flow relaxation
  *
  * IMPLEMENTATION NOTE: the current implementation always uses the "flow"
  * formulation, i.e., it unconditionally calls
  * generate_abstract_variables() on every commodity's SingleFlowDCRBlock
  * sub-Block (which must already exist in v_Block, one per commodity);
  * stvv (and f_BlockConfig->f_static_variables_Configuration) are read but
  * *not* actually used to select the formulation, so the "knapsack"
  * alternative described above is not yet available. */

 void generate_abstract_variables( Configuration * stvv = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// generate the mutual capacity Constraint of the MultiFlowDCRBlock
 /** This method first calls generate_abstract_constraints() on every
  * commodity's SingleFlowDCRBlock sub-Block (which constructs their own
  * per-commodity flow-conservation, bound and delay Constraint), and then,
  * unless this has already been done (see the HasMutual bit of AR),
  * generates the "mutual capacity" coupling Constraint MCs: a static
  * std::vector< FRowConstraint > with one entry per arc j, of the form
  * \f[
  *  -\infty \leq \sum_{ k = 0 }^{ NComm - 1 } r^k[ j ] \leq CapTot[ j ]
  * \f]
  * whose active Variable are the reserved-rate Variable r^k[ j ] of each
  * commodity's SingleFlowDCRBlock sub-Block (obtained via
  * SingleFlowDCRBlock::i2p_r()). The parameter stcc is currently ignored
  * by this method (it is only forwarded to the sub-Block). */

 void generate_abstract_constraints( Configuration * stcc = nullptr )
  override;

/*--------------------------------------------------------------------------*/
 /*!! not needed yet, the version of Block suffices so far
 void generate_objective( Configuration * objc = nullptr ) override;
 !!*/

/** @} ---------------------------------------------------------------------*/
/*---------- METHODS FOR PRINTING & SAVING THE MultiFlowDCRBlock -----------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for printing & saving the MultiFlowDCRBlock
 *  @{ */

 /// print the MultiFlowDCRBlock on an ostream with the given verbosity
 /** Print the MultiFlowDCRBlock on an ostream. So far vlvl is ignored and only very
  * basic information is printed.
  *
  * TODO: implement some verbosity level that produce output files in at
  *       least some of the single-file formats supported by load(); note
  *       that for multi-file formats, print( std::string & ) must be used.
  *
  * IMPLEMENTATION NOTE: currently the body of this method is entirely
  * commented out, hence it prints nothing at all. */

 void print( std::ostream & output , char vlvl = 0 ) const override;

/*--------------------------------------------------------------------------*/
/// extends Block::serialize( netCDF::NcGroup )
/** Extends Block::serialize( netCDF::NcGroup ) to the specific format of a
 * MultiFlowDCRBlock. See MultiFlowDCRBlock::deserialize(netCDF::NcGroup) for details of the
 * format of the created netCDF group; in short, besides calling
 * Block::serialize(), this writes the "NNodes", "NArcs", "NComm" and
 * "NCnst" dimensions, the "SN", "EN" and "Utot" variables (the graph
 * topology and the per-arc individual capacity), and the "U", "B" and "C"
 * matrix variables (indexed over commodity and, respectively, arc, node
 * and arc) holding the per-commodity data. */

 void serialize( netCDF::NcGroup & file ) const override;

/*--------------------------------------------------------------------------*/
 /// returns true if the mutual capacity Constraint are (approximately)
 /// satisfied
 /** Returns true if the current value of the reserved-rate Variable of all
  * the commodity sub-Block satisfies the mutual capacity Constraint MCs
  * (see generate_abstract_constraints()) within tolerance. The tolerance
  * and the type of violation (absolute vs. relative) are extracted, in
  * order, from:
  *
  * - fsbc, if it is a SimpleConfiguration< double > (only the tolerance,
  *   relative violation is assumed) or a
  *   SimpleConfiguration< std::pair< double , int > > (tolerance and,
  *   respectively, whether the violation is relative);
  *
  * - otherwise, f_BlockConfig->f_is_feasible_Configuration, tested against
  *   the same two types;
  *
  * - otherwise, the default tolerance 1e-6 and relative violation are used.
  *
  * Note that this method only checks the mutual capacity coupling
  * Constraint: it does *not* check the flow-conservation, bound or delay
  * feasibility of the individual commodities, which has to be separately
  * verified on each commodity's SingleFlowDCRBlock sub-Block. Also note
  * that \p useabstract is currently unused, as RowConstraint::is_feasible()
  * is always applied to the abstract Constraint MCs. */

 bool is_feasible( bool useabstract = false ,
                   Configuration * fsbc = nullptr ) override;

/** @} ---------------------------------------------------------------------*/
/*-------------- Methods for reading the data of the SingleFlowDCRBlock ----*/
/*--------------------------------------------------------------------------*/
/** @name Methods for reading the data of the MultiFlowDCRBlock
 *  @{ */

 /// get the number of nodes

 Index get_NNodes( void ) const { return( NNodes ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the number of arcs

 Index get_NArcs( void ) const { return( NArcs ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the number of commodities

 Index get_NComm( void ) const { return( NComm ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the vector of mutual (shared) arc capacities
 /** Returns a const reference to the vector CapTot of the *mutual* arc
  * capacities, i.e., the RHS of the coupling Constraint MCs (see
  * generate_abstract_constraints()). Note that this is *not* the protected
  * MultiVector U of per-commodity individual capacities. */

 c_Vec_double & get_U( void ) const { return( CapTot ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the vector of link (arc) delays, one entry per arc

 c_Vec_double & get_LinkDelays( void ) const { return( LinkDelays ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the vector of node delays, one entry per node

 c_Vec_double & get_NodeDelays( void ) const { return( NodeDelays ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the vector of the per-commodity Maximum Transmit Unit (MTU)

 c_Vec_double & get_MTU( void ) { return( MTU ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the vector of the per-commodity traffic burst

 c_Vec_double & get_FlowBurst( void ) { return( FlowBursts ); }

/*--------------------------------------------------------------------------*/
 /// returns true if the "flow" (as opposed to "knapsack") relaxation is used
 /** Returns true if the KnapsackRelaxation bit of AR is *not* set, i.e., if
  * the commodities are (as is currently always the case, see
  * generate_abstract_variables()) represented via SingleFlowDCRBlock
  * sub-Block rather than via BinaryKnapsackBlock sub-Block. */

 bool useFlowRelaxation( void ) const {
  return( ! ( AR & KnapsackRelaxation ) );
 }

/*--------------------------------------------------------------------------*/
 /// getting the current sense of the Objective (minimization by default)

 int get_objective_sense( void ) const override final {
  return( f_sense );
  }

/*--------------------------------------------------------------------------*/
 /// get the value of the reserved-rate variable r^k[ i ] of commodity k
 /** Returns the current value of the reserved-rate Variable r^k[ i ] of
  * commodity k on arc i, obtained from the corresponding
  * SingleFlowDCRBlock sub-Block (see SingleFlowDCRBlock::get_r()). */

  double get_rs( Index k , Index i ) const {
    return( static_cast< SingleFlowDCRBlock * >( v_Block[ k ] )->get_r( i ) );
  }

/*--------------------------------------------------------------------------*/
 /// get the flow of a given arc for a given commodity
 /** Given a commodity index k and an arc index ij, this function provides
  * the value of the associated variable x^k_ij. In the case of the knapsack
  * relaxation, the variables of the block are rescaled in such a way that
  * x \in [ 0 , 1 ]. In this case the functions get_flow provide the values
  * already rescaled wigth x^k_{ij} in [ 0 , u_ij ].
  *
  * IMPLEMENTATION NOTE: the current implementation always reads x^k_ij out
  * of the corresponding SingleFlowDCRBlock sub-Block (the "flow"
  * formulation); the AR-based dispatch outlined above, that would also
  * cover the "knapsack" formulation, is present but commented out. */

 double get_flow( Index k , Index i ) const {
  return( static_cast< SingleFlowDCRBlock * >( v_Block[ k ] )->get_x( i ) );
  /*
  if( ! ( AR & HasVar ) )
   return( 0 );

  if( ! ( AR & KnapsackRelaxation ) )
   return( static_cast< SingleFlowDCRBlock * >( v_Block[ k ] )->get_x( i ) );
  */
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the value of all flow variables associated to a given commodity
 /** Writes into \p fk (which must already be sized get_NArcs()) the value
  * of the flow Variable x^k_i for all the arcs i of commodity k. If the
  * abstract Variable have not been constructed yet (AR does not have the
  * HasVar bit set), \p fk is filled with zeroes instead. Note that, since
  * the "knapsack" formulation is currently never used (see
  * useFlowRelaxation()), the "if( ! ( AR & KnapsackRelaxation ) )" branch
  * is always taken whenever the Variable exist. */

 void get_flow( std::vector< double > & fk , Index k ) const {
  if( ! ( AR & HasVar ) ) {
   std::fill( fk.begin(), fk.end() , 0 );
   return;
   }

  if( ! ( AR & KnapsackRelaxation ) )
   static_cast<SingleFlowDCRBlock *>( v_Block[ k ] )->get_x( fk.begin() ,
						   Range( 0 , NArcs ) );
  }

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get a pointer to the ColVariable corresponding to the flow k , i
 /** Returns a pointer to the flow ColVariable x^k_i of commodity k on arc
  * i, obtained from the corresponding SingleFlowDCRBlock sub-Block (see
  * SingleFlowDCRBlock::i2p_x()); returns nullptr if the abstract Variable
  * have not been constructed yet (AR does not have the HasVar bit set), or
  * if the "knapsack" formulation is in use (which is currently never the
  * case, see useFlowRelaxation()). */

 ColVariable * get_flow_variable( Index k , Index i ) const {
  if( ! ( AR & HasVar ) )
   return( nullptr );

  if( ! ( AR & KnapsackRelaxation ) )
   return( static_cast< SingleFlowDCRBlock * >( v_Block[ k ] )->i2p_x( i ) );

  return( nullptr );
  }

/*--------------------------------------------------------------------------*/
 /// loads a MultiFlowDCRBlock out of a "traditional" SMS++ netCDF file
 /** Convenience method to load a MultiFlowDCRBlock out of a netCDF file
  * (as opposed to deserialize(), which works out of an already-open
  * netCDF::NcGroup): opens \p filename, reads the "SMS++_file_type"
  * attribute (currently unused beyond being read) and the "Block_0" group,
  * and calls deserialize() on it. */

 void load_nc4( std::string & filename ) {
  netCDF::NcFile f( filename, netCDF::NcFile::read );

  netCDF::NcGroupAtt gtype = f.getAtt( "SMS++_file_type" );

  int type = 0;
  gtype.getValues( &type );

  netCDF::NcGroup bg = f.getGroup( "Block_0" );

  deserialize( bg );

  //CmnIntlz();
  }


/** @} ---------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Protected methods for inserting and extracting
 *  @{ */


/*--------------------------------------------------------------------------*/
 /** called at the end of any constructor, does some initializations that are
  * common to them all: it is "protected" for allowing derived classes that
  * use the "void" constructor to call it.
  *
  * IMPLEMENTATION NOTE: currently unused (commented out), all the
  * initialization is done directly in load() / deserialize(). */

 //void CmnIntlz( void );

/** @} ---------------------------------------------------------------------*/
/*--------------------------- PROTECTED FIELDS  ----------------------------*/
/*--------------------------------------------------------------------------*/

 unsigned char AR;   ///< bit-wise coded: what abstract is there

 static constexpr unsigned char HasVar = 1;
 ///< first bit of AR == 1 if the formulation has been chosen already

 static constexpr unsigned char HasMutual = 2;
 ///< second bit of AR == 1 if the Mutual Constraints has been constructed

 static constexpr unsigned char KnapsackRelaxation = 4; 
 /**< third bit of AR == 1
   * - [1]: the standard knapsack formulation in which get_NArcs()
  *   BinaryKnapsackBlock sub-Block are constructed, one for each commodity,
  *   and the flow constraints are handled in the father MultiFlowDCRBlock;
  *
  * - [0]: the standard flow formulation in which get_NComm() SingleFlowDCRBlock
  *   sub-Block are constructed, one for each commodity, and the
  *   linking constraints are handled in the father MultiFlowDCRBlock;
  *
  * WHEN THE KNAPSACK RELAXATION IS CONSIDERED, THE PROVIDED FLOW SOLUTION
  * IS IN [ 0 , 1 ] TO OBTAIN THE SOLUTION OF THE INITIAL PROBLEM IS
  * NECESSARY TO RESCALE x^k_{ij} --> u_{ij} x^k_{ij}
  * the functions get_flow provides the value of the variable already
  * rescaled. 
  *
  * by default is considered the Flow relaxation
  */

 static constexpr unsigned char slc = 8;
 ///< fourth bit of AR == 1: true if we use the strong forcing constraints

 Index NNodes;         ///< Number of nodes
 Index NArcs;          ///< Number of arcs
 Index NComm;          ///< Number of commodities
 Index NCnst;          ///< Number of arcs with mutual capacity constraints

 MultiVector C;       ///< Matrix of the arc costs, one row per commodity
 MultiVector U;       ///< Matrix of the individual arc capacities, one
                       ///< row per commodity (not the mutual ones, see
                       ///< CapTot below)
 MultiVector B;       ///< Matrix of the node deficits, one row per
                       ///< commodity
 MultiVector I;       ///< Matrix of the variables integrality constraints
                       ///< (currently unused)

 Vec_double UTot;       ///< Vector of the "raw" (individual) total arc
                         ///< capacities, as read from the instance
 Vec_double CapTot;     ///< Vector of the actual *mutual* (shared) arc
                         ///< capacities enforced in the coupling
                         ///< Constraint MCs; computed from UTot by an
                         ///< oversubscription discount factor (see
                         ///< FACTOR in MultiFlowDCRBlock.cpp). Note that
                         ///< only UTot (not CapTot) is written/read by
                         ///< serialize()/deserialize(), so CapTot must be
                         ///< (re)computed after deserialize()

 Vec_double NodeDelays;     ///< Vector of the (fixed) per-node processing
                             ///< delays, one entry per node
 Vec_double LinkDelays;     ///< Vector of the (fixed) per-link propagation
                             ///< delays, one entry per arc

 Vec_double FlowBursts;     ///< Vector of the per-commodity traffic burst
                             ///< (network calculus arrival curve parameter)
 Vec_double FlowDeadlines;  ///< Vector of the per-commodity worst-case
                             ///< delay bound (deadline)
 Vec_double rho;            ///< Vector of the per-commodity sustained
                             ///< traffic rate
 Vec_double MTU;            ///< Vector of the per-commodity Maximum
                             ///< Transmit Unit

 Subset Startn;        ///< Topology of the graph: starting nodes
 Subset Endn;          ///< Topology of the graph: ending nodes

 Index StrtNme;        ///< The "name" of the first node

 std::vector< FRowConstraint > MCs;  ///< the static mutual capacity
                                      ///< Constraint, one per arc (see
                                      ///< generate_abstract_constraints())

 int f_sense = Objective::eMin;  ///< the sense of the Objective, as
                                  ///< returned by get_objective_sense()
                                  ///< (minimization by default)

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

 /// destroys the abstract representation and all the commodity sub-Block
 /** Clears the mutual capacity Constraint MCs, deletes all the commodity
  * SingleFlowDCRBlock sub-Block in v_Block, clears all the data structures
  * (C, U, B, I, UTot, CapTot, Startn, Endn, the delay/traffic vectors) and
  * resets AR to 0. Called by both the destructor and, to start from a
  * clean slate, by load() / deserialize(). */

 void guts_of_destructor( void );

/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE FIELDS ------------------------------*/
/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;  // insert it in the Block factory

/*--------------------------------------------------------------------------*/

 };  // end( class( MultiFlowDCRBlock ) )

/*--------------------------------------------------------------------------*/

/** @}  end( group( MultiFlowDCRBlock_CLASSES ) ) --------------------------*/
/*--------------------------------------------------------------------------*/

 }  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* MultiFlowDCRBlock.h included */

/*--------------------------------------------------------------------------*/
/*---------------------- End File MultiFlowDCRBlock.h ----------------------*/
/*--------------------------------------------------------------------------*/

