/*--------------------------------------------------------------------------*/
/*---------------------------- File MultiFlowDCRBlock.h ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the *concrete* class MultiFlowDCRBlock, which implements the
 * Block concept [see Block.h] for a Multicommodity Min Cost Flow problem.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Enrico Gorgone \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Francesco Demelas \n
 *         Laboratoire d'Informatique de Paris Nord \n
 *         Universite' Sorbonne Paris Nord \n
 *
 * Copyright &copy by Antonio Frangioni, Enrico Gorgone, Francesco Demelas
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
/*----------------------- MultiFlowDCRBlock-RELATED TYPES --------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public Types
 *
 * "Import" basic types from SingleFlowDCRBlock.
 *
 *  @{ */

 using CNumber = SingleFlowDCRBlock::CNumber;
 using c_RHSValue = RowConstraint::c_RHSValue;
 using Vec_CNumber = SingleFlowDCRBlock::Vec_CNumber;
 using FNumber = SingleFlowDCRBlock::FNumber;
 using Vec_FNumber = SingleFlowDCRBlock::Vec_FNumber;
 using c_Vec_FNumber = SingleFlowDCRBlock::c_Vec_FNumber;

 using FMultiVector = std::vector< Vec_FNumber >;
 using CMultiVector = std::vector< Vec_CNumber >;
 using MultiSubset = std::vector< Block::Subset >;

 using Vec_Bool = std::vector< bool >;

/** @}  end( types ) */
/*--------------------------------------------------------------------------*/
/*------------------------------- CLASSES ----------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup MultiFlowDCRBlock_CLASSES Classes in MultiFlowDCRBlock.h
 *  @{ */

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS MultiFlowDCRBlock ------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// Implementation of a simple MMCF Block concept.

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
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 *  @{ */

 /// constructor of MultiFlowDCRBlock
 /** Constructor of MultiFlowDCRBlock. It accepts a pointer to the father
  * Block, which can be of any type. */

 MultiFlowDCRBlock( Block *father = nullptr ) : Block( father ) , AR( 0 ) { }

/*--------------------------------------------------------------------------*/
 /// destructor of MultiFlowDCRBlock

 virtual ~MultiFlowDCRBlock() { guts_of_destructor(); }

/*@} -----------------------------------------------------------------------*/
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
  * (the "nuclear option") is issued. */

 void load( const std::string & input , char frmt = 0 ) override;

/*--------------------------------------------------------------------------*/
 /// load the MultiFlowDCRBlock out of an istream
 /** Load the MultiFlowDCRBlock out of an istream. Handles the two single-file
  * formats, i.e., Canad and PPRN.
  *
  * TODO: properly document the formats.
  *
  * If there is any Solver attached to this MultiFlowDCRBlock then a NBModification
  * (the "nuclear option") is issued. */

 void load( std::istream & input , char frmt = 0 ) override;

/*--------------------------------------------------------------------------*/
 /// extends Block::deserialize( netCDF::NcGroup )
 /** Extends Block::deserialize( netCDF::NcGroup ) to the specific format of
  * a MultiFlowDCRBlock.  */

 void deserialize( const netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/
 /// simplifies the problem
 /** Performs various pre-processing of the data, trying to make the instance
  * more easily solvable. The parameters to be given are the following:
  *
  * IncUk , DecUk   => (>= 0) upper bounds on the increase and decrease of the
  *                    mutual capacities: may be Inf<FNumber>() if unknown;
  *
  * IncUjk , DecUjk => (>= 0) same as above for single-commodity capacities;
  *
  * ChgDfct         => (>= 0) upper bound on the maximum change, in absolute
  *                    value, of the node deficits: it must be a finite
  *                    number, since it is used to generate "loose" but
  *                    finite individual capacities for arcs that have none;
  *
  * DecCsts         => (>= 0) upper bound on the decrease of arc Costs: must
  *                    be < Inf<CNumber>().
  *
  * Giving tight bounds (0 is the best, obviously) may cause the preprocessor
  * to find more redundant coupling constraints, to squeeze down individual
  * arc capacities, to remove more unused arcs and in general to do a better
  * preprocessing; for instance, IncUjk == 0 allows PreProcess() to declare
  * un-existent (set the cost to Inf<CNumber>()) any arc with 0 individual
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
  * ever). */

 void PreProcess( FNumber IncUk = 0 , FNumber DecUk = 0 ,
		  FNumber IncUjk = 0 , FNumber DecUjk = 0 ,
		  FNumber ChgDfct = 0 , CNumber DecCsts = 0 );

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
  */

 void generate_abstract_variables( Configuration * stvv = nullptr ) override;

/*--------------------------------------------------------------------------*/

 void generate_abstract_constraints( Configuration * stcc = nullptr )
  override;

/*--------------------------------------------------------------------------*/
 /*!! not needed yet, the version of Block suffices so far
 void generate_objective( Configuration * objc = nullptr ) override;
 !!*/

/** @} ---------------------------------------------------------------------*/
/*--------------- METHODS FOR PRINTING & SAVING THE MultiFlowDCRBlock --------------*/
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
  */

 void print( std::ostream & output , char vlvl = 0 ) const override;

/*--------------------------------------------------------------------------*/
/// extends Block::serialize( netCDF::NcGroup )
/** Extends Block::serialize( netCDF::NcGroup ) to the specific format of a
 * MultiFlowDCRBlock. See MultiFlowDCRBlock::deserialize(netCDF::NcGroup) for details of the
 * format of the created netCDF group. */

 void serialize( netCDF::NcGroup & file ) const override;

/*--------------------------------------------------------------------------*/

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

 c_Vec_FNumber & get_U( void ) const { return( CapTot ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 c_Vec_FNumber & get_LinkDelays( void ) const { return( LinkDelays ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 c_Vec_FNumber & get_NodeDelays( void ) const { return( NodeDelays ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 c_Vec_FNumber & get_MTU( void ) { return( MTU ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 c_Vec_FNumber & get_FlowBurst( void ) { return( FlowBursts ); }

/*--------------------------------------------------------------------------*/

 bool useFlowRelaxation( void ) const {
  return( ! ( AR & KnapsackRelaxation ) );
 }

/*--------------------------------------------------------------------------*/
 /// getting the current sense of the Objective, which is minimization

 int get_objective_sense( void ) const override final {
  return( f_sense );
  }

/*--------------------------------------------------------------------------*/

  double get_rs( Index k , Index i ) const {
    return( static_cast< SingleFlowDCRBlock * >( v_Block[ k ] )->get_r( i ) );
  }

/*--------------------------------------------------------------------------*/
 /// get the flow of a given arc for a given commodity 
 /** Given a commodity index k and an arc index ij, this function provides
  * the value of the associated variable x^k_ij. In the case of the knapsack
  * relaxation, the variables of the block are rescaled in such a way that
  * x \in [ 0 , 1 ]. In this case the functions get_flow provide the values
  * already rescaled wigth x^k_{ij} in [ 0 , u_ij ]. */

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

 ColVariable * get_flow_variable( Index k , Index i ) const {
  if( ! ( AR & HasVar ) )
   return( nullptr );

  if( ! ( AR & KnapsackRelaxation ) )
   return( static_cast< SingleFlowDCRBlock * >( v_Block[ k ] )->i2p_x( i ) );  

  return( nullptr );
  }

/*--------------------------------------------------------------------------*/

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
  * use the "void" constructor to call it. */

 //void CmnIntlz( void );

/* @} ----------------------------------------------------------------------*/
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

 static constexpr unsigned char addFixedCosts = 8; 

 static constexpr unsigned char slc = 8;
 ///< fourth bit of AR == 1: true if we use the strong forcing constraints

 Index NXtrV;          ///< Number of "extra" variables
 Index NXtrC;          ///< Number of "extra" constraints

 Subset IdxBeg;        ///< Description of "extra" constraints: start
 Subset CoefIdx;       ///< Description of "extra" constraints: indices
 Vec_CNumber CoefVal;  ///< Description of "extra" constraints: values

 Index NNodes;         ///< Number of nodes
 Index NArcs;          ///< Number of arcs
 Index NComm;          ///< Number of commodities
 Index NCnst;          ///< Number of arcs with mutual capacity constraints

 CMultiVector C;       ///< Matrix of the arc costs
 FMultiVector U;       ///< Matrix of the arc upper capacities
 FMultiVector B;       ///< Matrix of the node deficits
 FMultiVector I;       ///< Matrix of the variables integrality constraints

 Vec_FNumber UTot;       ///< Vector equal to the sum of mutual capacities
 Vec_FNumber CapTot;     ///< Vector of mutual capacities
 
 Vec_CNumber F;        ///< Vector of fixed costs
  
 Vec_FNumber NodeDelays;     ///< Vector of node delays
 Vec_FNumber LinkDelays;     ///< Vector of link delays

 Vec_FNumber FlowBursts;     ///< Vector of flow bursts
 Vec_FNumber FlowDeadlines;     ///< Vector of flow deadlines
 Vec_FNumber rho;     ///< Vector of rho
 Vec_FNumber MTU;     ///< Vector of MTU

 Subset Startn;        ///< Topology of the graph: starting nodes
 Subset Endn;          ///< Topology of the graph: ending nodes

 Index StrtNme;        ///< The "name" of the first node
 Subset NamesK;        /**< The dual multipliers relative to commodity K
			* start with NamesK[ k ] and end with NamesK[ k + 1 ]
			*/
 Subset Active;        /**< Set of the arcs for which a mutual capacity
			* constraint is defined */
 MultiSubset ActiveK;  ///< Like Active for individual capacities
 bool DrctdPrb;        ///< true if the problem is directed
 //std::vector<MCFType> PT;  ///< type of flow subproblem

 Vec_Bool CIsCpy;     ///< true for each row of C[] that is a copy of another
 Vec_Bool UIsCpy;     ///< true for each row of U[] that is a copy of another
 Vec_Bool BIsCpy;     ///< true for each row of B[] that is a copy of another

 std::vector< FRowConstraint > MCs;  ///< the static mutual capacity constrs
 
 int f_sense = Objective::eMin;

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

 void guts_of_destructor( void );
 
/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE FIELDS ------------------------------*/
/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;  // insert it in the Block factory

/*--------------------------------------------------------------------------*/

 };  // end( class( MultiFlowDCRBlock ) )

/*--------------------------------------------------------------------------*/

/*@}  end( group( MultiFlowDCRBlock_CLASSES ) ) ----------------------------*/
/*--------------------------------------------------------------------------*/

 }  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* MultiFlowDCRBlock.h included */

/*--------------------------------------------------------------------------*/
/*---------------------- End File MultiFlowDCRBlock.h ----------------------*/
/*--------------------------------------------------------------------------*/

