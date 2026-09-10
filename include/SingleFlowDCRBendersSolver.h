/*--------------------------------------------------------------------------*/
/*------------------ File SingleFlowDCRBendersSolver.h ---------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the SingleFlowDCRBendersSolver class, implementing a
 * Solver for Delay-Constrained Routing problems (DCR) relative to a Single
 * Flow, as set by SingleFlowDCRBlock, via a "Benders with nested Lagrange"
 * approach.
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
 * \copyright &copy; by Antonio Frangioni
 */ 
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __SingleFlowDCRBendersSolver
 #define __SingleFlowDCRBendersSolver
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "Solver.h"

#include "SingleFlowDCRBlock.h"

#include "BenBound.h"

#include "DCR.h"

#include "BlockSolverConfig.h"

#include "CDASolver.h"

/*--------------------------------------------------------------------------*/
/*-------------------------- NAMESPACE & USING -----------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{
/*--------------------------------------------------------------------------*/
/*------------------------------- CLASSES ----------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup SingleFlowDCRBendersSolver_CLASSES Classes in SingleFlowDCRBendersSolver.h
 *  @{ */

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS SingleFlowDCRBendersSolver --------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// Solver for SingleFlowDCRBlock via a "Benders with nested Lagrange" scheme
/** SingleFlowDCRBendersSolver implements the Solver concept [see Solver.h]
 * for SingleFlowDCRBlock by wrapping the (non-SMS++, "legacy") BenBound
 * class, which it publicly inherits from together with Solver.
 *
 * The DCR problem is a min-cost routing problem in which, besides the
 * flow-conservation and capacity constraints, a worst-case end-to-end
 * delay bound (computed via network calculus, as a function of the
 * reserved rate r[ i , j ] on each arc and of the "burst" of the flow)
 * must not exceed a given deadline. BenBound tackles the resulting
 * nonlinear (albeit convex, arc-wise separable) problem via a Benders-like
 * decomposition on r_min, the minimum reserved rate along the used path:
 * for a fixed value of r_min, the residual problem decouples into a
 * (Lagrangian-relaxable) shortest-path-like problem, solved by the nested
 * DCRLagrangianSolver (which in turn solves a sequence of shortest path
 * problems via SPT); the resulting piecewise-linear lower bound as a
 * function of r_min is then explored by a line-search (Benders master)
 * over r_min itself. See BenBound.h and DCRLagrangianSolver.h for the
 * details of this "legacy" solution algorithm.
 *
 * This class merely translates the data of a SingleFlowDCRBlock into the
 * format expected by BenBound::LoadProblem() (see set_Block() and
 * process_outstanding_Modification()), drives the solution process (see
 * compute()), and translates the results back into SMS++ terms (see
 * get_lb(), get_ub(), get_var_solution(), ...). */

class SingleFlowDCRBendersSolver : public Solver , public BenBound
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public Types
 *  @{ */

/** @} ---------------------------------------------------------------------*/
/*----------------- CONSTRUCTING AND DESTRUCTING SingleFlowDCRBendersSolver */
/*--------------------------------------------------------------------------*/
/** @name Constructing and destructing SingleFlowDCRBendersSolver
 *  @{ */

 /// constructor: does nothing special
 /** Void constructor: does nothing special besides default-constructing
  * both base classes, Solver and BenBound. */

 SingleFlowDCRBendersSolver( void ) : Solver() , BenBound()  { }

/*--------------------------------------------------------------------------*/
 /// destructor: does nothing special
 /** Destructor: currently empty; the Solver and BenBound base classes take
  * care of releasing their own resources (comprised any pending
  * Modification) via their own destructors. */

 virtual ~SingleFlowDCRBendersSolver() { }

/** @} ---------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *
 * SingleFlowDCRBendersSolver currently does not define any Solver
 * parameter of its own (no {get,set}_{num_,}*_par() method is overridden):
 * it only inherits the (empty) set of parameters of the base Solver class.
 *  @{ */

 /// set the (pointer to the) Block that the Solver has to solve
 /** Sets (or resets, if \p block == nullptr) the SingleFlowDCRBlock that
  * this Solver has to solve (an exception is thrown if \p block is not a
  * SingleFlowDCRBlock). Besides the usual bookkeeping (via
  * Solver::set_Block()), this [read_]locks the Block (unless already
  * owned), extracts source/sink node (the nodes with, respectively,
  * negative and positive deficit), traffic parameters (burst, rate,
  * deadline), MTU, and the per-arc topology/speed/capacity/delay/cost
  * data, translates them into the DCR::DCRFlow / DCR::DCRLink / DCR::DCRNode
  * format expected by BenBound, and calls BenBound::LoadProblem() with
  * them; the same translation is repeated (see
  * process_outstanding_Modification()) whenever the SingleFlowDCRBlock is
  * subsequently modified. If the call does nothing (\p block == f_Block),
  * it cowardly and silently returns. */

 void set_Block( Block * block ) override
 {

  if( block == f_Block )  // actually doing nothing
   return;                // cowardly and silently return

  Solver::set_Block( block );  // attach to the new Block

  if( block ) {  // this is not just resetting everything
   auto MCFB = dynamic_cast< SingleFlowDCRBlock * >( f_Block );
   if( ! MCFB )
    throw( std::invalid_argument(
		         "SingleFlowDCRBendersSolver:set_Block: block must be a SingleFlowDCRBlock" ) );

   bool owned = MCFB->is_owned_by( f_id );
   if( ( ! owned ) && ( ! MCFB->read_lock() ) )
    throw( std::logic_error( "cannot acquire read_lock on SingleFlowDCRBlock" ) );
   // load the new SingleFlowDCRBlock into the :BenBound object

   vector<double> B = MCFB->get_B();
   int source, sink;

   int nnodes = MCFB->get_NNodes();
   int narcs = MCFB->get_NArcs();

   for(int i = 0; i < nnodes; i++) {
    if( B[i] < 0 )
      source = i;
    if( B[i] > 0 )
      sink = i;
   }

   DCR::DCRFlow flows;
   flows = {};
   flows.sourcenode = source;
   flows.sinknode = sink;
   flows.burst = MCFB->get_FlowBurst();
   flows.rate = MCFB->get_rho();        
   flows.deadline = MCFB->get_FlowDeadline();   

   vector<double> u = MCFB->get_U(); 
   vector<double> c = MCFB->get_C();
   vector<Index> sn = MCFB->get_SN();
   vector<Index> en = MCFB->get_EN();
   vector<double> link_delay = MCFB->get_LinkDelays();
   vector<DCR::DCRLink> links(narcs);

   for (int i = 0; i < narcs; i++) {
      links[i] = {};
      links[i].startnode = sn[i]-1;
      links[i].endnode = en[i]-1;
      links[i].speed = u[i];
      links[i].capacity = u[i];
      links[i].delay = link_delay[i];
      links[i].cost = c[i];
   }

   vector<double> node_delay = MCFB->get_NodeDelays();

   vector<DCR::DCRNode> nodes(nnodes);
   for (int i = 0; i < nnodes; i++) {
      nodes[i] = {};
      nodes[i].delay = node_delay[i];
   }

   //DCR::DCRLink* link_ptr = links;
   //DCR::DCRNode* node_ptr = nodes;

   double mtu = MCFB->get_MTU();

   BenBound::LoadProblem(nnodes, narcs, flows, links, nodes, mtu);

   // once done, read_unlock the SingleFlowDCRBlock (if it was read-lock()-ed)
   if( ! owned )
    MCFB->read_unlock();

   // TODO: maybe log it
   }
  }  // end( set_Block )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

/** @} ---------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE Block ----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Solving the SingleFlowDCR 
 *  @{ */

 /// (try to) solve the SingleFlowDCR
 /** (Tries to) solve the SingleFlowDCRBlock currently attached to this
  * Solver via the "Benders with nested Lagrange" scheme implemented by
  * BenBound (see the class GENERAL NOTES). After acquiring the self-lock
  * and, if necessary, a read_lock on f_Block, processes any outstanding
  * Modification (see process_outstanding_Modification()), then calls
  * BenBound::DCRsetTime() / BenBound::Solve() / BenBound::DCRstopTime() to
  * actually run the algorithm. \p changedvars is currently unused.
  *
  * IMPLEMENTATION NOTE: the method contains a (large) commented-out block
  * that would run an external warm-start Solver (via BlockSolverConfig,
  * reading "MILPPar-initial.txt", "DCRCfgPC.txt" and a
  * "WarmStartCfg.txt"-selected Solver) to compute an initial bound on
  * r_min (stored in bound_PC) before calling BenBound; this warm-start is
  * currently disabled.
  *
  * Returns kOK unconditionally (i.e., irrespective of the actual status
  * reported by BenBound::getStat(), which can be queried via
  * has_var_solution()). */

 int compute( bool changedvars = true ) override
 {

  lock();  // first of all, acquire self-lock

  if( ! f_Block )           // there is no [SingleFlowDCRBlock] to solve
   return( kBlockLocked );  // return error 

  bool owned = f_Block->is_owned_by( f_id );       // check if already locked
  if( ( ! owned ) && ( ! f_Block->read_lock() ) )  // if not try to read_lock
   return( kBlockLocked );                         // return error on failure

  // while [read_]locked, process any outstanding Modification
  
  if( ! owned )             // if the [SingleFlowDCR]Block was actually read_locked
   f_Block->read_unlock();  // read_unlock it

  process_outstanding_Modification();

 auto f_Block_aux = f_Block;
/*
 const std::string warmstart_cfg_file = "MILPPar-initial.txt";

 auto warmstart_cfg = Configuration::deserialize( warmstart_cfg_file );
 auto warmstart_bsc = dynamic_cast< BlockSolverConfig * >( warmstart_cfg );

 auto dcrb_pc = dynamic_cast< BlockConfig * >
          ( Configuration::deserialize( "DCRCfgPC.txt" ) );

 // instantiate the warm-start Solver via the factory, picking name and
 // (optional) ComputeConfig from WarmStartCfg.txt. The Solver is NOT
 // registered on f_Block so it cannot interfere with PrimalProximalHeur
 // itself (which is already attached to f_Block).
 auto warmstart =
  dynamic_cast< CDASolver * >(
   Solver::new_Solver( warmstart_bsc->get_SolverName( 0 ) ) );
 
 if( warmstart_bsc->num_ComputeConfig() > 0 )
    if( auto cc = warmstart_bsc->get_SolverConfig( 0 ) )
      warmstart->set_ComputeConfig( cc );
 dcrb_pc->apply( f_Block );
 warmstart->set_Block( f_Block );
 warmstart->compute( changedvars );
 warmstart->get_var_solution();
 bound_PC = dynamic_cast< SingleFlowDCRBlock * >( f_Block )->get_rmin();

 std::cout << "bound_PC = " << bound_PC << std::endl;
 std::cout << "sol = " << warmstart->get_lb() << std::endl;
*/
  // ensure the timer exists (or reset it)
  BenBound::DCRsetTime( true );
  // then (try to) solve the SingleFlowDCR
  BenBound::Solve();
  BenBound::DCRstopTime();

  //std::cout << "LB = " << get_lb() << std::endl;
  //std::cout << "UB = " << get_ub() << std::endl;

  unlock();                  // unlock the mutex     

  // now give out the result: note that the vector MCFstatus_2_sol_type[]
  // starts from 0 whereas the first value of MCFStatus is -1 (= kUnSolved),
  // hence the returned status has to be shifted by + 1
  
  return( kOK );
  }

/** @} ---------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Accessing the found solutions (if any)
 *  @{ */

 /// returns the time (in seconds) spent so far in BenBound::Solve()

 double get_elapsed_time( void ) const override {
  return( this->BenBound::getTime() );
  }

/*--------------------------------------------------------------------------*/
 /// returns the best lower bound found so far by the Benders/Lagrange bound

 OFValue get_lb( void ) override {

  return( this->BenBound::getLB() );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns the best upper bound: BenBound's UB if DCR-feasible, else the
 /// (relaxed) heuristic value
 /** If the current best solution found by BenBound is (approximately)
  * feasible for the DCR delay constraint (see is_DCR_feasible()), returns
  * BenBound::getUB(); otherwise, since that solution does not actually
  * respect the delay deadline, the (possibly optimistic) relaxed value
  * BenBound::getHeurVal() is returned instead. */

 OFValue get_ub( void ) override {

  if( is_DCR_feasible() )
    return( this->BenBound::getUB() );
  else
    return( this->BenBound::getHeurVal() );
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns the objective value of the current reserved-rate solution
 /** Recomputes, from scratch, the objective value \f$ \sum_i C[ i ] r[ i ]
  * \f$ of the current reserved-rate solution stored in the
  * SingleFlowDCRBlock, rather than using BenBound::getObjVal() (which is
  * currently not used, see the commented-out line below). */

  OFValue get_var_value( void ) override {

    //return( this->BenBound::getObjVal() );

    double sum = 0.0;

    auto DCRB = dynamic_cast< SingleFlowDCRBlock * >( f_Block );
    auto costi =  DCRB->get_C();

    for (int i = 0; i < DCRB->get_NArcs() ; i++)
     sum += costi[i] * DCRB->get_r( i );

    return( sum );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns the number of Benders iterations performed so far

 Index get_BenIt( void ) {
  return( this->BenBound::getNumIterationBender() );
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns the total number of (nested) Lagrangian iterations performed

 Index get_LagIt( void ) {
  return( this->BenBound::getNumIterationLagr() );
}

/*--------------------------------------------------------------------------*/
 /// returns true if BenBound found a feasible (primal) solution
 /** Returns true if and only if BenBound::getStat() == BenBound::OK, i.e.,
  * a feasible solution was found and can be retrieved by
  * get_var_solution(). */

bool has_var_solution( void ) override {
  switch( this->BenBound::getStat() ) {
   case( BenBound::OK ): return( true ) ;
   default:                      return( false );
   }
  }

/*--------------------------------------------------------------------------*/
 /// writes the BenBound solution back into the SingleFlowDCRBlock
 /** Copies the per-arc reserved-rate solution found by BenBound (see
  * BenBound::getSolution()) into the r[] Variable of the SingleFlowDCRBlock
  * f_Block, deriving the corresponding "flow" ColVariable x[ i ] (1 if the
  * arc is used, i.e., its reserved rate is > 0, 0 otherwise) and the
  * theta[ i ] Variable (MTU / r[ i ] if the arc is used, 0 otherwise); \p
  * solc is currently unused. If f_Block is nullptr the method silently
  * does nothing. */

 void get_var_solution( Configuration * solc = nullptr ) override
 {
  if( ! f_Block )  // no [SingleFlowDCR]Block to write to
   return;         // cowardly and silently return

  //auto tsolc = dynamic_cast< SimpleConfiguration< int > * >( solc );
  //if( tsolc && ( tsolc->f_value == 2 ) )
  // return;
  
  auto DCRB = static_cast< SingleFlowDCRBlock * >( f_Block );
  int nnarc = DCRB->get_NArcs();

  double rmin = Inf<double>();
  double theta_min = Inf<double>();

  for( Index i = 0 ; i < nnarc ; ++i ){
    auto v = BenBound::getSolution( i );
    DCRB->set_r( i , v );
    if( v > 0.0 ){
      DCRB->set_x( i , 1 );
      rmin = std::min( rmin , v );
      DCRB->set_theta( i , DCRB->get_MTU() / v );
    } else {
      DCRB->set_x( i , 0 );
      DCRB->set_theta( i , 0 );
    }
  }
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns true if the current BenBound solution meets the DCR deadline
 /** Checks (with relative tolerance 1e-2) whether the current solution
  * found by BenBound satisfies the worst-case delay bound of the DCR
  * problem, computed via network calculus: for every arc j whose reserved
  * rate BenBound::getSolution( j ) is at least the flow's sustained rate
  * DCRB->get_rho() (i.e., every arc actually used on the routing path),
  * the per-arc worst-case delay contribution MTU / r[ j ] + ( MTU / U[ j ]
  * + LinkDelays[ j ] + NodeDelays[ start node of j ] ) is accumulated;
  * adding the burst term FlowBurst / r_min (r_min being the smallest
  * reserved rate along the path, see BenBound::getr_min()) yields the
  * overall worst-case end-to-end delay, which is compared against the
  * flow's deadline DCRB->get_FlowDeadline(). */

 bool is_DCR_feasible( void ) {

  auto DCRB = dynamic_cast< SingleFlowDCRBlock * >( f_Block );
  int nnarc = DCRB->get_NArcs();
  auto r_min = BenBound::getr_min();
  auto theta_min = DCRB->get_FlowBurst() / r_min;

  vector<double> LinkDelays = DCRB->get_LinkDelays();
  vector<double> NodeDelays = DCRB->get_NodeDelays();
  double lhs = 0.0;

  for( Index j = 0 ; j < nnarc ; ++j ) {
    auto v = BenBound::getSolution( j );
    if( v >= DCRB->get_rho() )
      lhs += DCRB->get_MTU() / v + ( DCRB->get_MTU() /  DCRB->get_U( j ) +  LinkDelays[ j ] + NodeDelays[ DCRB->get_SN( j ) - 1 ] );
  }

  lhs += theta_min;
  bool DCR_feasible = false;

  //std::cout << " pviol=" << ( lhs - DCRB->get_FlowDeadline() ) / DCRB->get_FlowDeadline() << std::endl;

  if( ( lhs - DCRB->get_FlowDeadline() ) / DCRB->get_FlowDeadline() <= 1e-2 ){
    DCR_feasible = true;
  }
  return( DCR_feasible );
}

/** @} ---------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/
/** @name Changing the data of the model
 *  @{ */

/*--------------------------------------------------------------------------*/
 /// processes any Modification issued by the SingleFlowDCRBlock since the
 /// last compute()
 /** Drains, under lock, the list v_mod of Modification collected by the
  * Solver since the last call (see Solver::process_outstanding_Modification
  * for the general mechanism); if at least one Modification is found, the
  * whole BenBound problem is rebuilt from scratch out of the current state
  * of the SingleFlowDCRBlock, exactly as done in set_Block() (topology,
  * costs, capacities, delays, flow burst/rate/deadline and MTU are all
  * re-read and passed to BenBound::LoadProblem()). No attempt is made to
  * incrementally update the BenBound data structures according to the
  * specific Modification received; any change at all triggers a full
  * reload. */

  void process_outstanding_Modification( void ) {

      bool reload = false;

      // note: since processing the Modification is fast, we don't bother with
      // being nice to other processes and do it all with v_mod under lock
       // try to acquire lock, spin on failure
      while( f_mod_lock.test_and_set( std::memory_order_acquire ) )
       ;
     
      // process all the Modifications
      for( auto mod : v_mod )
        if( auto tmod = mod.get() ) {
          reload = true;  // a reset must be done
          break;          // ignore all the remaining Modifications
        }
     
      v_mod.clear();  // all Modifications tackled, clear the list
     
      f_mod_lock.clear( std::memory_order_release );  // release lock
     
      if( reload ){

        auto MCFB = dynamic_cast< SingleFlowDCRBlock * >( f_Block );
        vector<double> B = MCFB->get_B();
        int source, sink;

        int nnodes = MCFB->get_NNodes();
        int narcs = MCFB->get_NArcs();

        for (int i = 0; i < nnodes; i++) {
          if (B[i] < 0 )
            source = i;
          if (B[i] > 0 )
            sink = i;
        }

        DCR::DCRFlow flows;
        flows = {};
        flows.sourcenode = source;
        flows.sinknode = sink;
        flows.burst = MCFB->get_FlowBurst();
        flows.rate = MCFB->get_rho();        
        flows.deadline = MCFB->get_FlowDeadline();   

        vector<DCR::DCRLink> links(narcs);
        vector<double> u = MCFB->get_U();
        vector<double> c = MCFB->get_C();
        vector<Index> sn = MCFB->get_SN();
        vector<Index> en = MCFB->get_EN();
        vector<double> link_delay = MCFB->get_LinkDelays();

        for (int i = 0; i < narcs; i++) {
            links[i] = {};
            links[i].startnode = sn[i]-1;
            links[i].endnode = en[i]-1;
            links[i].speed = u[i];
            links[i].capacity = u[i];
            links[i].delay = link_delay[i];
            links[i].cost = c[i];
        }

        vector<double> node_delay = MCFB->get_NodeDelays();

        vector<DCR::DCRNode> nodes(nnodes);
        for (int i = 0; i < nnodes; i++) {
            nodes[i] = {};
            nodes[i].delay = node_delay[i];
        }

        //DCR::DCRLink* link_ptr = links;
        //DCR::DCRNode* node_ptr = nodes;

        double mtu = MCFB->get_MTU();

        BenBound::LoadProblem(nnodes, narcs, flows, links, nodes, mtu);
      }
  }

/** @} ---------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

protected:

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*---------------------------- PROTECTED FIELDS  ---------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

private:

double bound_PC = 0.0;  ///< bound on r_min computed by an (optional) warm-start
                         ///< primal heuristic Solver; currently only set by
                         ///< the code commented out inside compute(), so it
                         ///< is unused in the present implementation

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

};  // end( class( SingleFlowDCRBendersSolver ) )

/*--------------------------------------------------------------------------*/

/** @}  end( group( SingleFlowDCRBendersSolver_CLASSES ) ) -----------------*/
/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* SingleFlowDCRBendersSolver.h included */

/*--------------------------------------------------------------------------*/
/*---------------- End File SingleFlowDCRBendersSolver.h -------------------*/
/*--------------------------------------------------------------------------*/





