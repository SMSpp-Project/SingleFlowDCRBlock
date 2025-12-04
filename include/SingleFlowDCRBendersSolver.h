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

//#include "MCFClass.h"

//#include "MILPSolver.h"

#include "BenBound.h"

#include "DCR.h"

#include "BlockSolverConfig.h"

/*--------------------------------------------------------------------------*/
/*-------------------------- NAMESPACE & USING -----------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{
  
 //using namespace MCFClass_di_unipi_it;
 //using Index = Block::Index;
 
 //class SingleFlowDCRBendersSolverState;  // forward declaration of SingleFlowDCRBendersSolverState

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
/// Solver for SingleFlowDCRBlock

class SingleFlowDCRBendersSolver : public Solver , public BenBound {

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
 /** Void constructor: does nothing special, except verifying that the
  * template argument derives from BenBound. */

 SingleFlowDCRBendersSolver( void ) : Solver() , BenBound()  { }

/*--------------------------------------------------------------------------*/
 /// destructor: it has to release all the Modifications

 virtual ~SingleFlowDCRBendersSolver() { }

/** @} ---------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *
 * Parameter-wise, SingleFlowDCRBendersSolver maps the parameters of [CDA]Solver
**/

 /// set the (pointer to the) Block that the Solver has to solve

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
   flows.burst = MCFB->get_FlowBurst()[0];
   flows.rate = MCFB->get_rho()[0];        
   flows.deadline = MCFB->get_FlowDeadline()[0];   

   DCR::DCRLink links[narcs];
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

   DCR::DCRNode nodes[nnodes];
   for (int i = 0; i < nnodes; i++) {
      nodes[i] = {};
      nodes[i].delay = node_delay[i];
   }

   //DCR::DCRLink* link_ptr = links;
   //DCR::DCRNode* node_ptr = nodes;

   vector<double> mtu = MCFB->get_MTU();

   BenBound::LoadProblem(nnodes, narcs, flows, links, nodes, mtu[0]);
   // once done, read_unlock the SingleFlowDCRBlock (if it was read-lock()-ed)
   if( ! owned )
    MCFB->read_unlock();

   //std::cout << MCFB->get_U(0) << std::endl;
   //std::cout << "\n--------------------" << std::endl;

   //MCFB->print(std::cout);

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

 int compute( bool changedvars = true ) override
 {

  lock();  // first of all, acquire self-lock

  if( ! f_Block )           // there is no [SingleFlowDCRBlock] to solve
   return( kBlockLocked );  // return error 

  bool owned = f_Block->is_owned_by( f_id );       // check if already locked
  if( ( ! owned ) && ( ! f_Block->read_lock() ) )  // if not try to read_lock
   return( kBlockLocked );                         // return error on failure

  // while [read_]locked, process any outstanding Modification
/*
  if( ! f_dmx_file.empty() ) {  // if so required
   // output the current instance (after the changes) to a DMX file
   std::ofstream ProbFile( f_dmx_file , ios_base::out | ios_base::trunc );
   if( ! ProbFile.is_open() )
    throw( std::logic_error( "cannot open DMX file " + f_dmx_file ) );

   WriteMCF( ProbFile );
   ProbFile.close();
   }
*/

  process_outstanding_Modification();
  
  if( ! owned )             // if the [SingleFlowDCR]Block was actually read_locked
   f_Block->read_unlock();  // read_unlock it

  // ensure the timer exists (or reset it)
  BenBound::DCRsetTime( true );

  // then (try to) solve the SingleFlowDCR
  BenBound::DCRstartTime();
  BenBound::Solve();
  //std::cout << "SOLUTION: " << get_lb() << "\n";
  BenBound::DCRstopTime();
  //std::cout << "TIME: " << BenBound::getTime() << std::endl;

  unlock();                  // unlock the mutex     

  //if ( get_lb() <= 0.0 )
  //  std::cout << "ERROR!" << std::endl;

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

 double get_elapsed_time( void ) const override {
  return( this->BenBound::getTime() );
  }
 
/*--------------------------------------------------------------------------*/

 OFValue get_lb( void ) override {  
/*
  vector<double> soluzione;
  vector<double> costi;
  double sum = 0.0;

  auto MCFB = dynamic_cast< SingleFlowDCRBlock * >( f_Block );
  costi =  MCFB->get_C();
  
  soluzione.resize(MCFB->get_NArcs());
  for (int i = 0; i < MCFB->get_NArcs() ; i++){
     soluzione[i] = BenBound::getSolution( i );
     sum += costi[i] * soluzione[i];
     //if(soluzione[i] > 0)
      //std::cout << "soluzione[" << i << "] = " << soluzione[i] << " ";
  }
*/  
  //if( std::abs(sum-getObjVal()) > 2)
  //std::cout << get_ub() << "," << getLB() << std::endl;
  
  //return( sum ); 
  //std::cout << getObjVal() << std::end;
  
  //if(getLB() == -Inf<double>())
    //std::cout << "getLB() == -Inf<double>()" << std::endl;
    //return( getObjVal() );

  auto lb = this->BenBound::getLB();
  auto ub = this->BenBound::getUB();
  //std::cout<<lb<<","<<ub<<std::endl;
  //std::cout<<std::abs(ub-lb)/std::abs(ub)<<std::endl;

  if(std::abs(ub-lb)/std::abs(ub) > 1e-4){ //1e-2 
  ///if (lb <= 0.0){
    ///std::cout << lb << "," << ub << std::endl;
    ///std::cout << "Resolution with MILPSolver" << std::endl;
    auto bsc2 = dynamic_cast< BlockSolverConfig * >(Configuration::deserialize( "MILPPar1.txt" ) );
    auto SCFB = dynamic_cast< SingleFlowDCRBlock * >( f_Block );
    bsc2->apply( SCFB );
    ////Solver * slvr = SCFB->get_registered_solvers().front();
    ////int rtrn2 = slvr->compute();
    ////std::cout << "new_lb=" << slvr->get_lb() << std::endl;
    ////if(std::abs(slvr->get_lb()-ub)/ub > 1e-2) {std::cout<<std::abs(slvr->get_lb()-ub)/ub<<std::endl;}
    ////return(slvr->get_lb());  
    ////return(ub);
  }

  ////std::cout << "Resolution with BenderSolver" << std::endl;
  return( std::min(lb,ub) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 OFValue get_ub( void ) override { 

  //std::cout << BenBound::getUB() << " " << getLB() << std::endl;

  //return( getObjVal() );
/*
  if (this->BenBound::getLB() < 0.0){
    auto bsc2 = dynamic_cast< BlockSolverConfig * >(Configuration::deserialize( "MILPPar1.txt" ) );
    auto SCFB = dynamic_cast< SingleFlowDCRBlock * >( f_Block );
    bsc2->apply( SCFB );
    Solver * slvr = SCFB->get_registered_solvers().front();
    int rtrn2 = slvr->compute();
    return(slvr->get_ub());
  }
*/
  auto lb = this->BenBound::getLB();
  auto ub = this->BenBound::getUB();
  return( std::max(lb,ub) ); 

  vector<double> soluzione;
  vector<double> costi;
  double sum = 0.0;

  auto MCFB = dynamic_cast< SingleFlowDCRBlock * >( f_Block );
  costi =  MCFB->get_C();
  
  soluzione.resize(MCFB->get_NArcs());
  for (int i = 0; i < MCFB->get_NArcs() ; i++){
     soluzione[i] = BenBound::getSolution( i );
     sum += costi[i] * soluzione[i];
  }
  //std::cout << sum << "," << getObjVal() << std::endl;

  double DCR = 0.0;
  bool solution_is_feasible = true;
  auto min_r = Inf< double >();

  for (int i = 0; i < MCFB->get_NArcs() ; i++)
    if( min_r < soluzione[i] and soluzione[i] > 1e-6)
      min_r = soluzione[i];

  DCR += 1.0/min_r * MCFB->get_FlowBurst()[0];
  for (int i = 0; i < MCFB->get_NArcs() ; i++)
    if(soluzione[i] > 1e-6){
      DCR += 1.0/soluzione[i] * MCFB->get_MTU()[0] + MCFB->get_MTU()[ 0 ] / MCFB->get_U()[ i ] 
          + MCFB->get_LinkDelays()[ i ] + MCFB->get_NodeDelays()[ MCFB->get_SN( i ) - 1 ];
    }
  
  if(DCR > MCFB->get_FlowDeadline()[ 0 ]){
    solution_is_feasible = false;
    std::cout << "solution NOT feasible: " << (DCR - MCFB->get_FlowDeadline()[ 0 ])/DCR << std::endl;
  }

  if( solution_is_feasible )
    return( sum ); 
  else {
    std::cout << BenBound::getUB() << " " << get_lb() << std::endl;
    return( this->BenBound::getUB() ); 
  }  

}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

  OFValue get_var_value( void ) override { 

    return( getObjVal() );
    
    vector<double> soluzione;
    vector<double> costi;
    double sum = 0.0;

    auto MCFB = dynamic_cast< SingleFlowDCRBlock * >( f_Block );
    costi =  MCFB->get_C();
    
    soluzione.resize(MCFB->get_NArcs());
    for (int i = 0; i < MCFB->get_NArcs() ; i++){
      soluzione[i] = BenBound::getSolution( i );
      sum += costi[i] * soluzione[i];
    }
    
    return( sum ); 

  }

/*--------------------------------------------------------------------------*/

bool has_var_solution( void ) override {
  switch( this->BenBound::getStat() ) {
   case( BenBound::OK ):
   case( BenBound::Infeasible ): return( true );
   default:                      return( false );
   }
  }

/*--------------------------------------------------------------------------*/

 void get_var_solution( Configuration * solc = nullptr ) override
 {
  if( ! f_Block )  // no [SingleFlowDCR]Block to write to
   return;         // cowardly and silently return

  auto tsolc = dynamic_cast< SimpleConfiguration< int > * >( solc );
  if( tsolc && ( tsolc->f_value == 2 ) )
   return;
  
  auto DCRB = static_cast< SingleFlowDCRBlock * >( f_Block );
  int nnarc = DCRB->get_NArcs();

  //std::vector<double> v(nnarc, 0.0);
  //std::vector<double> vx(nnarc, 0.0);

  for( Index i = 0 ; i < nnarc ; ++i ){
    auto v = BenBound::getSolution( i );
    DCRB->set_r( i, v );
    if( v > 0.0 ){
      //std::cout << i << "," << v[i] << std::endl;
      DCRB->set_x( i, 1 );
    } else {
      DCRB->set_x( i, 0 );
    }
  }
}

/** @} ---------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/
/** @name Changing the data of the model
 *  @{ */

/*--------------------------------------------------------------------------*/
  
  void process_outstanding_Modification( void ) {

      bool reload = false;

      // note: since processing the Modification is fast, we don't bother with
      // being nice to other processes and do it all with v_mod under lock
       // try to acquire lock, spin on failure
      while( f_mod_lock.test_and_set( std::memory_order_acquire ) )
       ;
     
      // process all the Modifications
      for( auto mod : v_mod )
        if( auto tmod = dynamic_cast< C05FunctionModLinRngd * >( mod.get() ) ) {
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
        flows.burst = MCFB->get_FlowBurst()[0];
        flows.rate = MCFB->get_rho()[0];        
        flows.deadline = MCFB->get_FlowDeadline()[0];   

        DCR::DCRLink links[narcs];
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
            //std::cout << c[i] << std::endl;
        }

        vector<double> node_delay = MCFB->get_NodeDelays();

        DCR::DCRNode nodes[nnodes];
        for (int i = 0; i < nnodes; i++) {
            nodes[i] = {};
            nodes[i].delay = node_delay[i];
        }

        //DCR::DCRLink* link_ptr = links;
        //DCR::DCRNode* node_ptr = nodes;

        vector<double> mtu = MCFB->get_MTU();

        BenBound::LoadProblem(nnodes, narcs, flows, links, nodes, mtu[0]);
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

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

};  // end( class( SingleFlowDCRBendersSolver ) )
}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* SingleFlowDCRBendersSolver.h included */

/*--------------------------------------------------------------------------*/
/*---------------- End File SingleFlowDCRBendersSolver.h -------------------*/
/*--------------------------------------------------------------------------*/





