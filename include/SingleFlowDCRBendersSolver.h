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

#include "BenBound.h"

#include "DCR.h"

#include "BlockSolverConfig.h"

/*--------------------------------------------------------------------------*/
/*-------------------------- NAMESPACE & USING -----------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{
  
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

  process_outstanding_Modification();
  
  if( ! owned )             // if the [SingleFlowDCR]Block was actually read_locked
   f_Block->read_unlock();  // read_unlock it

  // ensure the timer exists (or reset it)
  BenBound::DCRsetTime( true );

  // then (try to) solve the SingleFlowDCR
  BenBound::DCRstartTime();
  BenBound::Solve();
  BenBound::DCRstopTime();

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

 double get_elapsed_time( void ) const override {
  return( this->BenBound::getTime() );
  }
 
/*--------------------------------------------------------------------------*/

 OFValue get_lb( void ) override {  
  return( this->BenBound::getLB() );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 OFValue get_ub( void ) override { 
  auto DCRB = static_cast< SingleFlowDCRBlock * >( f_Block );

  //for (int i = 0; i < DCRB->get_NArcs() ; i++)
    //std::cout << DCRB->get_r( i ) << std::endl;    

  //if( DCRB->is_feasible() )
    //std::cout << "feasible" << std::endl;

  if( this->get_var_value() < 1e200 && this->BenBound::getUB() > 1e200 )
   return( this->get_var_value() );
  
  //std::cout << this->get_var_value() << " " << BenBound::getUB() << std::endl;
  return( std::max( this->get_var_value() , this->BenBound::getUB() ) );
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

  OFValue get_var_value( void ) override { 

    ///!return( getObjVal() );
    
    double sum = 0.0;

    auto DCRB = dynamic_cast< SingleFlowDCRBlock * >( f_Block );
    auto costi =  DCRB->get_C();
    
    for (int i = 0; i < DCRB->get_NArcs() ; i++)
     sum += costi[i] * DCRB->get_r( i );
  
    return( sum ); 
  }

/*--------------------------------------------------------------------------*/

bool has_var_solution( void ) override {
  switch( this->BenBound::getStat() ) {
   case( BenBound::OK ): return( true ) ;
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
      //std::cout << i << "," << v << std::endl;
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





