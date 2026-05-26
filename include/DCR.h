/*--------------------------------------------------------------------------*/
/*--------------------------- File DCR.h -----------------------------------*/
/*--------------------------------------------------------------------------*/

/** @file
 * Header file for the abstract (pure virtual) base class DCR, which
 * defines a standard interface for solvers of Delay Constrained Routing
 * (DCR) problems.
 * 
 * \version 1.00
 *
 * \date April - 2013
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Laura Galli \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy 2013 by Antonio Frangioni, 
 * 	   Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 */ 

/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef DCR_H
#define DCR_H


/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "OPTUtils.h"

using namespace std;	
//using namespace OPTtypes_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*--------------------------- CLASS DCR ------------------------------------*/
/*--------------------------------------------------------------------------*/

/** This class defines a standard abstract interface for solvers of
 *  Delay Constrained Routing (DCR) problems. A DCR problem consists of
 *  an (Integer) Multi-Commodity Flow part plus a max-delay constraint part.
 *  The Multi-Commodity Flow part can be single/multi flow and single/multi path.
 *  The max-delay part is defined via network calculus and different delay formula
 *  can be used according to different network traffic shapers (e.g., 
 *  Strictly-Rate-Proportional, Weakly-Rate-Proportional, Frame-Based).  
 *  */

class DCR 
{
	
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/	
	public:
	
/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/

/** DCR solver status*/
	enum DCRStatus 
	{ 
	  OK = 0, 	///<  solver found a feasible solution
	  Stopped,      ///<  solver stopped 
	  Infeasible,   ///<  problem infeasible
	  Unbounded,    ///<  problem unbounded
	  Error         ///<  solver error
	};
	
/** DCR delay formula*/	
	enum DCRDelay
	{
	    SRP = 0, ///< Strictly Rate Proportional
	    WRT,     ///< Weakly Rate Proportional
	    FB	     ///< Frame Based
	};

/** DCR network flow */	
	struct DCRFlow
	{
		int sourcenode;     ///< flow source
		int sinknode;       ///< flow sink
		double burst;       ///< flow burst
		double rate;        ///< flow rate
		double deadline;    ///< flow deadline
		double *costs;      ///< arc-flow costs
		double *caps;       ///< arc-flow individual capacity
	};
	
/** DCR network link */
	struct DCRLink
	{
		int startnode;		///< link start-node
		int endnode;		///< link end-node
		double speed;       ///< link speed
		double capacity;	///< link mutual capacity 
		double delay;       ///< link delay
		double cost;        ///< link cost
	};
	
/** DCR network node */
	struct DCRNode
	{
		double delay; 		///< node delay
	};
	
/** Very small class for DCR exceptions
 */
	class DCRException : public std::exception 
	{
		public:
		DCRException( const char *const msg = 0 ) { errmsg = msg; }

		const char* what( void ) const throw () { return( errmsg ); }

		private:
		const char *errmsg;
	};
	

/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructors
    @{ */		

	DCR( void )
	
	{	
		numNodes = numLinks = numFlows = 0;		
		MTU = 0;
		dtype = SRP; //i.e., SRP=strictly rate proportional delay
		
		optEps = fsbEps = 0;  
		
		log = 0; 
		verbosity = 0; 
	
		timer = 0;  
		
		tlimit = 0;	
		
		Psol = false;
		UBsol = false;
	}
	/**< Constructor of the base class: gives some default values to the data
   structure of the base class, that can be changed later in the constructors
   of the derived classes.*/
	
/** @} */ 

/*--------------------------------------------------------------------------*/
/*--------------------------- SET METHODS-----------------------------------*/
/*--------------------------------------------------------------------------*/

/** @name Other initializations
    @{ */
    
	virtual void DCRsetVerbosity( int lvl = 0 )
	{
			verbosity = lvl;
	}
    
    /**< Sets the level of verbosity (the higher the more verbose). 
     * \param lvl int expressing the level of verbosity */
    
/*--------------------------------------------------------------------------*/
  
	//FIXME: would be better to use an *ostream* type to make it more general  
    virtual void DCRsetLog( ostream *outs ) 
    {
	log = outs;
    }  
    
    /**< Creates a log outstream. 
     * \param outs ostream pointer*/ 

/*--------------------------------------------------------------------------*/

	virtual void DCRsetTime( bool timeON = true )
	{
		if( timeON )
			if( timer ) timer->ReSet();
			else timer = new OPTtimers();
		else
			delete timer; 
    }
    
    /**< If timeON is true sets or resets the timer, if false deletes the timer.
     * \param timeON bool value */

/*--------------------------------------------------------------------------*/

	virtual void DCRstartTime( void )
	{
		if( timer == 0 )
			throw( DCR::DCRException( "DCR::DCRstartTimer(): failed to start timer" ) ); 
		
		timer->Start();
    }
    
    /**< If timer was set or re-set, starts or re-starts timer ticking. */

    
/*--------------------------------------------------------------------------*/

	virtual void DCRstopTime( void )
	{
		if ( timer == 0 )
			throw( DCR::DCRException( "DCR::DCRstopTimer(): failed to stop timer" ) );

		timer->Stop();
	}
    
    /**< If timer was set, stops timer ticking. */

/*--------------------------------------------------------------------------*/

	virtual void DCRsetTimeLimit( long secs )
	{
		tlimit = secs;
	}	
	
	 /**< Sets a time limit for the solver.
	  * \param secs long expressing timelimit in seconds */	

/*--------------------------------------------------------------------------*/

	virtual void DCRsetOptEps( double OE = 0 )
	{
		optEps = OE;
    }
    
    /**< In many cases, only an "approximate" solution of the problem is possible;
   alternatively, only an "approximate" solution may be required for the
   purposes of the caller (in order to save time).
   The exact meaning of "approximate" is solver-dependent, but the more
   common ways in which this happens are

   - either the value of the solution is not exactly optimal;

   - or the constraints are not exactly satisfied.

   SetOptEps() tells that any solution that is OE-optimal w.r.t. the value of
   the objective function can be considered optimal. 
   * \param OE double expressing optimality tolerance 
    */

/*--------------------------------------------------------------------------*/

	virtual void DCRsetFsbEps( double FE = 0)
	{
		fsbEps = FE;
    }
    
    /**< In many cases, only an "approximate" solution of the problem is possible;
   alternatively, only an "approximate" solution may be required for the
   purposes of the caller (in order to save time).
   The exact meaning of "approximate" is solver-dependent, but the more
   common ways in which this happens are

   - either the value of the solution is not exactly optimal;

   - or the constraints are not exactly satisfied.

   SetFsbEps() tells that any solution where the violation of the 
   constraints is not larger than FE can be considered feasible.
   * \param FE double expressing feasibility tolerance 
*/
   
/** @} */ 
    
/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR SOLVING THE PROBLEM ---------------------*/
/*--------------------------------------------------------------------------*/
/** @name Solving the problem
    @{ */
    
	
	virtual DCRStatus DCRsolve( void ) = 0;

/**< DCR solver, returns a DCRStatus value */
	
/** @} */ 
	
/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading solver output 
  @{ */
  
	virtual double DCRgetObj( void ) = 0;
	
	/**< Returns objective function value in the current solution, if any
	 * (see Psol and UBsol).*/
	
/*--------------------------------------------------------------------------*/
	
	virtual int DCRgetUBSolNPaths( int k ) = 0;
	
	/**< Returns the number of paths of flow k in the current (mixed-)integer solution, 
	 * if any (see UBsol).
	 * \param k flow index
	  */ 
	
	virtual int DCRgetUBSolPath( int k, int p, int *X, double *R) = 0;
	
	/**< Returns the number of hops in path p of flow k in the current (mixed-)integer solution, 
	 * if any (see UBsol). 
	 * X is an array of int of size (n-1) allocated by the user to store
	 * the indices of the links in path p. R is an array of double
	 * of size (n-1) allocated by the used to store the corresponding rates. 
	 * \param k flow index
	 * \param p path index
	 * \param X pointer to an array of int of size (numNodes-1)
	 * \param R pointer to an array of double of size (numNodes-1) */
	 
	 virtual int DCRgetPSolNPaths( int k ) = 0;
	 
	 /**< Returns the number of (possibly "splitted") "paths" of flow k in the current *continuous* solution, 
	  * if any (see Psol).
	  * Note that in the case of a continuous solution the concept of path no longer exists,
	  * because the flow can be splitted. Yet, since we also consider multi-path versions of the problem
	  * it makes sense to have different paths also for the continuous case. 
	  * \param k flow index*/
	  
	  virtual void  DCRgetPSolPath(int k, int p, double *X, double *R) = 0;
	  
	  /**< X is an array of double of size numLinks allocated by the user to store
	   * the values of the path (x_ij) variables in the current *continuous* solution,
	   * if any (see Psol).
	   * R is an array of double of size numLinks allocated by the user to store
	   * the values of the rates (r_ij) variables in the current *continuous* solution.
	   * \param k flow index
	   * \param p path index
	   * \param X pointer to an array of double of size numLinks 
	   * \param R pointer to an array of double of size numLinks */
	   
	    virtual void getNewSol()
	    {
		Psol = false;
		UBsol = false;    
	    }
	    /**< Moves to the next solution if multiple solutions are given.
	     *   If no more solutions are available, sets all solution flags to zero. */
	
/*--------------------------------------------------------------------------*/

	virtual double DCRgetTime( void )
	{		
		return( timer ? timer->Read() : 0 );
	}
	
	/**< Returns the elapsed time. If the clock is ticking, returns the *total*
    time since the last Start() without stopping the clock; otherwise,
    returns the total elapsed time of all the past runs of the clock since
    the last reset. */
    
/*--------------------------------------------------------------------------*/


 	virtual void DCRgetTime( double &t_us , double &t_ss )
	{
		t_us = t_ss = 0;
		if( timer ) timer->Read( t_us , t_ss ); 
	}
	
/**< As DCRgetTime ( void ) but *adds* user and system time to tu and ts. 
 * \param t_us double expressing user time in seconds
 * \param t_ss double expressing system time in seconds*/
	
/** @} */ 

/*--------------------------------------------------------------------------*/
/*----------------------------------LOAD METHOD ----------------------------*/
/*--------------------------------------------------------------------------*/

/** @name Loading the data of the problem
    @{ */
    		
	virtual void  DCRloadProblem(int nnodes, int nlinks, int nflows, 
	DCRFlow *flows, DCRLink *links, DCRNode *nodes, double MTU, DCRDelay deltype) = 0;
	
	/**< Reads data from memory and creates a DCR problem. 
	 * \param nnodes int expressing the number of nodes in the network
	 * \param nlinks int expressing the number of links in the network
	 * \param nflows int expressing the number of flows in the network
	 * \param flows pointer to an array of DCRFlow 
	 * \param links pointer to an array of DCRLink
	 * \param nodes pointer to an array of DCRNode
	 * \param MTU double expressing maximum transmit unit
	 * \param deltype DCRdelay enum expressing type of delay formula to be used*/ 
	
/** @} */ 
	
/*--------------------------------------------------------------------------*/
/*----------------------------------GET METHODS-----------------------------*/
/*--------------------------------------------------------------------------*/

/** @name Reading the data of the problem
    @{ */
		
	virtual int DCRgetNumNodes( void ) const
	{
		return( numNodes );
	}
	
	/**< Returns number of nodes in the network. */
/*--------------------------------------------------------------------------*/

	virtual int DCRgetNumLinks( void ) const
	{
		return( numLinks );
	}
	
	/**< Returns number of links in the network. */

/*--------------------------------------------------------------------------*/

	virtual int DCRgetNumFlows( void ) const
	{
		return( numFlows );
	}
	
	/**< Returns number of flows in the network. */

/*--------------------------------------------------------------------------*/
	
	virtual double DCRgetMTU( void ) const
	{
		return( MTU );
	}
	
	/**< Returns Maximum Transmit Unit (MTU) of the network. */

/** @} */ 

/*--------------------------------------------------------------------------*/
/*----------------------------------MODIFIERS-------------------------------*/
/*--------------------------------------------------------------------------*/

//TODO: change costs, demand, capacities, open and close arcs.

/** @name Changing the data of the problem
    @{ */

virtual void DCRcloseArcs( int * whch , int na) = 0;

virtual void DCRcloseArcs( int k , int * whch , int na) = 0;

/** @} */ 

/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

/** @name Destructor
    @{ */
    
	virtual ~DCR( void )
	{
		delete timer; 
	}
	
    /**< Frees up dinamically allocated memory */

/** @} */ 

/*--------------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/	

/** @name Protected fields of the class
    @{ */	
    
	protected:
	
	int numNodes; ///< number of nodes in the network
	int numLinks; ///< number of links in the network
	int numFlows; ///< number of flows in the network
	
	double MTU; ///< maximum transmit unit
	
	double optEps;  ///< optimality tolerance
	double fsbEps;  ///< feasibility tolerance
	
	ostream *log; ///< log output
	int verbosity; ///< verbosity level
	
	OPTtimers *timer;  ///< timer
	long tlimit; ///< time limit
	
	DCRDelay dtype; ///< network delay formula used
	
	bool Psol; ///< true if continuous solution was found, false otherwise
	bool UBsol; ///< true if (mixed)-integer solution was found, false otherwise
};

/** @} */ 

#endif
	
/*--------------------------------------------------------------------------*/
/*-------------------------- End File DCR.h --------------------------------*/
/*--------------------------------------------------------------------------*/
	
