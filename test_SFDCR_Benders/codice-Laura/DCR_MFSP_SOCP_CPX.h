/*--------------------------------------------------------------------------*/
/*--------------------------- File DCR_MFSP_SOCP_CPX.h----------------------*/
/*--------------------------------------------------------------------------*/

/** @file
 * Header file for the class DCR_MFSP_SOCP_CPX, which implements
 * a solver for Multi-Flow Single-Path (MFSP) DCR problems. We use two
 * Second-Order Cone formulations: "perspective" and "bigM" (the former
 * being stronger than the latter). We use CPLEX MIQCP solver.
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

#ifndef DCR_MFSP_SOCP_CPX_H
#define DCR_MFSP_SOCP_CPX_H

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "DCR.h"
#include <cplex.h>

/*--------------------------------------------------------------------------*/
/*--------------------------- CLASS DCR_MFSP_SOCP_CPX-----------------------*/
/*--------------------------------------------------------------------------*/

/** This class implements a solver for Multi-Flow Single-Path (MFSP)
 *  DCR problems using a Second Order Cone Programming (SOCP) 
 *  formulation and using CPLEX MIQCP solver. 
 *  Two SOCP formulations are implemented:
 *  - perspective
 *  - big-M
 *  The former is stronger than the latter.
 *  */

class DCR_MFSP_SOCP_CPX : public DCR
{

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/	

	public:
	
/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructors
    @{ */
    	
	DCR_MFSP_SOCP_CPX( void );
	
	/**< Constructor: gives some default values to all the data
   members of the class.*/
   
 /** @} */ 
	
/*--------------------------------------------------------------------------*/
/*----------------------------------READ DATA ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Loading the data of the problem
    @{ */
    	
	virtual void DCRloadProblem(int nnodes, int nlinks, int nflows, 
	DCRFlow *flows, DCRLink *links, DCRNode *nodes, double MTU, DCRDelay deltype);
	
 /** @} */		
/*--------------------------------------------------------------------------*/
/*----------------------------------WRITERS---------------------------------*/
/*--------------------------------------------------------------------------*/
	
	void DCRwrite( char *filename );
	
	/**< Write the model in LP format to a file named filename. 
	 *  This method does not appear in the base-class interface.
	 * This function can be called only *after* the problem has been loaded.
	 * \param filname name of the LP file*/
	
/*--------------------------------------------------------------------------*/
/*----------------------------------SOLVE-----------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Solving the problem
    @{ */
    
    //FIXME: consider more status codes ???	
	virtual DCRStatus DCRsolve( void );
	
	/**< Solves the problem using a SOCP formulation using CPLEX MIQCP solver.
	 * This function can be called only *after* the problem has been loaded.
	 *  */

/** @} */ 
	
/*--------------------------------------------------------------------------*/
/*----------------------------------GET RESULTS-----------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading solver output 
  @{ */
  	
	virtual double DCRgetObj( void );
	
	/**< Returns the optimal value.
	 *  This function can be called only *after* the problem has been solved.
	  */
	
	virtual int DCRgetNodes( void );	
	
	/**< Returns the number of nodes explored in the Branch-and-Bound tree. 
	 * This method does not appear in the base-class interface.
	 * This function can be called only *after* the problem has been solved.
	  */ 
	
	virtual int DCRgetUBSolNPaths( int k );
	
	/**< Returns 1 if a (mixed-)integer solution exists, since we only
	 *  solve *single-path* instances (see UBsol). 
	 * This function can be called only *after* the problem has been solved.
	 * \param k flow index
	  */ 
	
	virtual int DCRgetUBSolPath( int k, int p, int *X, double *R);
	
	/**< Returns the number of hops in path 0 of flow k in the current (mixed-)integer solution, 
	 * if any (see UBsol). It only works for p=0, since we only solve *single-path* instances. 
	 * X is an array of int of size (n-1) allocated by the user to store
	 * the indices of the links in path p. R is an array of double
	 * of size (n-1) allocated by the used to store the corresponding rates.
	 * This function can be called only *after* the problem has been solved. 
	 * \param k flow index
	 * \param p path index (only 0)
	 * \param X pointer to an array of int of size (numNodes-1)
	 * \param R pointer to an array of double of size (numNodes-1) */
	 
	 virtual int DCRgetPSolNPaths( int k );
	 
	 /**< Returns the 1 if a *continuous* solution exists, since we only
	  *  solve *single-path* instances (see Psol). 
	  * Note that in the case of a continuous solution the concept of path no longer exists,
	  * because the flow can be splitted. Yet, since we also consider multi-path versions of the problem
	  * it makes sense to have different paths also for the continuous case. 
	  * This function can be called only *after* the problem has been solved.
	  * \param k flow index*/
	  
	  virtual void  DCRgetPSolPath(int k, int p, double *X, double *R);
	  
	  /**< X is an array of double of size numLinks allocated by the user to store
	   * the values of the path (x_ij) variables in the current *continuous* solution,
	   * if any (see Psol). It only works for p=0, since we only solve *single-path* instances.
	   * R is an array of double of size numLinks allocated by the user to store
	   * the values of the rates (r_ij) variables in the current *continuous* solution.
	   * This function can be called only *after* the problem has been solved.
	   * \param k flow index
	   * \param p path index (only 0)
	   * \param X pointer to an array of double of size numLinks 
	   * \param R pointer to an array of double of size numLinks */
	
/** @} */ 
			
/*--------------------------------------------------------------------------*/
/*----------------------------------SETTERS---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
    @{ */
    
    //FIXME use relative tolerance
	virtual void DCRsetOptEps( double OE = 0 );
	/**< Overrides DCR base-class method to set CPLEX optimality tolerance.
	 *  This function can be called only *after* the problem has been loaded.  
     * \param OE double expressing tolerance to set*/
    
    //FIXME use relative tolerance
	virtual void DCRsetFsbEps( double FE = 0);
	/**< Overrides DCR base-class method to set CPLEX feasibility tolerance.  
	 * This function can be called only *after* the problem has been loaded.
     * \param FE double expressing tolerance to set*/
     
     //FIXME: the base class takes an "ostream" though...
	virtual void DCRsetLog(char *filename);
	/**< Overrides DCR base-class method to set CPLEX log file. 
	 * This function can be called only *after* the problem has been loaded. 
     * \param filename a string containing the name of the log file*/
     
     virtual void DCRsetTimeLimit( double );
	/**< Overrides DCR base-class method to set time-limit. 
	 * This function can be called only *after* the problem has been loaded.
	 * */
     
	virtual void DCRsetCuts( int );
	/**< This method does not appear in the base-class interface. 
	 * This function can be called only *after* the problem has been loaded.
	 * */
	
	virtual void DCRsetNodeLimit( int );
	/**< This method does not appear in the base-class interface. 
	 * This function can be called only *after* the problem has been loaded.
	 * */
	
	virtual void DCRsetHeur( int );
	/**< This method does not appear in the base-class interface. 
	 * This function can be called only *after* the problem has been loaded.
	 * */
	
	virtual void DCRsetQCPType( void );
	/**< Makes SOCP model continuous (i.e., x_ij path variables in [0,1]). 
	 * This method does not appear in the base-class interface. 
	 * This function can be called only *after* the problem has been loaded.
	 * */
	
	virtual void DCRsetMIQCPType( void );
	/**< Makes SOCP model mixed-integer (i.e., x_ij path variables in {0,1}). 
	 * This method does not appear in the base-class interface.
	 * This function can be called only *after* the problem has been loaded.
	 * */
	
	virtual void DCRsetNodeFile( int );	
	/**< This method does not appear in the base-class interface. 
	 * This function can be called only *after* the problem has been loaded.
	 * */

	virtual void DCRsetMIQCPStrat(int );
	/**< This method does not appear in the base-class interface. 
	 * This function can be called only *after* the problem has been loaded.
	 * */
	
	virtual void DCRsetNodeSelection( int );
	/**< This method does not appear in the base-class interface. 
	 * This function can be called only *after* the problem has been loaded.
	 * */	
	
	virtual void DCRsetModel(char model);
	 /**< Sets the type of SOCP model to use: either perspective ('p') or bigM ('b'). 
	  * This method does not appear in the base-class interface.
	  * This function can be called only *after* the problem has been loaded.
     * \param model char value can be 'p' for perspective model or 'b' for bigM model*/
	
	/** @} */ 	
	
/*--------------------------------------------------------------------------*/
/*----------------------------------MODIFIERS-------------------------------*/
/*--------------------------------------------------------------------------*/

/** @name Changing the data of the problem
    @{ */

virtual void DCRcloseArcs( int * whch , int na);
/**< Closes arcs of the network, for all flows.
 * This function can be called only *after* the problem has been loaded.
 * \param whch indices of the arcs to be closed
 * \param na number of arcs to be closed (for each flow)
 * */

virtual void DCRcloseArcs( int k , int * whch, int na);
/**< Closes arcs of the network for a specific flow.
 * This function can be called only *after* the problem has been loaded.
 * \param which indices of the arcs to be closed
 * \param k index of the flow to which the process applies
 * \param na number of arcs to be closed for flow k
 * */

/** @} */ 
	
/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Destructor
    @{ */
    		
	~DCR_MFSP_SOCP_CPX();

    /**< Frees up dinamically allocated memory */

/** @} */ 
	
/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/	
/** @name Private part of the class
    @{ */
    	
	private:
	
	//network data
	DCRNode *Nodes;
	DCRLink *Links;
	DCRFlow *Flows;
	
	//CPLEX data
	CPXENVptr env ;
	CPXLPptr  lp ;	
	
	//SOCP formulation flag (0=perspective, 1=bigM)
	int flagBigM;
	
	//variable offsets
	int r_offset;
	int rmin_offset;
	int s_offset;
	int t_offset;	
	int rp_offset;
	int theta_offset;
	int FOFFSET;	
	
	//solution info
	double *Xsol;
	double *Rsol;	
	double objval;

	//private methods
	//create SOCP model methods
	void addCols( int flow ); 
	void addFlowRows( int flow ); 
	void addRateBounds( int flow ); 
	void addDelayRows( int flow );
	void addDelayRowsBigM( int flow ); 
	void addCapacityRows( void );
	int getNumInLinks(int node);
	int getNumOutLinks(int node);
	bool isLinkInNode(int node, int arc);
	bool isLinkOutNode(int node, int arc);
	bool isFlowSource(int node, int flow);
	bool isFlowSink(int node, int flow);
	//memory management methods
	void copyDataArrays(DCRFlow *flows, DCRLink *links, DCRNode *nodes);
	void clean_up( void );
	//solution management methods
	void getSol();
	
	/** @} */ 
	
};



#endif

/*--------------------------------------------------------------------------*/
/*-------------------------- End File DCR_MFSP_SOCP_CPX.h ------------------*/
/*--------------------------------------------------------------------------*/
