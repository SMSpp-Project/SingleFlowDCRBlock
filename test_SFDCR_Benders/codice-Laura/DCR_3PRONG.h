/*--------------------------------------------------------------------------*/
/*--------------------------- File DCR_3PRONG.h-------------------------------*/
/*--------------------------------------------------------------------------*/

/** @file
 * Header file for the class DCR_3PRONG to solve Single-Flow Single-Path (SFSP)
 * Delay Constrained Routing (DCR) problems using the 3-pronged heuristic
 * described in Frangioni et al. (2013).
 * 
 * \version 1.00
 *
 * \date July - 2013
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

#ifndef DCR_3PRONG_H
#define DCR_3PRONG_H

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "DCR.h"
#include "DCR_SPT.h"
#include "DCR_MFSP_SOCP_CPX.h"

/*--------------------------------------------------------------------------*/
/*--------------------------- CLASS DCR_3PROG --------------------------------*/
/*--------------------------------------------------------------------------*/

/** This class implements the 3-pronged heuristic algorithm
 *  described in Frangioni et al. (2013).
 *  This algorithm works as follows:
 * 	i.)   run ERA-I algorithm, if infeasible RETURN infeasible.
 *  ii.)  run ERA-H algorithm, if feasible RETURN solution.
 *  iii.) run SOCP algorithm, RETURN solution. 
 *  
 *  */
 
class DCR_3PRONG : public DCR
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
    	
	DCR_3PRONG( void );
	
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
/*----------------------------------SOLVE-----------------------------------*/
/*--------------------------------------------------------------------------*/	
	/** @name Solving the problem
    @{ */
    	
	virtual DCRStatus DCRsolve( void );
	
	/**< Solves the problem using two heuristics. */

/** @} */ 
	
/*--------------------------------------------------------------------------*/
/*----------------------------------GET RESULTS-----------------------------*/
/*--------------------------------------------------------------------------*/	
/** @name Reading solver output 
  @{ */
  	
	virtual double DCRgetObj( void );
	/**< Returns best upper-bound value found (may not be optimal in general,
	 * because this is an heuristic  algorithm) */
	
	virtual int DCRgetUBSolNPaths( int k );
	
	/**< Returns 1 if a (mixed-)integer solution exists, since we only
	 *  solve single-flow *single-path* instances (see UBsol). 
	 * \param k flow index
	  */ 
	
	virtual int DCRgetUBSolPath( int k, int p, int *X, double *R);
	
	/**< Returns the number of hops in path 0 of flow k in the current (mixed-)integer solution, 
	 * if any (see UBsol). It only works for p=0, since we only solve single-flow *single-path* instances. 
	 * X is an array of int of size (n-1) allocated by the user to store
	 * the indices of the links in path p. R is an array of double
	 * of size (n-1) allocated by the used to store the corresponding rates. 
	 * \param k flow index
	 * \param p path index (only 0)
	 * \param X pointer to an array of int of size (numNodes-1)
	 * \param R pointer to an array of double of size (numNodes-1) */
	 
	 virtual int DCRgetPSolNPaths( int k );
	 
	 /**< This function is not supported because heuristics
	  * cannot solve the continuous relaxation of the problem.*/
	  
	  virtual void  DCRgetPSolPath(int k, int p, double *X, double *R);
	  
	  /**< This function is not supported because heuristics
	  * cannot solve the continuous relaxation of the problem.*/
	
/** @} */ 

	
/*--------------------------------------------------------------------------*/
/*----------------------------------MODIFIERS-------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Changing the data of the problem
    @{ */

virtual void DCRcloseArcs( int * whch, int na );
/**< Closes arcs of the network, for all flows.
 * \param whch indices of the arcs to be closed*/

virtual void DCRcloseArcs( int k , int * whch, int na );
/**< Closes arcs of the network for a specific flow.
 * \param which indices of the arcs to be closed
 * \param k index of the flow to which the process applies*/

/** @} */ 	

/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Destructor
    @{ */
    		
	~DCR_3PRONG();
	
 /**< Frees up dinamically allocated memory */
/** @} */ 	
	
/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/	
/** @name Private part of the class
    @{ */
    	
	private:
	
	//network data
	DCRNode *Nodes;
	DCRLink *Links;
	DCRFlow *Flows;
	
	//solution info
	int *Xsol;
	double *Rsol;	
	double objval;
	int nhops;	
	DCR_MFSP_SOCP_CPX* socp;//FIXME: and GUROBI ???
	DCR_SPT* eraI;
	DCR_SPT* eraH;
	
	
	//private methods
	//memory management methods
	void copyDataArrays(DCRFlow *flows, DCRLink *links, DCRNode *nodes);
	void clean_up( void );
	int DCRgetLink(int from, int to);

};

#endif

/*--------------------------------------------------------------------------*/
/*-------------------------- End File DCR_3PRONG.h ----------------------------*/
/*--------------------------------------------------------------------------*/
