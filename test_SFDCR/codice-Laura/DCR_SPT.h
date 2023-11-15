/*--------------------------------------------------------------------------*/
/*--------------------------- File DCR_SPT.h--------------------------------*/
/*--------------------------------------------------------------------------*/

/** @file
 * Header file for the class DCR_SPT, which implements
 * two heuristics for Single-Flow Single-Path (SFSP) DCR problems. 
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

#ifndef DCR_SPT_H
#define DCR_SPT_H

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "DCR.h"

/*--------------------------------------------------------------------------*/
/*--------------------------- CLASS DCR_SPT --------------------------------*/
/*--------------------------------------------------------------------------*/

/** This class implements two heuristics for Single-Flow Single-Path (SFSP)
 *  DCR problems.
 *  */
 
class DCR_SPT : public DCR
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
    	
	DCR_SPT( void );
	
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
		
	virtual int DCRgetUBSolNPaths( int k );
	
	/**< Returns 1 if a (mixed-)integer solution exists, since we only
	 *  solve *single-path* instances (see UBsol). 
	 * \param k flow index
	  */ 
	
	virtual int DCRgetUBSolPath( int k, int p, int *X, double *R);
	
	/**< Returns the number of hops in path 0 of flow k in the current (mixed-)integer solution, 
	 * if any (see UBsol). It only works for p=0, since we only solve *single-path* instances. 
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
/*----------------------------------SETTERS---------------------------------*/
/*--------------------------------------------------------------------------*/
	/** @name Other initializations
    @{ */
		 
	virtual void DCRsetHeur(char);
	 /**< Sets the type of heuristic to use.
	  * This method does not appear in the base-class interface.*/
	  
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
    		
	~DCR_SPT();
	
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
	
	//can be either 1=ERA-I or 2=ERA-H
	int heur;
		
	//solution info
	int *Xsol;
	double *Rsol;	
	double objval;
	int nhops;
	
	//private methods
	int DCRgetLink(int from, int to);
	void DCRheurERAI();
	void DCRheurERAH();
	//memory management methods
	void copyDataArrays(DCRFlow *flows, DCRLink *links, DCRNode *nodes);
	void clean_up( void );
	//solution management methods
	
	/** @} */ 

};

#endif

/*--------------------------------------------------------------------------*/
/*-------------------------- End File DCR_SPT.h ----------------------------*/
/*--------------------------------------------------------------------------*/
