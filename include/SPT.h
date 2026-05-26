/*--------------------------------------------------------------------------*/
/*--------------------------- File SPT.h ----------------------------*/
/*--------------------------------------------------------------------------*/

/*<SPT per grafi diretti con singola sorgente*/

/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef SPT_H
#define SPT_H

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/
#include<vector>

using namespace std;
/*--------------------------------------------------------------------------*/
/*--------------------------- CLASS LineSearch -----------------------------*/
/*--------------------------------------------------------------------------*/

class SPT
{


/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/	

	public:

/*--------------------------------------------------------------------------*/
/*--------------------------- Inf() and Eps() ------------------------------*/
/*--------------------------------------------------------------------------*/
/** Very small class to simplify extracting the "+ infinity" value for a
    basic type; just use Inf<type>(). */

 template <typename T>
  class Inf {
   public:
  Inf() {}
  operator T() { return( std::numeric_limits<T>::max() ); }
  };

/*--------------------------------------------------------------------------*/
/** Very small class to simplify extracting the "machine epsilon" for a
    basic type; just use Eps<type>(). */

 template <typename T>
  class Eps {
   public:
  Eps() {}
  operator T() { return( std::numeric_limits<T>::epsilon() ); }
  };

	
/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/

   enum Status 
   { 
	  OK, 	///<  solver found a feasible solution
	  Infeasible,   ///<  problem infeasible
	  Unbounded,    ///<  problem unbounded
	  Error
   };

/*--------------------------------------------------------------*/

   struct SPTLink 
   {
     int startnode;
     int endnode;
     double cost;

     double rstar; //non molto utile in SP, ma evita un doppio controllo 
                   //in ingresso e uscita per la soluzione del lagrangiano.
   };
	
/*archi semplificati per SP, si richiede solo il costo dell'arco, oltre a testa e coda*/

/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/

	SPT( void );

/*--------------------------------------------------------------------------*/
/*---------------------------------- METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

	virtual void Solve();
	/*<risolve SP con l'algoritmo di Bellman*/

	virtual int getNHops();
	/*<restituisce il numero di archi del cammino ottimo*/

	virtual Status getStatus();
	/*<restituisce lo status del problema, ovvero se è stato risolto o è Unbounded*/

    virtual vector<int> getPath();
	/*<restituisce il cammino ottimo dalla sorgente al pozzo, 
	 * sotto forma di indici dei link coinvolti.
	 * NB. se viene usato su un grafo ridotto questi non coincideranno con quelli
	 * del grafo completo*/

	virtual vector<double> getLabel();
	/*<restituisce la soluzione duale di SP.*/

	virtual void LoadProblem(int nnodes, int nlinks, vector<SPTLink> links, int sourcenode, int sinknode);

	virtual void updCosts(vector<SPTLink> links);

/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

	~SPT();

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/	


	private:

	vector<SPTLink> Links;
	int numNodes;
	int numLinks;
	int nhops;
	int s; //sourceindex.
	int t; //sinkindex.
	vector<int> Sol;  //indici del cammino s-t soluzione di SP.
	vector<double> DualSol; //soluzione duale di SP (etichette).
	Status stat; //per controllare la presenza di cicli di costo negativo.


	void copyDataArrays(vector<SPTLink> links);
	int getLink(int from, int to);
	/*<trova, se ve ne è uno, l'arco con la testa = to e la coda = from*/

    void clean_up(void);


};

#endif

/*--------------------------------------------------------------------------*/
/*-------------------------- End File SPT.h ----------------------------*/
/*--------------------------------------------------------------------------*/
