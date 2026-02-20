/*--------------------------------------------------------------------------*/
/*--------------------------- File BenBound.h ----------------------------*/
/*--------------------------------------------------------------------------*/


/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef BENBOUND_H
#define BENBOUND_H

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "DCRLagrangianSolver.h" //serve?
#include "DCR.h"

#include <vector>


using namespace std;
//using namespace OPTtypes_di_unipi_it;
/*--------------------------------------------------------------------------*/
/*--------------------------- CLASS Benders -----------------------------*/
/*--------------------------------------------------------------------------*/

class BenBound
{


/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/	

	public:
	
/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/

   enum BndrStat
   { 
	  OK, 	///<  solver found a feasible solution
	  Infeasible,   ///<  problem infeasible
	  Unbounded,    ///<  problem unbounded
	  Error
   };

/*--------------------------------------------------------------------------*/


   struct SubInterval //è necessario salvare i due tagli ottimi di volta in volta?
   {
     double rmin;
     double inter; //ascissa dell'intersezione tra le rette, candidata a nuovo valore di rmin nel sottointervallo..
     double interVal; //lower bound dato dall'intersezione delle rette.
     double Val; //upper bound per ObjVal ottenuto dal sottointervallo.

     int bestCutpos;// segnala la posizione del taglio migliore per i casi in cui infeasflag sia 1.

     int branchedflag; //ci dice se non c'è più speranza di trovare l'ottimo in questo sottointervallo
     int solflag; //ci dice se per \bar{r}_{min}=inter avevamo già risolto il lagrangiano.
     int infeasflag; //ci dice se nel primo intervallo ci sono solo tagli positivi

     DCRLagrangianSolver::LinearCut pCut;
     DCRLagrangianSolver::LinearCut mCut; //coppia di tagli che definisce la soluzione ottima.
     vector<DCRLagrangianSolver::LinearCut> Cuts; //per salvare i tagli attualmente in gioco nel sottointervallo.
   };

/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/

	BenBound(void);

/*--------------------------------------------------------------------------*/
/*---------------------------------- METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/


	void LoadProblem(int nnodes, int nlinks, DCR::DCRFlow flow, DCR::DCRLink *links, DCR::DCRNode *nodes, double mtu);
    //FIXME: aggiungere un metodo pubblico per modificare solo il flusso
    // ed uno per modificare myparam.
    void LoadProblem(DCR::DCRFlow flow);

	void Solve();

    double getObjVal();
    /*<restituisce il valore ottimo trovato*/

    double getApproxVal();
    /*<restituisce l'approssimazione per questo valore ottimo*/

    double getr_min();
    /*<r_min ottimo*/

    double getHeurVal();
    /*<valore della migliore soluzione primale ammissibile trovata finora*/

    int getNumIterationBender();
    /*<numero di punti visitati durante la soluzione del Bender*/

    int getNumIterationLagr();
    /*<somma totale del numero di iterazioni del lagrangiano*/

    bool IsConvexInteration();
    /*<ritorna true se c'è stata almeno un'iterazione convessa*/

    vector<double> getyPlotData();

    vector<double> getxPlotData();

    double getTime() const;

    BndrStat getStat();

    void DCRsetTimeLimit(long secs);

    void DCRsetTime(bool timeON);
 
    void DCRstartTime();

    void DCRstopTime();

    double getSolution(int i);
    /*<soluzione primale ammissibile trovata finora*/

    double getUB();

    double getLB();

    int get_Links(){ return(numLinks); };

    double* get_caps(){ return(caps);}

/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

	~BenBound();

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/	


	private:
    
    DCR::DCRFlow Flow;
    DCR::DCRNode* Nodes;
    DCR::DCRLink* Links;
    int numNodes;
    int numLinks;
    double MTU;
    double eps;
    double ObjVal;
    int solvedflag;
    double rmin; //rmin in concomitanza della soluzione ottima.
    double HeurVal; //valore dell'euristica
    double ApproxVal; //valore dell'approssimazione.
    double mystep;
    double BestUB;
    double BestLB;
    double crit_capc;
    int counter_ite_Ben; //counter per numero di iterazioni di Bender.
    int counter_ite_Lag; //counter per numero totale di iterazioni del Lagrangiano.
    bool is_convex_iteration; //flag per capire se si è trovato almeno un taglio convesso

    double SOL_VALUE;
    vector<double> SOLUTION;

    OPTtimers *timer;  ///< timer
    long tlimit; ///< time limit

    //dati per la stampa della funzione
    vector<double> xPlot;
    int numPoints;

    double* caps; //vettore delle capacità degli archi.

    BndrStat BenStat;

    vector<double> SPLabels; //duali di SP
    DCRLagrangianSolver lagSol;
    vector<SubInterval> Q; //sottointervalli di V considerati.

    double myparam; //paramentro in (0,1) che sceglie il punto nei casi di non ammissibilità
    
    void Inizial();
    double* Limitrmin();
    void capSort(int sx, int dx);
    int Distrib(int sx, int pivot, int dx);
    void Swap(int a, int b);
    

    int LineSearch();//restituisce la posizione in Q dell'intersezione ottima tra i tagli
    int whchbest(int i); // metodo di supporto alla LS.
    void UpdCut(double alpha, double beta, int i); // metodo di supporto alla LS.
    

    void copyDataArray(DCR::DCRFlow flow, DCR::DCRLink *links, DCR::DCRNode *nodes);
    void clean_up();




};

#endif

/*--------------------------------------------------------------------------*/
/*-------------------------- End File BenBound.h ----------------------------*/
/*--------------------------------------------------------------------------*/