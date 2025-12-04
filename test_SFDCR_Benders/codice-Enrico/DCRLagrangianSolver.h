/*--------------------------------------------------------------------------*/
/*--------------------------- File DCRLagrangianSolver.h ----------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef DCR_LAGRANGIAN_SOLVER_H
#define DCR_LAGRANGIAN_SOLVER_H

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "DCR.h"
#include "SPT.h"

#include <vector>

using namespace std;

/*--------------------------------------------------------------------------*/
/*--------------------------- CLASS DCRLagrangianSolver -----------------------------*/
/*--------------------------------------------------------------------------*/

class DCRLagrangianSolver
{
 
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/	
	public:


/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/

  enum LAGStatus 
   { 
    OK,   ///<  solver found a feasible solution
    Infeasible,   ///<  problem infeasible
    Unbounded,    ///<  problem unbounded
    Error
   };

/*--------------------------------------------------------------------------*/

  struct LinearCut
  {
     	double m; //coefficiente angolare del taglio
	    double q; //intercetta del taglio
  };

/*--------------------------------------------------------------------------*/  

  struct Cut_Val  //struttura usata per riottimizzazione
  {
      LinearCut c; //taglio
      vector<double> RSol;  //soluzione associata
      int RSolsize; //sua dimensione

      double rmin; //rmin al momento della generazione del taglio, utile per aggiustarne la pendenza
  };

/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/

   DCRLagrangianSolver(void);

/*--------------------------------------------------------------------------*/
/*---------------------------------- METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

  
   void Solve();

   LAGStatus getStatus();
   /*<restituisce lo status del problema*/

   double getLambda();
   /*<restituisce il punto di intersezione tra i vincoli*/

   double getInterVal();
   /*<restituisce il valore dell'intersezione tra i vincoli nel punto lambda (serve?)*/

   double getObjVal();
   /*<restituisce il valore della funzione nel punto lambda*/

   vector<int> getPath();
   /*<restituisce il cammino ottimo, da chiamare solo per lambda ottimo nullo*/

   int getnumHops();
   /*<restituisce il numero di archi della soluzione ottima*/

   double getHeurVal();
   /*<restituisce il valore della f.o. per soluzioni ammissibili*/

   double getOptBeta();
   /*<pendenza in concomitanza dell'ottimo*/

   vector<double> getRSol();
   /*<restituisce il vettore degli r_ij della soluzione ottima*/

   vector<double> getRSolCosts();
   /*<restituisce il vettore degli f_ij relativi agli r_ij della soluzione ottima*/  

   vector<int> getRedGraphPos(); 
   /*<restituisce il grafo ridotto, come sue posizioni nel grafo totale, per il controllo delle condizioni di Bellman*/

   int getCardRedGraph();
   /*<restituisce la cardinalità del grafo ridotto*/

   vector<int> getComplGraph();
   /*<restituisce il complementare del grafo ridotto, come sue posizioni nel grafo totale, per il controllo sulla non convessità*/

   vector<double> getCheckSol();
   /*<restituisce la combinazione tra le sol pos e neg per il check di cplex*/

   vector<int> getCheckSolPos();
   /*<e le sue posizioni nel grafo*/

   vector<int> getCheckSolNeg();
   /*<e le sue posizioni nel grafo*/

   vector<double> getSPLabels();
   /*<restituisce i potenziali di SP per la soluzione ottima*/

   int getNumIte();
   /*<restituisce il numero di iterazioni del lagrangiano per un flusso*/

   //vector<Cut_Val> getCuts();
   /*<restituisce l'insieme dei tagli attualmente considerati (serve?)*/

   double DCRgetTime( void );

   bool isFeasible ();
   /*<inizializza senza risolvere, restituisce 1 se il primale è vuoto*/

   bool is0opt();
   /*<ci dice se lambda = 0 è ottimo usando al pià due SPT, da chiamare post caricamento dati, utile per reopt in Benders*/

   //FIXME:dovrei richiedere anche la precisione desiderata? (implementato con precisione di macchina)
   void LoadProblem (int nnodes, int nlinks, DCR::DCRFlow flow, DCR::DCRLink *links, DCR::DCRNode *nodes, double mtu, double r_min);
   /*<trasmette i dati del problema a spt, da invocare prima di Solve() e dei Get Method, 
    * in caso contrario daranno errori, o valori di deafault.*/

   void LoadProblem(DCR::DCRFlow flow);
   /*<altro tipo di caricamento dati, riapre gli archi modificati (vedere il suo FIXME) e permette di modificare il
    *flusso, da usare solo dopo aver usato il carimento di tutti i dati..*/

   void updrmin(double r_min);
    /*<permette di modificare il valore di r_min senza ricaricare gli alri dati*/

   void closeArcs(vector<int> arcs, int na);
  /*<permette di chiudere gli archi passati per indice: la loro capacità viene posta a 0*/

   void openArcs(int * arcs, int na);
  /*<permette di riaprire gli archi passati per indice: la loro capacità viene 
   *riportata al suo valore iniziale */

   void DCRsetTimeLimit(long secs);

   void DCRsetTime(bool timeON);
 
   void DCRstartTime();

   void DCRstopTime();
/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

	~DCRLagrangianSolver();

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/	
	
	private:
    
    void copyDataArray(DCR::DCRFlow flow, DCR::DCRLink *links, DCR::DCRNode *nodes);

    void getLimitVal();
    /*<risolve SP in un caso pessimo e restituisce un valore limite*/

    double getCost(double r, int lindex);
    /*<restituisce il costo per SP a seconda del tipo di arco.*/

    void Inizial();
    /*<operazioni di inizializzazione, modifica lo status
     * se la pendenza in 0 è negativa, o se il primale è vuoto.*/

    void UpdCut(double alpha, double beta);
    /*<aggiunge il nuovo taglio ottenuto dal punto e
     * aggiorna il nuovo taglio ottimo.*/

    void setReducedGraph();
    /*<considera solo gli archi con capacità non maggiore di r_min*/ 

    void setSPTcosts();
    /*<traduce i dati in input per SPT Solve, che richiede solo il vettore
       degli archi e i costi.*/

    int Reopt();
    /*<tenta di saltare la fase di inizializzazione riscaldando la 
     * line search, sfruttando le informazioni generate alle iterazioni precedenti*/


/*--------------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/	
	
	protected:

	LinearCut pCut; //uno dei tagli che definisce la sol ottima.
	LinearCut mCut; //l'altro taglio che la definisce.
  vector<Cut_Val> Cuts; //insieme dei tagli attualmente in gioco. 
  double InterVal; //valore dell'intersezione.
  double lambda;    //punto di intersezione.
  double ObjVal;   //valore della funzione in lambda.
  double LimitVal; //valore pessimo per DCR.
  double HeurVal; //valore in un punto in cui ci sia ammissibilità per il vincolo sul ritardo.
  double optBeta; //pendenza del lagr in concomitanza dell'ottimo.

  int num_ite; //numero di iterazioni del lagrangiano 

  //variabili per salvare i dati in ingresso e per trasmetterli a SPT.
  double rmin;
  double MTU;
  DCR::DCRNode * Nodes;
  DCR::DCRLink * Links; 
  DCR::DCRFlow  Flow; //solo uno può essere accettato, dare errore altrimenti.
  int numNodes,numLinks;
  int cardRedGraph; //numero di archi del grafo ridotto.
  vector<double> ModCaps; //vettore delle capacità modificate su cui lavoriamo
  //int * closedArcs; //vettore degli archi che sono stati chiusi
  //int nclar; //commentati perché potrebbe servire a riaprirne, ma non è davvero necessario

  //variabili per i dati in uscita da SPT.
  int nhops; //numero di archi del cammino ottimo.
  vector<DCR::DCRLink> RedGraLinks; //archi DCR del grafo ridotto. 
  vector<int> RedGraPos; //posizioni del grafo ridotto in quello totale.
  vector<int> ComplGraph; //complementare del grafo ridotto.
  vector<int> OptPath; //posizioni del cammino ottimo.
  int solvedflag;
  int inizialflag;
  int reoptflag;
  //double * XSol;
  int numHops;
  vector<double> RSol;
  vector<double> RSolCosts; //f_ij associati agli r*_ij.
  vector<double> SPLabels; //potenziale del pozzo della soluzione di SP.

  int maxCutSize; //massimo numero di tagli che salveremo per la riottimizzazione.

  ///VARIABILI PER I CHECH SU CPLEX utilizzabili per migliorare l'euristica
  vector<double> solpos;
  vector<double> solneg;
  double betaneg;
  double betapos;
  int nonposflag;
  vector<int> checkSolPos;
  vector<int> checkSolNeg;

  double eps; //precisione di macchina
  SPT spt;     //oggetto spt per la soluzione della LineSearch.
  vector<SPT::SPTLink> Linksp; //archi semplificati per SP.
  SPT::Status spstat; //status che dichiara quale esito ha avuto l'SP.
  LAGStatus lagstat;  //status che dichiara quale esito ha avuto in lagrangiano.

  OPTtimers *timer;  ///< timer
  long tlimit; ///< time limit

  void clean_up();

  

};

#endif

/*--------------------------------------------------------------------------*/
/*-------------------------- End File DCRLagrangianSolver.h ---------------------------*/
/*--------------------------------------------------------------------------*/