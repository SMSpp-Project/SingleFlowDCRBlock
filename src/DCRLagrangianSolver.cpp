/*--------------------------------------------------------------------------*/
/*-------------------------- File DCRLagrangianSolver.cpp ------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------------- IMPLEMENTATION -------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "DCR.h"
#include "SPT.h"
#include "DCRLagrangianSolver.h"

#include <string.h>
#include <fstream>
#include <vector>
#include <cmath>

using namespace std;

/*--------------------------------------------------------------------------*/
/*--------------------- IMPLEMENTATION OF DCRLagrangianSolver---------------*/
/*--------------------------------------------------------------------------*/

/*<Risolve tramite line search il duale lagrangiano di DCR, rispetto al vincolo sul ritardo.
 *in particolare tiene una copia della struttura del grafo, su cui permette di aprire e chiudere
 *gli archi tramite il vettore ModCaps.
 *Crea a partire da questa una versione ridotta del grafo in cui elimina gli archi
 *con capacità inferiori ad r_min (passato in input al momento del caricamento), di questi archi
 *costruisce una versione semplificata con cui si risolverà SP*/

/*--------------------------------------------------------------------------*/
/*--------------------------- PUBLIC METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*----------------------------- CONSTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

  DCRLagrangianSolver:: DCRLagrangianSolver(void)
  {  
     eps = 1e-10;//Eps<double>();
     //cout<<"precisione: "<<eps<<endl;
     ObjVal = - Inf<double>();
     InterVal = 0; 
     lambda = 0;
     HeurVal = 0;

     maxCutSize = 50;

     LimitVal = 0;
     cardRedGraph = 0; 

     num_ite = 0;

     lagstat = Error;
     solvedflag = 0;
     inizialflag = 0;
     reoptflag = 0;
     
     timer = 0;
     RedGraLinks.resize(1);
     RedGraPos.resize(1);
     Linksp.resize(1);
     Links = 0;
     Nodes = 0;
     ModCaps.resize(1);
     //XSol = 0;
     numHops = 0;
     RSol.resize(1);
     RSolCosts.resize(1);
     SPLabels.resize(1);
     Cuts.resize(1);
     solneg.resize(1);
     solpos.resize(1);
     checkSolPos.resize(1);
     checkSolNeg.resize(1);
     OptPath.resize(1);

  }

/*--------------------------------------------------------------------------*/
/*----------------------------------LOAD DATA ------------------------------*/
/*--------------------------------------------------------------------------*/

   void DCRLagrangianSolver::LoadProblem (int nnodes, int nlinks, DCR::DCRFlow flow, DCR::DCRLink *links, DCR::DCRNode *nodes, double mtu, double r_min)
   {

      clean_up();

      /*una volta ripuliti i dati di eventuali chiamate precedenti
       *salviamo i nuovi dati e rimettiamo a zero le varie flag*/
      numNodes = nnodes;
      numLinks = nlinks;
      MTU = mtu;
      rmin = r_min;
      ObjVal = - Inf<double>();
      lambda = 0;
      HeurVal = 0;

      num_ite = 0;

      inizialflag = 0;
      solvedflag = 0;
      reoptflag = 0;
      nonposflag = 0;
      copyDataArray(flow, links, nodes);

      RSOLS.resize(nlinks);

   } 

/*--------------------------------------------------------------------------*/
   
  void DCRLagrangianSolver::LoadProblem (DCR::DCRFlow flow)
   {  

    ////FIXME: dovrei riaprire tutti gli archi che sono stati chiusi?
    //attualmente lo fa, modificare il secondo for commentando l'assegnamento su ModCaps
    //se non lo si desidera.
      int i;

    ///non potendo chiamare clean_up per 
    //perdere la struttura del grafo, ripuliamo esclusivamente le cose necessarie:
      for(i = 0; i < Cuts.size(); i++)
       {
         Cuts[i].RSol.clear();
       }

      Cuts.clear(); 
      

      ///inserimento dei nuovi dati
      ObjVal = - Inf<double>();
      lambda = 0;
      HeurVal = Inf<double>();

      inizialflag = 0;
      solvedflag = 0;
      reoptflag = 0;

      num_ite = 0;

      ModCaps.resize(numLinks);

      Flow = flow;

      for(i = 0; i < numLinks; i++)
       {
        Links[i].cost = Flow.costs[i];
        ModCaps[i] = Links[i].capacity; //riapre tutti gli archi
       }
   }
/*--------------------------------------------------------------------------*/
/*----------------------------------MODIFIERS-------------------------------*/
/*--------------------------------------------------------------------------*/

   void DCRLagrangianSolver::updrmin(double r_min)

   {
     
    // if(rmin-r_min > 0)
    //cout<<"è calato!"<<endl;

    if(rmin != r_min || solvedflag == 0) //altrimenti avevamo già risolto il problema per questo rmin
    {    
     lambda = 0;
     inizialflag = 0;

     solvedflag = 0;
     reoptflag = 1;  //proviamo reopt nel solve

     ObjVal = - Inf<double>();
     HeurVal = Inf<double>();

   /*  for(int i = 0; i < Cuts.size(); i++)
     cout<<"pendenza del taglio: "<<i<<" = "<<Cuts[i].c.m;

   cout<<endl;*/

    //eliminiamo gli archi per SP e la struttura
    /* ModCaps.clear();
     SPLabels.clear();
     RedGraLinks.clear();
     RedGraPos.clear();
     RSol.clear();
     RSolCosts.clear();
     Linksp.clear();*/

     rmin = r_min;
    }
    //else reoptflag = 0;
   // cout<<"lambda dopo updrmin"<<lambda<<endl;
    //si potrebbe richiedere else solvedflag = 1, ma dovrebbe essere già rimasto tale
    // dall'iterazione precedente.
   }

/*--------------------------------------------------------------------------*/

   void DCRLagrangianSolver::closeArcs(vector<int> arcs, int na)
   {
     int i;

     inizialflag = 0;

     for(i = 0; i < na; i++)
       ModCaps[arcs[i]] = 0;

  //evitiamo di modificare le capacità in Links, perché in tal caso
  //perderemmo l'informazione su quanto quella capacità fosse.
  //Sarebbe in tal caso un problema il riaprire gli archi chiusi.
  //Il vettore ModCaps tiene traccia di ogni modifica alle capacità.
   }

/*--------------------------------------------------------------------------*/

   void DCRLagrangianSolver::openArcs(int * arcs, int na)
   {
     int i;
    
     inizialflag = 0;

     for(i = 0; i < na; i++)
      ModCaps[arcs[i]] = Links[arcs[i]].capacity;
   }

/*--------------------------------------------------------------------------*/
/*--------------------------- GET METHODS-----------------------------------*/
/*--------------------------------------------------------------------------*/

/*<L'idea dietro al seguente metodo è che per verificare l'ammissibilità del lagrangiano,
 * basta la fase d'inizializzazione del problema, quindi se si volesse sapere solo questo, 
 * nel caso la risposta sia positiva, non vorremmo generare tutta l'informazione data
 * dal risolvere il problema.*/
   bool DCRLagrangianSolver::isFeasible(void)
   {
    Inizial();

    inizialflag = 1; //segnala che il problema è già stato inizializzato, nel caso
    //poi si voglia risolverlo con gli stessi dati, si può evitare di ripetere questa fase.

    if(lagstat == 0)
      return 0;
    else return 1;
   }

  /*--------------------------------------------------------------------------*/

//si risolve spt per lambda = 0 e lambda = eps, dal confronto dei valori della f.o. si conclude se lambda è ottimo

  bool DCRLagrangianSolver::is0opt(void)
   {
     vector<int> path;
     double beta_0 = Flow.burst / rmin - Flow.deadline;
     double alpha = 0;        //inizializzate   
     double beta = beta_0;    //al valore costante. 
     double nulobj = 0;
     double epsobj = 0;
     int i;

     lambda = 0;

     setReducedGraph();

     setSPTcosts();

     spt.updCosts(Linksp);
     spt.Solve(); //risolviamo SP con questo nuovo lambda.
       
     nhops = spt.getNHops();
         
     path.resize(nhops);
     path = spt.getPath(); //path prende gli indici del cammino ottimo di sp.

     for(i = 0; i < nhops; i++) //calcoliamo i valori di beta ed alpha.
      {
        beta += MTU / Linksp[path[i]].rstar + MTU / RedGraLinks[path[i]].speed + RedGraLinks[path[i]].delay + Nodes[RedGraLinks[path[i]].startnode].delay;
        alpha += RedGraLinks[path[i]].cost * Linksp[path[i]].rstar;
      }

    if(beta <= 0) 
      {
       //ci salviamo il valore di questa soluzione ammissbile:
       if(HeurVal >= alpha)  HeurVal = alpha;

       return 0; //allora lambda = 0 è ottimo
      }
    else //altrimenti proviamo per lambda = eps;
    {
      nulobj = alpha;
      lambda = 1e-6;

      alpha = 0;
      beta = beta_0;

      setSPTcosts();

      spt.updCosts(Linksp);
      spt.Solve(); //risolviamo SP con questo nuovo lambda.
       
      nhops = spt.getNHops();
         
      path.resize(nhops);
      path = spt.getPath(); //path prende gli indici del cammino ottimo di sp.

      for(i = 0; i < nhops; i++) //calcoliamo i valori di beta ed alpha.
       {
        beta += MTU / Linksp[path[i]].rstar + MTU / RedGraLinks[path[i]].speed + RedGraLinks[path[i]].delay + Nodes[RedGraLinks[path[i]].startnode].delay;
        alpha += RedGraLinks[path[i]].cost * Linksp[path[i]].rstar;
       }

      epsobj = alpha + lambda * beta;

      //cout<< " valore in eps "<<epsobj<< " valore in 0 "<<nulobj<<endl;

      if(epsobj > nulobj) return 1;
      else return 0;
    }

   }

  /*--------------------------------------------------------------------------*/

   double DCRLagrangianSolver::getLambda(void)
   {
    //cout<<"il lambda passato dal lagrangiano è"<<lambda<<endl;
    return(lambda);
   }

/*--------------------------------------------------------------------------*/


   vector<double> DCRLagrangianSolver::getCheckSol(void)
   {
    double molt;
    vector<double> checkSol;
    int i;


    if (nonposflag == 0) molt = betaneg / (betaneg - betapos);
    else  molt = 0;

    //molt = 0;
    cout<<"molt = "<<molt<<endl;
   
    cout<<"sol neg: ";
    for(i=0; i<solneg.size(); i++)
     {cout<<solneg[i]<<"-";}
    cout<<endl;

    cout<<"sol pos: ";
    for(i=0; i<solpos.size(); i++)
     {cout<<solpos[i]<<"-";}
    cout<<endl;

    checkSol.resize(solneg.size());
    
   // cout<<"soluzione di taglia"<< checkSol.size();

    for( i = 0; i < checkSol.size(); i++)
       checkSol[i] = (1-molt)*solpos[i] + (molt)*solneg[i];

    return(checkSol);
   }


/*--------------------------------------------------------------------------*/

   vector<int> DCRLagrangianSolver::getCheckSolPos()
   {
    return(checkSolPos);
   }


/*--------------------------------------------------------------------------*/

   vector<int> DCRLagrangianSolver::getCheckSolNeg()
   {
    return(checkSolNeg);
   }


/*--------------------------------------------------------------------------*/

   double DCRLagrangianSolver::getInterVal(void) //non dovrebbe essere davvero necessaria
   {
   	return(InterVal);
   }

/*--------------------------------------------------------------------------*/

   double DCRLagrangianSolver::getOptBeta()
   {
    return(optBeta);
   }

/*--------------------------------------------------------------------------*/

   double DCRLagrangianSolver::getObjVal(void)
   {
    return(ObjVal);
   }

/*--------------------------------------------------------------------------*/

   vector<double> DCRLagrangianSolver::getSPLabels(void)
   {
    return(SPLabels);
   }

/*--------------------------------------------------------------------------*/

   int DCRLagrangianSolver::getnumHops()
   {
    return(numHops);
   }

/*--------------------------------------------------------------------------*/  


   vector<int> DCRLagrangianSolver::getPath()
   {
    return(OptPath);
   } 


/*--------------------------------------------------------------------------*/

   vector<double> DCRLagrangianSolver::getRSol()
   {
    return(RSol);
   }

/*--------------------------------------------------------------------------*/

   vector<double> DCRLagrangianSolver::getRSolCosts()
   {
    return(RSolCosts);
   }


/*--------------------------------------------------------------------------*/  


   vector<int> DCRLagrangianSolver::getRedGraphPos()
   {
    return(RedGraPos);
   } 

/*--------------------------------------------------------------------------*/  

   int DCRLagrangianSolver::getCardRedGraph()
   {
    return(cardRedGraph);
   }


/*--------------------------------------------------------------------------*/  


   vector<int> DCRLagrangianSolver::getComplGraph()
   {
    return(ComplGraph);
   } 


/*--------------------------------------------------------------------------*/  

   double DCRLagrangianSolver::getHeurVal()
   {
    return(HeurVal);
   }

/*--------------------------------------------------------------------------*/ 
  /* vector<DCRLagrangianSolver::Cut_Val> DCRLagrangianSolver::getCuts(void)
   {
   	return(Cuts);
   }*/


/*--------------------------------------------------------------------------*/   


   double DCRLagrangianSolver::DCRgetTime( void )
    {   
     return( timer ? timer->Read() : 0 );
    }
  

/*--------------------------------------------------------------------------*/

   DCRLagrangianSolver::LAGStatus DCRLagrangianSolver::getStatus(void)
   {
    return(lagstat);
   }

/*--------------------------------------------------------------------------*/

   int DCRLagrangianSolver::getNumIte(void)
   {
    return(num_ite);
   }

/*--------------------------------------------------------------------------*/
/*--------------------------- SET TIME METHODS-----------------------------------*/
/*--------------------------------------------------------------------------*/

 void DCRLagrangianSolver::DCRsetTime( bool timeON )
  {
    if( timeON )
      if(timer) timer->ReSet();
      else timer = new OPTtimers();
    else
      delete timer; 
    }
    
    /**< If timeON is true sets or resets the timer, if false deletes the timer.
     * \param timeON bool value */

/*--------------------------------------------------------------------------*/

 void DCRLagrangianSolver::DCRstartTime( void )
  {
    if( timer == 0 )
      {
        cout<<"DCRLagrangianSolver::DCRstartTimer(): failed to start timer"; 
        exit(1);
      }
    
    timer->Start();
    }
    
    /**< If timer was set or re-set, starts or re-starts timer ticking. */

    
/*--------------------------------------------------------------------------*/

 void DCRLagrangianSolver::DCRstopTime( void )
  {
    if ( timer == 0 )
     {
       cout<<"DCRLagrangianSolverL::DCRstopTimer(): failed to stop timer";
       exit(1);
     }

    timer->Stop();
  }
    
    /**< If timer was set, stops timer ticking. */

/*--------------------------------------------------------------------------*/

 void DCRLagrangianSolver::DCRsetTimeLimit( long secs )
  {
    tlimit = secs;
  } 
  
   /**< Sets a time limit for the solver.
    * \param secs long expressing timelimit in seconds */ 

/*--------------------------------------------------------------------------*/
/*----------------------------------SOLVE-----------------------------------*/
/*--------------------------------------------------------------------------*/

   void DCRLagrangianSolver::Solve(void)
   {
    int i;
    int opt = 1; //avverte se la funzione reopt ha sostituito l'inizializzazione
    vector<int> path;


    double releps; //prende il max tra objval e 1, per la precisione relativa.
    
    double alpha;
    double beta;
    double beta_o = Flow.burst / rmin - Flow.deadline; //per evitare di ricalcolare il valore
    
    reoptflag = 0; //DEBUG: spegne la riottimizzazione
    //inizialflag == 0;

    if(reoptflag == 1)
    {
      opt = Reopt();     
    }
    
    //cout<<"opt = "<<opt<<" inizialflag = "<<inizialflag<<endl;
    if(opt == 1)
    {

      if(inizialflag == 0)//potrebbe essere già stato inizializzato dal metodo isFeasible
        {
          Inizial();
          inizialflag = 1;
        }
    }
    else 
      {
        setReducedGraph();
        //setSPTcosts();
        //spt.LoadProblem(numNodes, cardRedGraph, Linksp, Flow.sourcenode, Flow.sinknode);
      }
    //cout<<"inizializzato"<<endl;
    //cout<<"lagstat ="<<lagstat<<endl;

    if(lagstat == OK && solvedflag != 1)
    {
     do
       {         
        ////variabili per salvare il taglio
        ///generate a ogni iterazione perché allungare il loro tempo di vita
        //comporta un deterioramento dei dati precedentemente immessi nel vettore Cuts.
        
        DCRLagrangianSolver::Cut_Val cut; //taglio e sol da inserire nei tagli salvati

        //vector<double> tempRsol; //soluzioni generate da salvare insieme ai tagli
         //cout<<"beta_0 = "<<beta<<"  "<<betainizial<<"  "<<Flow.burst / rmin - Flow.deadline<<endl;
        //intercetta e pendenza del nuovo taglio
        alpha = 0;       //inizializzate   
        beta = beta_o; //al valore costante. 

         lambda = (mCut.q - pCut.q) / (pCut.m - mCut.m);  //lambda prende il punto d'intersezione tra i tagli.
           
         if(lambda < 0) //aggiungiamo alla LS anche il vincolo lineare \lambda = 0;
          {lambda = 0; /*cout<<"sto imponendo lambda = 0";*/}

         //cout<<"all'interno di lagrange come evolve lambda?"<<lambda<<endl;
         InterVal = mCut.q + lambda * mCut.m; //vecchia retta (è indifferente quale) valutata nel nuovo punto.
        //cout<<"InterVal = "<<InterVal<<endl;
        //vector<double> tempRsol;

         //cout<<"1beta_0 ="<<beta_o<<endl;
         setSPTcosts(); //nuovo lambda, nuovi costi.
         spt.updCosts(Linksp);
         spt.Solve(); //risolviamo SP con questo nuovo lambda.
       
         nhops = spt.getNHops();
         
         path.resize(nhops);
         cut.RSol.resize(nhops);
        //tempRsol.resize(nhops);
         
       
         // cout<<"2beta_0 ="<<beta_o<<endl;
         path = spt.getPath(); //path prende gli indici del cammino ottimo di sp.

         for(i = 0; i < nhops; i++) //calcoliamo i nuovi valori di beta ed alpha.
          {
            beta += MTU / Linksp[path[i]].rstar + MTU / RedGraLinks[path[i]].speed + RedGraLinks[path[i]].delay + Nodes[RedGraLinks[path[i]].startnode].delay;
            alpha += RedGraLinks[path[i]].cost * Linksp[path[i]].rstar;
            //tempRsol[i] = Linksp[path[i]].rstar;
            cut.RSol[i] = Linksp[path[i]].rstar;

          }

         optBeta = beta;

         if(beta < -1e-20 && alpha <= HeurVal)
           {HeurVal = alpha;} //salviamo la migliore soluzione primale ammissibile
         
          /*
         /////salviamo i dati per il CHECK su CPLEX!!!!
         if(beta < 0)
         { 
           solneg.resize(nhops);
           
           checkSolNeg.resize(nhops);

           for(i = 0; i < nhops; i++)
           {
            solneg[i] = Linksp[path[i]].rstar;
            checkSolNeg[i] = RedGraPos[path[i]];
           }
           betaneg = beta;
         }
         else
         {
           solpos.resize(nhops);
           checkSolPos.resize(nhops);

           for(i = 0; i < nhops; i++)
            {
              solpos[i] = Linksp[path[i]].rstar;
             checkSolPos[i] = RedGraPos[path[i]];
           }
           betapos = beta;

         }
         ///fine del salvataggio, commentare in futuro per non appesantire in memoria
        */
         cut.c.m = beta;
         cut.c.q = alpha;

         cut.RSolsize = nhops;  
         cut.rmin = rmin;     
        
         if(Cuts.size() < maxCutSize)
           Cuts.push_back(cut); //salviamo il nuovo taglio

         UpdCut(alpha, beta); //aggiorniamo il valore di uno dei tagli ottimi.

         ObjVal = alpha + lambda * beta; //nuova retta nel punto, quindi val della funzione obiettivo.

         if(InterVal > 1) releps = InterVal;
         else releps = 1;
         
         //cout<<"in solve differenza assoluta"<<(InterVal - ObjVal)<<endl;
         //cout<<"soglia relativa"<<eps*releps<<endl<<endl;
         num_ite++;

      }while((InterVal - ObjVal) > (eps * releps)); ////differenza > epsilon * max{1,InterVal}
    
    //cout<<" eps = "<<eps<<endl;
    //cout<<"in solve differenza relativa"<<(InterVal - ObjVal)/InterVal<<endl;
     //cout<<endl<<endl<<"beta = "<<beta<<endl<<endl;
     //cout<<endl<<endl<<"alpha = "<<alpha<<endl<<endl;
    
    //cout<<"numero iterazioni della line search"<<counter<<endl;
    //cout<<"oltre l'iterazione?"<<endl;

    SPLabels.resize(numNodes);

    vector<double> labels;
    
    labels = spt.getLabel(); //salviamo i potenziali.

    for(i = 0; i < numNodes; i++)
        SPLabels[i] = labels[i];

   //cout<<"stampi?"<<endl;

    RSol.resize(nhops);
    RSolCosts.resize(nhops);

    for(i = 0; i < numLinks; i++) 
      RSOLS[i] = 0.0;
    
    for(i = 0; i < nhops; i++)  //salviamo le coppie (r_ij,f_ij) della soluzione ottima
      {
        RSol[i] = Linksp[path[i]].rstar;
        //std::cout << i << "," << RSol[i] << std::endl;
        RSolCosts[i] = Linksp[path[i]].cost;
      }

    for(int i = 0; i < numLinks; i++) 
      for(int j = 0; j < nhops; j++)
        if(i == path[j])
          RSOLS[i] = RSol[j];
        

    numHops = nhops;
    solvedflag = 1;
    }//se non era stato risolto e non è considerato già unfeasible
   //if(solvedflag == 1)
   //cout<<"numero di tagli salvati"<<Cuts.size()<<endl;
     //cout<<"lambda ottimo: "<<lambda<<" risolto con stato: "<<lagstat<<endl;

   path.clear();
   //cout<<"non ho fatto l'iterazione"<<endl;

   }

/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

   DCRLagrangianSolver::~DCRLagrangianSolver()
   {
	  clean_up();
   }

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

   void DCRLagrangianSolver::copyDataArray(DCR::DCRFlow flow, DCR::DCRLink *links, DCR::DCRNode *nodes)
   {
    int i;
    Links = new DCR::DCRLink[numLinks];
    ModCaps.resize(numLinks);    
    Nodes = new DCR::DCRNode[numNodes];

    Flow = flow;

    for(i = 0; i < numLinks; i++)
    {
      Links[i].startnode = links[i].startnode;
      Links[i].endnode = links[i].endnode;
      Links[i].speed = links[i].speed;
      Links[i].capacity = links[i].capacity;
      Links[i].delay = links[i].delay;
      Links[i].cost = links[i].cost; //Flow.costs[i];
      ModCaps[i] = Links[i].capacity;
    }
   
    //cout<<"una capacità="<<Links[0].capacity<<endl;
    //cout<<"un costo="<<Links[0].cost<<endl;
    for(i = 0; i < numNodes; i++)
    {
      Nodes[i].delay = nodes[i].delay;
    }
   }


/*--------------------------------------------------------------------------*/

   void DCRLagrangianSolver::getLimitVal(void)
   {
    int i;

    LimitVal = 0;

    for(i = 0; i < cardRedGraph; i++)
      {
        if(RedGraLinks[i].cost > 0)
          LimitVal += RedGraLinks[i].cost * RedGraLinks[i].capacity; //saturiamo tutti gli archi di costo positivo
      }                                                              //e ne sommiamo i costi.
   }
 /*upper bound al valore di ogni soluzione ammissibile*/  

/*--------------------------------------------------------------------------*/

   double DCRLagrangianSolver::getCost(double r, int lindex)
   {    
    double c;
    double barl = MTU / RedGraLinks[lindex].speed + RedGraLinks[lindex].delay + Nodes[RedGraLinks[lindex].startnode].delay;
    c = lambda * barl + RedGraLinks[lindex].cost * r + lambda * MTU / r;

    return c; 
   }
/*<restituisce i costi per l'archi lindex di SP in funzione di r^*_ij, passato in input come r*/ 

/*--------------------------------------------------------------------------*/

   int DCRLagrangianSolver::Reopt(void)
   {

    /*FIXME: trova dei valori per lambda per cui l'intersezione tra i tagli
     stia sotto il valore della funzione obiettivo.
     Siamo sicuri che i tagli siano salvati in modo sicuro? (sembrerebbe di sì)*/

     int isnegcut = 1; //flag per capire se ci sono tagli negativi
     int negpos; //posizione del taglio negativo in Cuts
     int isposcut = 1; //lo stesso per tagli positivi
     int pospos;
     double epsrel; //fattore per la precisione relativa
     int i, k;
     int delflag = 0; //flag per la cancellazione di soluzioni non ammissibili
     double someinter;
     double value;
     double approx;
     double min = Inf<double>();
     int minpos;
     int counterdeleted = 0;
     int maxiter = 0;
     
     //cout<<"tagli presenti"<<Cuts.size()<<endl;
     if(Cuts.size() <= 1) return 1; //solo un taglio o meno è troppo poco per riottimizzare

     for(i = 0; i < Cuts.size(); i++) //controlliamo quali soluzioni sono ancora ammissibili
     {
      //cout<<"nhops della sol del taglio "<<i<<" = "<<Cuts[i].RSolsize<<endl;
       for(k = 0; k < Cuts[i].RSolsize; k++)
        {
          if(Cuts[i].RSol[k] < rmin) //basta accada anche per un solo r_ij
          {
            //cout<<Cuts[i].RSol[k]<<"-"<<rmin<<endl;
            delflag = 1; //in tal caso non sarà ammissibile
          }

        }
       if(delflag == 1)
        {
         counterdeleted++;
         Cuts.erase(Cuts.begin()+i); //e verrà cancellata
         i--;
         delflag = 0;
        }
       else //se il taglio è valido ne controlliamo la pendenza
        //FIXME: si potrebbe aggiungere una sensibilità
       {
        //correggiamo la pendenza del taglio con il nuovo valore di rmin
        Cuts[i].c.m = Cuts[i].c.m - Flow.burst / Cuts[i].rmin + Flow.burst / rmin;
        //e ci segnamo che la correzione per questo taglio è avvenuta
        Cuts[i].rmin = rmin;

        if(Cuts[i].c.m < 0) 
          {
            isnegcut = 0;
            negpos = i;
          }
        else 
          {
            if(Cuts[i].c.m >= 0) //dovrebbe essere superfluo
            {
              isposcut = 0;
              pospos = i;
            }
          }

       }
     }
     //cout<<"tagli cancellati"<<counterdeleted<<endl;
     if(Cuts.size() <= 1) return 1;
     //cout<<"ne sono rimasti "<<Cuts.size()<<endl;

    if(isnegcut == 1 || isposcut == 1) return 1; //se non sono sopravvissuti tagli positivi o negativi rinunciamo
  
  ////line search per la ricerca dell'intersezione ottima tra i tagli ancora ammissibili
 
    pCut.m = Cuts[pospos].c.m;
    pCut.q = Cuts[pospos].c.q;
    mCut.m = Cuts[negpos].c.m;
    mCut.q = Cuts[negpos].c.q;

    do
    {
      lambda = (mCut.q - pCut.q) / (pCut.m - mCut.m);
      
      ////FIXME: capire perché succede!!!!!!!
      if(lambda < 0 || abs(pCut.m - mCut.m) < 1e-20)
       {//cout<<"Bug in reopt"<<endl;
        lambda = 0;
        return 1;
       }
      approx = mCut.q + lambda * mCut.m;

      for(i = 0; i < Cuts.size(); i++) //per calcolare il valore della funzione per quel lambda
       {                               //vediamo per quale taglio si ottiene il valore più basso.
         someinter = Cuts[i].c.q + lambda * Cuts[i].c.m;

         if(min >= someinter)
           {
            min = someinter;
            minpos = i;
           } 
       }

      value = Cuts[minpos].c.q + lambda * Cuts[minpos].c.m;

      UpdCut(Cuts[minpos].c.q, Cuts[minpos].c.m);

      if(approx > 1) epsrel = approx;
      else epsrel = 1;

      //cout<<"ecco il loop"<<endl;
      maxiter ++;
      //cout<<"approx - value"<<approx - value<<endl;
    }while(approx - value > eps*epsrel && maxiter <= 30);
    //FIXME: va in loop con garr199904, lanciata da BenBound, al flusso 78
   // cout<<"approx || value : "<<approx<<"||"<<value<<endl;
    //cout<<"niter = "<<maxiter<<endl;
    if(maxiter >= 30) 
      {
        cout<<"guarda che è andato in loop il reopt"<<endl;
        return 1;
      }
   // cout<<"finisce il reopt con valore = "<<value<<endl;
    //cout<<"tagli pubblici:"<<pCut.m<<"||"<<mCut.m<<endl;
    return 0;
   }

/*--------------------------------------------------------------------------*/

   void DCRLagrangianSolver::Inizial(void)
   {
    int i, j;
    vector<int> path;
   // vector<double> tempRsol;
    double beta_0 = Flow.burst / rmin - Flow.deadline;
    double alpha = 0;        //inizializzate   
    double beta = beta_0;    //al valore costante. 

    //cout<<"lambda = "<<lambda<<endl;
    //cout<<"beta = "<<beta<<endl;
    setReducedGraph(); //lavoriamo sul grafo ridotto.
    
    setSPTcosts();

    //cout<<"cardRedGraph="<<cardRedGraph<<endl;
    //cout<<Flow.sourcenode<<Flow.sinknode<<endl;

    spt.LoadProblem(numNodes, cardRedGraph, Linksp, Flow.sourcenode, Flow.sinknode);
    spt.Solve();
    
    spstat = spt.getStatus(); 
    //cout<<"spstat è "<<spstat<<endl;
    if(spstat  == 0) //il grafo ridotto non ha cicli di costo negativo, né è sconnesso.
    { 
      getLimitVal();

      nhops = spt.getNHops();
      path.resize(nhops);
      path = spt.getPath(); //path prende gli indici del cammino ottimo di sp

     // cout<<"path nel lagrangiano e rispettivi rstar: ";

      for(i = 0; i < nhops; i++) //calcoliamo i valori di beta ed alpha.
        {
          beta += MTU / Linksp[path[i]].rstar + MTU / RedGraLinks[path[i]].speed + RedGraLinks[path[i]].delay + Nodes[RedGraLinks[path[i]].startnode].delay;
          alpha += RedGraLinks[path[i]].cost * Linksp[path[i]].rstar;
         // cout<<"barl in lagr per l'arco "<<i<<" è "<<MTU / RedGraLinks[path[i]].speed + RedGraLinks[path[i]].delay + Nodes[RedGraLinks[path[i]].startnode].delay<<endl;
         // cout<<"ed il suo rstar è "<<Linksp[path[i]].rstar<<endl;
        }
    // cout<<endl;
     //cout<<"alpha= "<<alpha<<endl;
     //cout<<"beta = "<<beta<<endl;

      UpdCut(alpha,beta); //aggiorniamo il primo taglio ottimo.

      /////verifichiamo se 0 è già il lambda ottimo
      if(beta <= 0) 
        { 
          DCRLagrangianSolver::Cut_Val cut;
         //cout<<"pendenza negativa"<<endl;
          lagstat = OK;
          InterVal = alpha;
          ObjVal = alpha;
          lambda = 0;
          optBeta = beta;

          HeurVal = alpha; // ci salviamo il valore della f.o. in concomitanza di beta<0

          cut.RSol.resize(nhops);
          path.resize(nhops);
          OptPath.resize(nhops);

          path = spt.getPath();


          RSol.resize(nhops);
          RSolCosts.resize(nhops);

          /*
          ////SALVIAMO I DATI PER IL CHECK SU CPLEX
          solneg.resize(nhops);
          checkSolNeg.resize(nhops);

          for(i = 0; i < nhops; i++)
           {
            solneg[i] = Linksp[path[i]].rstar;
            checkSolNeg[i] = RedGraPos[path[i]];
           }

          betaneg = beta;
          nonposflag = 1;
          //fine salvataggio
          */

          for(i = 0; i < nhops; i++)  //salviamo le coppie (r_ij,f_ij) della soluzione ottima
          {
           RSol[i] = Linksp[path[i]].rstar;
           cut.RSol[i] = Linksp[path[i]].rstar;
           RSolCosts[i] = Linksp[path[i]].cost;
           OptPath[i] = path[i];
          }

          for(i = 0; i < numLinks; i++) 
            RSOLS[i] = 0.0;
               
          for(int i = 0; i < numLinks; i++) 
            for(int j = 0; j < nhops; j++)
              if(i == path[j])
                RSOLS[i] = RSol[j];
          
          cut.RSolsize = nhops;
          cut.c.m = beta;
          cut.c.q = alpha;

          cut.rmin = rmin;
          
          if(Cuts.size() < maxCutSize)
          Cuts.push_back(cut); //salviamo il taglio con pendenza negativa.

          SPLabels.resize(numNodes);
          vector<double> labels;
    
          labels = spt.getLabel(); //salviamo i potenziali.

          for(i = 0; i < numNodes; i++)
            SPLabels[i] = labels[i];

          solvedflag = 1;
        }
      else  //////ricerca del taglio con pendenza negativa
        {
          //cout<<"ricerca del taglio con pendenza negativa"<<endl;
         lambda = 1;
        
         while(beta > 0 && ObjVal <= LimitVal && spstat == 0)
           {
            alpha = 0;           //inizializzate a ogni iterazione 
            beta = beta_0;

            setSPTcosts();
            spt.updCosts(Linksp);
            spt.Solve();
            spstat = spt.getStatus();
            nhops = spt.getNHops();

            path.resize(nhops);
            path = spt.getPath(); //path prende gli indici del cammino ottimo di sp
            //cout<<path[0]<<path[1]<<endl;
            for(i = 0; i < nhops; i++) //calcolo valori di beta ed alpha.
              {  
                 beta += MTU / Linksp[path[i]].rstar + MTU / RedGraLinks[path[i]].speed + RedGraLinks[path[i]].delay + Nodes[RedGraLinks[path[i]].startnode].delay;
                 alpha += RedGraLinks[path[i]].cost * Linksp[path[i]].rstar;
              } 

            //cout<<"beta = "<<beta<<endl;
            //cout<<"ObjVal ="<<ObjVal<<endl;

            ObjVal = alpha + lambda * beta; 
            lambda = 2 * lambda; //si possono scegliere anche successioni diverse
          } // while ricerca del beta negativo
          
          DCRLagrangianSolver::Cut_Val cut;        
          //tempRsol.resize(nhops);
          cut.RSol.resize(nhops);

          for(i = 0; i < nhops; i++)
             {
              //tempRsol[i] = Linksp[path[i]].rstar;
              cut.RSol[i] = Linksp[path[i]].rstar;
             }
          
          /*/
          ////SALVIAMO I DATI PER IL CHECK SU CPLEX
          solneg.resize(nhops);
          checkSolNeg.resize(nhops);

          for(i = 0; i < nhops; i++)
           {
            solneg[i] = Linksp[path[i]].rstar;
            checkSolNeg[i] = RedGraPos[path[i]];
           }

          betaneg = beta;
          nonposflag = 1;
          //fine salvataggio
          */

          cut.c.m = beta;
          cut.c.q = alpha;         

          /*for(i = 0; i < nhops; i++)
           cut.RSol[i] = tempRsol[i];*/

          cut.RSolsize = nhops;
          cut.rmin = rmin;

          if(Cuts.size() < maxCutSize)
          Cuts.push_back(cut); //salviamo il taglio con pendenza negativa.

          UpdCut(alpha, beta); //aggiorniamo il secondo taglio ottimo.
          //cout<<"LimitVal="<<LimitVal<<endl;
          //cout<<"ObjVal="<<ObjVal<<endl;

          if(ObjVal > LimitVal)   {lagstat = Infeasible; } // allora il lagrangiano è vuoto. ?? ma se siamo partiti con un
          						  //taglio, almeno un punto ci sarà, potrebbe significare fun crescente
          else lagstat = OK;
        }
          //cout<<"beta = "<<beta<<endl;
    } // se SP ha funzionato
    else lagstat  = Infeasible;
    path.clear();

   }


/*--------------------------------------------------------------------------*/
 
 ////FIXME:sarebbe meglio farlo con un solo scorrimento, ma viene fatto solo una volta
   void DCRLagrangianSolver::setReducedGraph()
   {
     int i, j = 0, h = 0;
     cardRedGraph = 0; //per le chiamate successive
    
     for(i = 0; i < numLinks; i++)
     {
       if(rmin <= ModCaps[i]) //questi sono quelli su cui lavoreremo (r_min<=c_ij)
         cardRedGraph ++;
     }
    
 //FIXME: mi serve davvero RedGraLinks, adesso che salvo le sue posizioni?
     Linksp.resize(cardRedGraph);
     RedGraLinks.resize(cardRedGraph);
     RedGraPos.resize(cardRedGraph);
     //ComplGraph.resize(numLinks - cardRedGraph);

     for(i = 0; i < numLinks; i++)
     {
       if(rmin <= ModCaps[i]) //altrimenti  SP non deve considerarli (i.e. x_ij=0).
       {
         RedGraLinks[j].startnode = Links[i].startnode;
         RedGraLinks[j].endnode = Links[i].endnode;
         RedGraLinks[j].capacity = Links[i].capacity;
         RedGraLinks[j].speed = Links[i].speed;
         RedGraLinks[j].delay = Links[i].delay;
         RedGraLinks[j].cost = Links[i].cost;
         Linksp[j].startnode = Links[i].startnode;
         Linksp[j].endnode = Links[i].endnode;
         RedGraPos[j] = i;
         j++;
       }
      /* else
       {
         ComplGraph[h] = i;
         h++;
       }*/
     }
   }

/*--------------------------------------------------------------------------*/

   void DCRLagrangianSolver::setSPTcosts()
   {
     int i;
     double sqr;
     
     for(i = 0; i < cardRedGraph; i++)
     {
      ///scelta dei costi in base ai vari casi, si veda paragrafo 2.
      //cout<<"set STPcost"<<endl;
      if(RedGraLinks[i].cost < 0) 
        {
          Linksp[i].rstar = RedGraLinks[i].capacity;
          Linksp[i].cost = getCost(RedGraLinks[i].capacity, i);
        }
      else
       { 
         sqr = sqrt(lambda * MTU / RedGraLinks[i].cost);
 
         if(sqr < rmin)  
           { 
              Linksp[i].rstar = rmin;
              Linksp[i].cost = getCost(rmin, i); 
           }
         else
          {
            if(RedGraLinks[i].capacity < sqr)  
             { 
                Linksp[i].rstar = RedGraLinks[i].capacity;
                Linksp[i].cost = getCost(RedGraLinks[i].capacity, i);
             }
            else 
             { //cout<<"scelgo"<<sqr<<endl;
                Linksp[i].rstar = sqr;
                Linksp[i].cost = getCost(sqr, i);
             }//caso rmin < sqr < c_ij
          }
       } //fine suddivisione casi per costi.
     }//for all reduced graph arcs
    }

/*--------------------------------------------------------------------------*/
 
   void DCRLagrangianSolver::clean_up()
   {
     int i;

     if(Links != 0) {delete [] Links; Links = 0;}
     if(Nodes != 0) {delete [] Nodes; Nodes = 0;}
     //if(timer != 0) {delete [] timer; timer = 0;}
     ModCaps.clear();
     SPLabels.clear();
     RedGraLinks.clear();
     RedGraPos.clear();
     RSol.clear();
     RSolCosts.clear();
     Linksp.clear();
     solpos.clear();
     solneg.clear();
     checkSolPos.clear();
     checkSolNeg.clear();
     OptPath.clear();
     RSOLS.clear();

     for(i = 0; i < Cuts.size(); i++)
       {
         Cuts[i].RSol.clear();
       }
     
     Cuts.clear(); 
   }  

/*--------------------------------------------------------------------------*/
/*----------------------------------MODIFIERS-------------------------------*/
/*--------------------------------------------------------------------------*/


   void DCRLagrangianSolver::UpdCut(double alpha, double beta)
   {
    if(beta>=0)
    	{
        pCut.m = beta;
        pCut.q = alpha;
      }

    else
      {
        mCut.m = beta;
        mCut.q = alpha;
      }
   }
/*<modifica uno dei due tagli ottimi che definiscono la soluzione al momento a seconda della pendenza*/


/*--------------------------------------------------------------------------*/
/*-------------------------- End File DCRLagrangianSolver.cpp -------------------------*/
/*--------------------------------------------------------------------------*/
