/*--------------------------------------------------------------------------*/
/*-------------------------- File BenBound.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------------- IMPLEMENTATION -------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "BenBound.h"
#include "DCRLagrangianSolver.h"
#include "DCR.h"
#include "OPTUtils.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <ctime>
#include <cstdlib>

using namespace std;

/*--------------------------------------------------------------------------*/
/*--------------------- IMPLEMENTATION OF BenBound--------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------------- PUBLIC METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*----------------------------- CONSTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

  BenBound::BenBound()
  {
  	eps = 1e-6;
  	BenStat = Error;

    timer = 0;

  	Links = 0;
  	Nodes = 0;
  	caps = 0;
  	solvedflag = 0;

    myparam = 0.995; //FIXME:creare un metodo pubblico per l'upd (before it was set to 0.995)
    
    SPLabels.resize(1);
  	Q.resize(1);
  }

/*--------------------------------------------------------------------------*/
/*----------------------------------LOAD DATA ------------------------------*/
/*--------------------------------------------------------------------------*/

  void BenBound::LoadProblem (int nnodes, int nlinks, DCR::DCRFlow flow, DCR::DCRLink *links, DCR::DCRNode *nodes, double mtu)
   {
      clean_up();

      numNodes = nnodes;
      numLinks = nlinks;
      MTU = mtu;
      ObjVal = Inf<double>();
      HeurVal = Inf<double>();
      BestUB = Inf<double>();
      BestLB = -Inf<double>();
      ApproxVal = -Inf<double>();
      solvedflag = 0;
      counter_ite_Ben = 0;
      counter_ite_Lag = 0;
      is_convex_iteration = false;

      SPLabels.resize(numNodes);

      copyDataArray(flow, links, nodes);

      SOLUTION.resize(numLinks);

      lagSol.LoadProblem(numNodes, numLinks, Flow, Links, Nodes, MTU, 1);

   } 

/*--------------------------------------------------------------------------*/

  void BenBound::LoadProblem (DCR::DCRFlow flow)
   {
      int i;

      for(int i = 0; i < Q.size(); i++)
        Q[i].Cuts.clear();
      Q.clear();

      ObjVal = Inf<double>();
      HeurVal = Inf<double>();
      BestUB = Inf<double>();
      BestLB = -Inf<double>();
      ApproxVal = -Inf<double>();
      solvedflag = 0;
      counter_ite_Ben = 0;
      counter_ite_Lag = 0;
      is_convex_iteration = false;

      Flow = flow;

      for(i = 0; i < numLinks; i++)
        Links[i].cost = Flow.costs[i];

      lagSol.LoadProblem(numNodes, numLinks, Flow, Links, Nodes, MTU, 1);

   } 

//metodo per il caricamento, a parità di topologia, del solo nuovo flusso da servire.

/*--------------------------------------------------------------------------*/
/*--------------------------- GET METHODS-----------------------------------*/
/*--------------------------------------------------------------------------*/

   double BenBound::getObjVal()
   {
     //return(SOL_VALUE);
     return(ObjVal);
   }

   
/*--------------------------------------------------------------------------*/


   double BenBound::getr_min()
   {
     return(rmin);
   }


/*--------------------------------------------------------------------------*/


   double BenBound::getUB()
   {
     //std::cout << "UB=" << BestUB << " " << " " << ObjVal << 
     //   " " << " " << HeurVal << " " << ApproxVal<< std::endl;
     //return(std::max(BestUB,HeurVal));
     return(BestUB);
   }

/*--------------------------------------------------------------------------*/

   double BenBound::getLB()
   {
     //std::cout << "LB=" << BestLB << std::endl;
     //std::cout << abs(ApproxVal-HeurVal) << std::endl;
     //std::cout<<"BenStat="<<BenStat<<std::endl;
     //if(abs(ApproxVal-HeurVal)>eps*abs(ObjVal))
      //BestLB = -1e+301;
     return(BestLB);
     //return(ObjVal);
   }


/*--------------------------------------------------------------------------*/


   BenBound::BndrStat BenBound::getStat()
   {
     return(BenStat);
   }

/*--------------------------------------------------------------------------*/

  double BenBound::getHeurVal()
   {
    return(HeurVal);
   }

/*--------------------------------------------------------------------------*/

  double BenBound::getSolution( int i )
   { 
    return(SOLUTION[ i ]);
   }

/*--------------------------------------------------------------------------*/

  double BenBound::getApproxVal()
   {
    return(ApproxVal);
   }

/*--------------------------------------------------------------------------*/

  int BenBound::getNumIterationBender()
   {
    return(counter_ite_Ben);
   }

/*--------------------------------------------------------------------------*/

  int BenBound::getNumIterationLagr()
   {
    return(counter_ite_Lag);
   }

/*--------------------------------------------------------------------------*/

  bool BenBound::IsConvexInteration()
   {
    return(is_convex_iteration);
   }

/*--------------------------------------------------------------------------*/

   vector<double> BenBound::getxPlotData()
   {
     int i;
     numPoints = 1000;

     xPlot.resize(numPoints);

     Inizial();

     for(i = 1; i < numPoints + 1; i++)
     {
       xPlot[i-1] = (i*(Q[1].rmin - Q[1].inter)/numPoints)+Q[1].inter; //salviamo la divisione in 100 punti dell'intervallo
     }
     
    return(xPlot);
   }


/*--------------------------------------------------------------------------*/

   vector<double> BenBound::getyPlotData()
   {
     int i;
     double lambda;
     vector<double> yPlot;
     vector<int> RedGraph; //grafo ridotto che ci faremo passare dal Lagrangian solver
     int cardRedGraph;

     yPlot.resize(numPoints);
     
     for(i = 0; i < numPoints; i++)
     {
      lagSol.updrmin(xPlot[i]);

      lagSol.Solve();

      lambda = lagSol.getLambda();
      SPLabels = lagSol.getSPLabels();
      cardRedGraph = lagSol.getCardRedGraph();
      RedGraph.resize(cardRedGraph);
      RedGraph = lagSol.getRedGraphPos();

      yPlot[i] = lambda * Flow.burst / xPlot[i] - lambda * Flow.deadline + SPLabels[Flow.sinknode]; //valore della fo in xplot[i]
     }
    
    return yPlot;
   }


/*--------------------------------------------------------------------------*/

  double BenBound::getTime( void ) const
   {   
     return( timer ? timer->Read() : 0 );
   }


/*--------------------------------------------------------------------------*/
/*--------------------------- SET TIME METHODS-----------------------------------*/
/*--------------------------------------------------------------------------*/

 void BenBound::DCRsetTime( bool timeON )
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

 void BenBound::DCRstartTime( void )
  {
    if( timer == 0 )
      {
        /*cout<<"CPX_LAG::DCRstartTimer(): failed to start timer";*/
        exit(1);
      }
    
    timer->Start();
    }
    
    /**< If timer was set or re-set, starts or re-starts timer ticking. */

    
/*--------------------------------------------------------------------------*/

 void BenBound::DCRstopTime( void )
  {
    if ( timer == 0 )
     {
       /*cout<<"CPX_LAG::DCRstopTimer(): failed to stop timer";*/
       exit(1);
     }

    timer->Stop();
  }
    
    /**< If timer was set, stops timer ticking. */

/*--------------------------------------------------------------------------*/

 void BenBound::DCRsetTimeLimit( long secs )
  {
    tlimit = secs;
  } 
  
   /**< Sets a time limit for the solver.
    * \param secs long expressing timelimit in seconds */ 

/*--------------------------------------------------------------------------*/
/*----------------------------------SOLVE-----------------------------------*/
/*--------------------------------------------------------------------------*/

   void BenBound::Solve()
   {
    int i,j;
    int Qsize;
    double approx;
    int p = 1; //posizione di inserimento nella nuova suddivisione in Q
    int nonConvexflag; //flag per distinguere il caso di tagli convessi
    int triedflag = 0; //flag per capire se abbiamo già provato il salto sull'ottimo
    double lambda;
    int noLSneeded = 0;

    DCRLagrangianSolver::LAGStatus lgstat;
    DCRLagrangianSolver::LinearCut lc; //left cut
    DCRLagrangianSolver::LinearCut rc; //right cut

    double HeurApprox;

    double sqr;
    double barl;
    double SPcost;
    double releps;
    double zetam;
    double zetap;
    double zetaunico;
    double d;
    double tildez;
    double discrim; //se <=0 quell'arco dà contributo nullo
    double ms; //appoggio per i conti, valore parziale della pendenza

    Inizial();

   // cout<<"inizializzato"<<endl;
    Qsize = Q.size(); // = 2, valori limite.

    /*cout<<"limiti: "<<Q[0].rmin<<"-"<<Q[1].rmin<<endl;*/

    //before next line was commented 
    //solvedflag = 0;

    if(solvedflag == 0)
     {
      if(BenStat == 0)
       {
         do
          {
    //         for(int count_feas = Q[0].rmin; count_feas < Q[1].rmin; count_feas += 1000){
    //           lagSol.updrmin(count_feas);
    //           cout<<" : "<<lagSol.isFeasible()<<endl;
    //           }

    //         cout<<"num of iteration"<< counter_ite_Ben<<endl;
    // cout<<"lagstat = "<<lgstat<<" lambda = " <<lambda<<" Q[p]inter = "<< Q[p].inter<< "valore = "<<Q[p].Val<<endl;
    // cout<<"Objval = "<<ObjVal<<" Heurval = "<<HeurVal<<endl;

            if(Q[p].solflag == 0)
             { 
               if(Q[p].infeasflag == 1) //aggiustiamo le cose per il lagrangiano trovando un nuovo punto di ammissibilità
                { 
                  //cout<<"entro in infeasflag ==1"<<endl;
                   double temppoint; //restringiamo l'intervallo con l'ultimo valore non ammissibile come limite inferiore
                   double oldinterx;
                   int feas = 1;
                   
                   oldinterx = Q[1].inter;

                   //cerchiamo l'altro punto di ammissibilità
                   
                   while(feas != 0)
                    {
                      temppoint = mystep * Q[0].rmin + (1 - mystep) * Q[1].inter; //vediamo chi vari in base alla feasibility

                      lagSol.updrmin(temppoint); 
                      feas = lagSol.isFeasible();

                      if (feas != 0) Q[0].rmin = temppoint;
                      else Q[1].inter = temppoint;
                    }
                     
                    Q[1].interVal = Q[1].Cuts[Q[1].bestCutpos].q + Q[1].inter * Q[1].Cuts[Q[1].bestCutpos].m;//FIXME: in realtà non è sicuro che sia bestCutpos la pos ottima
                    Q[1].solflag = 0;
                    Q[1].infeasflag = 0;

                    //se non mi muovo più mi accontento.

                    if(Q[1].inter < 1) releps = 1;
                    else releps = Q[1].inter;

                    //before next two lines were commented
                    if(abs(Q[1].inter - oldinterx) < eps * releps/100)
                      noLSneeded = 1;

                }//se infeasflag = 1

              lagSol.updrmin(Q[p].inter); //il candidato ci mostri le sue potenzialità :)
              lagSol.Solve();
              counter_ite_Lag += lagSol.getNumIte();

              lgstat = lagSol.getStatus();
              lambda = lagSol.getLambda();

              HeurApprox = lagSol.getHeurVal();

              if(HeurApprox < HeurVal) //se abbiamo trovato una soluzione ammissibile migliore, la salviamo
                HeurVal = HeurApprox;

              SPLabels = lagSol.getSPLabels();

              if(lgstat == 0) //dovrebbe essere garantito da limitrmin e dal check poco sopra la riottimizzazione
              {
                
                lc.m = 0;
                lc.q = 0;

               ///////CALCOLO NUOVI TAGLI

               ///Taglio sinistro
               
               //calcolo del valore tenendo conto del contributo di ciascun arco
               //anche le intercette vengono sommate di arco in arco, visto che alcune causano discontinuità.
                  
                sqr = sqrt(lambda * MTU / Links[1].cost);

                for(i = 0; i < numLinks; i++)
                 {
                  if(Links[i].capacity > Q[0].rmin)
                  {  
                    sqr = sqrt(lambda * MTU / Links[i].cost); //sqrt(lambda * L / f_ij)

                    //definiamo le duali di SP che sono a +\infty
                    if(SPLabels[Links[i].endnode] > 1e200) SPLabels[Links[i].endnode] = 0;
                    if(SPLabels[Links[i].startnode] > 1e200) SPLabels[Links[i].startnode] = 0;

                    d = -SPLabels[Links[i].endnode] + SPLabels[Links[i].startnode]; //d_i-d_j
                    
                    if(d < -1e-10)//altrimenti contributo sicuramente nullo                
                    {
                      if(lambda == 0) //primo caso semplificato
                      {
                        zetaunico = (-d) / Links[i].cost;  //(d_j - d_i)/f_ij

                        if(zetaunico > Q[0].rmin) //se questo non vale avremo contributo nullo
                        {
                          if(zetaunico <= Links[i].capacity)
                          {
                            if(Q[p].inter -zetaunico <= eps*Q[p].inter)  ms = Links[i].cost;  //caso già convesso
                            else  ms = (Links[i].cost * Q[0].rmin + d)/(Q[0].rmin - Q[p].inter);  //facciamo partire la retta approssimante da un punto noto garantito sopra a 0, ovvero Q[0].rmin
                            
                            lc.m += ms;
                            lc.q += - ms * Q[p].inter;
                          }
                          else
                          {
                            if(Q[1].rmin <= Links[i].capacity)  //caso rmin<c_ij<zeta (non dovrebbe succedere)
                              {ms = Links[i].cost;} //caso già convesso
                            else //caso con rmin > c_ij, e zeta > c_ij, qui c'è discontinuità! Aggiriamola //FIXME: non dovrebbe essere la retta tra (c,\phi(c)) ed (r_min,\phi(r_min))?
                              {ms = (Links[i].cost * Q[0].rmin + d)/(Q[0].rmin - Q[p].inter);}
                            lc.m += ms;
                            lc.q += - ms * Q[p].inter;
                            
                          } //zetaunico > c_ij
                        } //zetaunico > Q[0].rmin
                     } //lambda = 0
                    else //caso generale
                     {
                      barl = MTU / Links[i].speed + Links[i].delay + Nodes[Links[i].startnode].delay;
                      discrim = pow(lambda * barl  + d , 2) - 4*Links[i].cost * lambda * MTU; 
              
                      if(discrim > 1e-10) //allora sono definiti zeta_\pm
                        {
                         zetap = (-lambda * barl - d + sqrt(discrim))/(2 * Links[i].cost); //radice più grande
                          
                         if(zetap > sqr) //altrimenti sono gli zeri del massimo
                         {
                            if(sqr < Links[i].capacity)
                            { 

                              if(sqr < Q[p].inter)
                              {
                                ms = (Links[i].cost * sqr + lambda * MTU / sqr + lambda * barl + d)/(sqr - Q[p].inter); //la pendenza è quella della scorciatoia
                                lc.m += ms;

                                if(ms >= Links[i].cost) //cerchiamo il punto di tangenza del nuovo taglio per definirne l'intercetta
                                {

                                  if(zetap < Links[i].capacity)
                                  tildez = zetap;

                                  else tildez = Links[i].capacity;
    
                                }
                                else  tildez = sqrt((lambda * MTU)/(Links[i].cost - ms));
                      
                                lc.q += Links[i].cost * tildez + lambda * MTU / tildez + lambda * barl + d  - tildez * ms;
                              }
                              
                              //cout << zetap << endl;
                            }
                          else //quindi se sqrt >= c_ij
                           {
                            if(Links[i].capacity < Q[p].inter)
                            {
                              //if(Links[i].capacity - Q[p].inter < - 1e-6*Q[p].inter) // nell'ordine giusto, ma non troppo vicini
                             // {
                               zetam = (-lambda * barl - d - sqrt(discrim))/(2 * Links[i].cost); //radice più piccola
                    
                               if(zetam < Links[i].capacity && zetam > Q[0].rmin) //FIXME: perché mi chiedo se zeta- > Q[0].rmin??
                                { 
                                   // cout<<"sqrt = "<<sqr<<" c_ij = "<<Links[i].capacity<<" rmin = "<<Q[p].inter<<" zetamm = "<<zetam<<endl;
                                  ms = (Links[i].cost * Links[i].capacity + lambda * MTU / Links[i].capacity + lambda * barl + d)/(Links[i].capacity - Q[p].inter);
                                //cout<<"ms senza disc "<<ms<<endl;
                                  lc.m += ms;
                                  lc.q += - ms * Q[p].inter;
                                  
                                  //cout << zetam << endl;
                                } //se questo non succede gli archi avranno di nuovo contributo nullo

                             }//c_ij < r_min
                            
                            }//sqrt >= c_ij
                         }//se sono gli zeri del minimo.

                        }//se sono definiti i punti in cui si annulla il contributo dell'arco
                         //altrimenti il contributo sarà nullo.
                      }//caso lambda > 0

                   }//se d < 0

                  }//solo gli archi con capacità più grandi della minima ammissibile!
                  }//per tutti gli archi.

                //infine dobbiamo sommare i contributi della funzione g(z).

                ms = - lambda * Flow.burst / pow(Q[p].inter,2);
                lc.m += ms;
                lc.q += lambda * Flow.burst / Q[p].inter - lambda * Flow.deadline + SPLabels[Flow.sinknode] - Q[p].inter * ms; //g(r_min) - r_min * ms
                  
                //taglio destro

                rc.m = - lambda * Flow.burst / pow(Q[p].inter,2);
                rc.q = lambda * Flow.burst / Q[p].inter - lambda * Flow.deadline + SPLabels[Flow.sinknode] - Q[p].inter * rc.m;

                // cout<< " i tagli prodotti: lc.m = "<<lc.m<<" lc.q = "<<lc.q<<endl;
                // cout<<" rc.m = "<<rc.m<<" rc.q = "<<rc.q<<endl;

                Q[p].Val = rc.q + Q[p].inter * rc.m; //prendiamo il valore nel punto. 
                
                Q[p].solflag = 1;

                ObjVal = Q[p].Val; 

                if(Q[p].Val < BestUB)
                  BestUB = Q[p].Val;

                /*<la convessità vi è quando lc.m <= rc.m, ora rc.m è negativo e lc.m è ottenuto da lui sommando cose positive.
                   Quindi in effetti l'unico caso in cui vi può essere convessità è quando sono uguali.*/
                
                if(lc.m == rc.m && lc.q == rc.q) nonConvexflag = 1; //FIXME:uguaglianza a meno di approssimazioni
                else nonConvexflag = 0;

                /////////CASO NON CONVESSO

                if(nonConvexflag == 0)//aggiungiamo il nuovo punto in Q e i tagli da questo prodotti
                 {
                   BenBound::SubInterval I;
          
                   I.Cuts.resize(Q[p].Cuts.size());

                   for(i = 0; i < Q[p].Cuts.size(); i++)
                     {
                      I.Cuts[i].m = Q[p].Cuts[i].m; //inizializziamo con i tagli del sottointervallo in cui era contenuto
                      I.Cuts[i].q = Q[p].Cuts[i].q;
                     }
         
                   I.rmin = Q[p].inter; 
                   I.Val = Q[p].Val;
                   I.inter = Q[p].inter; //ad ora la soluzione è proprio questa!
                   I.interVal = Q[p].interVal;
                   I.pCut.m = Q[p].pCut.m; //ci copiamo anche i tagli che la definiscono
                   I.pCut.q = Q[p].pCut.q; //da cui far partire la LS.
                   I.mCut.m = Q[p].mCut.m;
                   I.mCut.q = Q[p].mCut.q;
                   I.solflag = 1;
                   I.branchedflag = 0;
                   I.infeasflag = 0;

                   if(Q[p].Cuts.size() == 1) 
                    {
                      Q[p].mCut.q = rc.q; 
                      Q[p].mCut.m = rc.m;

                      I.mCut.q = lc.q;
                      I.mCut.m = lc.m;
                    } //se eravamo nel caso definito da un solo taglio, adesso abbiamo anche l'altro
                 
                   Q.insert(Q.begin()+p, I);
                   Qsize++;

                  //aggiungiamo lc a tutti i punti a sx:
 
                  for(i = 1; i <= p; i++) //non in pos 0, visto che non rappresenta davvero un subinterval
                    {
                      Q[i].Cuts.push_back(lc);
                    }

                  //aggiungiamo rc a tutti i punti a dx:
       
                  for(i = p+1; i < Qsize; i++)
                   {
                     Q[i].Cuts.push_back(rc);
                   }

                  } //if non convex.
        
                  //////CASO CONVESSO non si aggiunge un nuovo punto in Q, solo i tagli prodotti in ogni posizione preesistente
                else
                  { 
                    /*cout<<"caso convesso"<<endl;*/
                    is_convex_iteration = true;
                     
                    if(Q[p].Cuts.size() == 1) 
                    {
                      Q[p].mCut.q = rc.q; 
                      Q[p].mCut.m = rc.m;
                    } //se eravamo nel caso definito da un solo taglio, adesso abbiamo anche l'altro

                    for(i = 1; i < Qsize; i++)
                     {
                        Q[i].Cuts.push_back(lc); //anche in posizione p, non avendolo aggiunto prima
                     }
                  } // if convex
              }//se il lagrangiano trova una soluzione        
              else //altrimenti vorrà dire che stiamo sbagliando qualcosa
               { 
                 /*cout<<"Errore in BenBound::Solve"<<endl;*/
                 exit(1);
               }
            }//se non avevamo già risolto per quel valore

           //next line originally not in the code
           ///!if(Q[p].solflag == 1){noLSneeded = 1;}

           if(noLSneeded == 0) //if(noLSneeded == 0)
           {
            p = LineSearch();
            //if(Q.size()==2||Q.size()==3) cout<<"P = "<<p<<endl;
           }

           if(Q[p].interVal > BestLB)
             BestLB = Q[p].interVal;

           ObjVal = Q[p].Val; 

           approx = Q[p].interVal; //miglioriamo le nostre stime con i valori dati dalla LS

           if(approx > 1) releps = approx; //precisione relativa
           else releps = 1;

           counter_ite_Ben++; //numero di punti visitati
           
           //if(Q.size()==2||Q.size()==3) cout<<"numero di punti = "<<Q.size()<<" with rmin = "<<rmin<<endl;
           
          //}while(solvedflag == 0 && abs(ObjVal - approx) > eps * releps); //(before on the code)
          }while(counter_ite_Ben <= 100 && solvedflag == 0 && (abs(ObjVal - approx) > eps * releps/100));

         ObjVal = approx; //in ogni caso dobbiamo restituire un LB!

         rmin = Q[p].inter;
         //cout << "Q[p].Val = " << Q[p].Val << endl;
         //cout << "rmin = " << rmin << endl;

       }//se non ci sono errori 
     }//se non era già risolto il problema
     
     SOLUTION = lagSol.getRSOLS();
   }

/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

   BenBound::~BenBound()
   {
    clean_up();
   }

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/


   void BenBound::copyDataArray(DCR::DCRFlow flow, DCR::DCRLink *links, DCR::DCRNode *nodes)
   {

    int i;
    Links = new DCR::DCRLink[numLinks];
    Nodes = new DCR::DCRNode[numNodes];

    Flow = flow;

    for(i = 0; i < numLinks; i++)
    {
      Links[i].startnode = links[i].startnode;
      Links[i].endnode = links[i].endnode;
      Links[i].speed = links[i].speed;
      Links[i].capacity = links[i].capacity;
      Links[i].delay = links[i].delay;
      Links[i].cost = links[i].cost;//Flow.costs[i];
    }
   
    for(i = 0; i < numNodes; i++)
    {
      Nodes[i].delay = nodes[i].delay;
    }
   }


/*--------------------------------------------------------------------------*/

   int BenBound::LineSearch()
   {
     int i, j ;
     double interx;
     double interxApprox;
     double interxVal;
     double someinter;
     int maxpos;
     double max; //in realtà poi sarà il minimo della funzione massimo
     double min;
     int minpos;
     double releps;
     int counter = 0;

     int corrflag = 0;
     int minconflag = 0;
     int isize;
     int branchcounter = 0;

     for(i = 1; i < Q.size(); i++) //su ogni sottointervallo a parte quello in posizione 0, facciamo una LS
       { 
         //poscounter = 0;
         isize = Q[i].Cuts.size();
         max = -Inf<double>(); //per sicurezza lo rinizializziamo.
         counter = 0;

             if(Q[i].branchedflag == 0)
              {
               //facciamo un po' di controlli per capire se siamo nel caso con tutte rette positive:
               if(Q[i].pCut.m > 0 && Q[i].mCut.m > 0 && Q[i].Cuts[isize-1].m > 0) //allora sono tutti positivi
               {
                 maxpos = whchbest(i);
                 interx = Q[i-1].rmin;
                 interxVal = Q[i].Cuts[maxpos].q + interx * Q[i].Cuts[maxpos].m;
                 if(abs(Q[i-1].inter - Q[i-1].rmin) > 0) //se non era il suo valore ottimo si inizializza
                      {
                        Q[i].solflag = 0;
                        Q[i].Val = Inf<double>();
                      }
                      
                 if(i == 1) //altrimenti ci serve di conservare il suo valore
                   { 

                     Q[i].bestCutpos = maxpos;
                     Q[i].infeasflag = 1;
                     Q[i].interVal = Q[i].pCut.q + Q[i-1].rmin * Q[i].pCut.m;
                   }
                  else
                   {
                     Q[i].inter = Q[i-1].rmin;
                     Q[i].interVal = Q[i].pCut.q + Q[i-1].rmin * Q[i].pCut.m;
                   }
               } 
               else
               {
                if(Q[i].pCut.m > 0 && Q[i].mCut.m > 0 && Q[i].Cuts[isize-1].m <= 0) 
                  UpdCut(Q[i].Cuts[isize-1].q, Q[i].Cuts[isize-1].m,i); //se l'ultimo è negativo aggiorniamo
                
                do //finalmente la LS! //poi in ogni caso partiamo con la LS
                 {
                  //poscounter = 0;
                  interx = (Q[i].pCut.q - Q[i].mCut.q) / (Q[i].mCut.m - Q[i].pCut.m); //calcoliamo la nuova intersezione
      
                  interxApprox = Q[i].mCut.q + interx * Q[i].mCut.m; //e l'altezza sui vecchi tagli
                  //double test = Q[i].pCut.q + interx * Q[i].pCut.m;
    
                  for(j = 0; j < isize; j++) //cerchiamo il taglio che garantisce l'altezza massima in quel punto
                   {
                    someinter =  Q[i].Cuts[j].q + interx * Q[i].Cuts[j].m;

                    if(max <= someinter)
                      {
                        max = someinter;
                        maxpos = j;
                      }            
                   }
   
                  interxVal = Q[i].Cuts[maxpos].q + interx * Q[i].Cuts[maxpos].m; //trovando il vero valore dell'approssimazione
                  
                  UpdCut(Q[i].Cuts[maxpos].q, Q[i].Cuts[maxpos].m, i); //aggiorniamo quindi uno dei tagli che definisce la soluzione

                  counter++;

                  if(abs(interxVal) > 1) releps = abs(interxVal);
                  else releps = 1;
                     
                 //}while(abs(interxVal - interxApprox) > 1e-3);
                 }while(abs(interxVal - interxApprox) > eps*releps/100 && counter<=100); //tanto sono LS esatte, volendo si può aggiungere un releps
              
             //a questo punto dobbiamo controllare di essere rimasti all'interno del sottointervallo,
             //altrimenti prenderemo come valore l'estremo più vicino.

                if(interx <= Q[i].rmin && interx > Q[i-1].rmin) //se siamo dentro l'ottimo è quello trovato
                 {
                   if(abs(Q[i].inter - interx) >  eps * interx/100) //nuovo valore, rinizializziamo tutto.
                    {                   
                      Q[i].solflag = 0;
                      Q[i].Val = Inf<double>();
                    }

                   Q[i].inter = interx;
                   Q[i].interVal = interxVal;

                 }
                else
                 {
                   if(interx <= Q[i-1].rmin) //se sto nell'intervallo precedente
                   { 
                      if(abs(Q[i-1].inter - Q[i-1].rmin) > eps * Q[i-1].rmin/100) //se non era il suo valore ottimo si inizializza
                      {
                        Q[i].solflag = 0;
                        Q[i].Val = Inf<double>();
                      }
                      
                      if(i == 1) //altrimenti ci serve di conservare il suo valore
                        { 
                          Q[i].bestCutpos = maxpos;
                          Q[i].infeasflag = 1;
                          Q[i].interVal = Q[i].pCut.q + Q[i-1].rmin * Q[i].pCut.m;
                        }
                      else
                        {
                          Q[i].inter = Q[i-1].rmin;
                          Q[i].interVal = Q[i].pCut.q + Q[i-1].rmin * Q[i].pCut.m;
                        }
                   }
          
                   else //se sto nell'intervallo successivo
                   {
                     if(abs(Q[i].inter - Q[i].rmin) > eps * Q[i].rmin/100) //se non era il suo valore ottimo si inizializza
                      {
                        Q[i].solflag = 0;
                        Q[i].Val = Inf<double>();
                      }

                      Q[i].inter = Q[i].rmin;
                      Q[i].interVal = Q[i].mCut.q + Q[i].inter * Q[i].mCut.m;
                   }
     
                  }
                }//se non erano tutti positivi
             }//se non era già escluso che l'ottimo fosse in questo intervallo
       }//per tutti i sottointervalli
       
    //svolte tutte le LS, dobbiamo capire chi sia il minimo tra i minimi e restituire la sua posizione
    // lui sarà il valore più promettente, per cui calcoleremo la soluzione del DL.
    // Potrebbe essere fatto all'interno, per ora messo qui per evitare errori.
    // Mentre scorriamo possiamo anche fare del branching: se per un certo sottointervallo
    // il valore trovato dalla LS è >= di ObjVal, quindi di quello ritenuto ad ora l'ottmo,
    // potremo eliminare quel sottointervallo

       min = Inf<double>();

       for(i = 1; i < Q.size(); i++)
       {
         if(min >= Q[i].interVal)
         {
           min = Q[i].interVal;
           minpos = i;
         }
         
         if(Q[i].interVal >= BestUB)
            {
              Q[i].branchedflag = 1;
              branchcounter++;
            }
       }

       //next line not in the original code
       //if(Q[minpos].interVal < BestLB)
        BestLB = Q[minpos].interVal;
       ApproxVal = Q[minpos].interVal;

       if(corrflag == 0) return minpos;
       else return minconflag;
   }


/*--------------------------------------------------------------------------*/  

   void BenBound::UpdCut(double alpha, double beta, int i)
   {
      if(Q[i].pCut.m != beta || Q[i].pCut.q != alpha) 
      {
         if(Q[i].mCut.m != beta || Q[i].mCut.q != alpha) //se era già uguale a uno dei due, è inutile aggiornarlo
         {
          if(beta >= 0)
           {
             Q[i].pCut.q = alpha;
             Q[i].pCut.m = beta;
           }
          else
           {
             Q[i].mCut.q = alpha;
             Q[i].mCut.m = beta;
           }
         }
      }
   }
/*<modifica uno dei due tagli ottimi che definiscono la soluzione al momento a seconda della pendenza*/

/*--------------------------------------------------------------------------*/   

   int BenBound::whchbest(int i)
   {

      int j;
      double max = -Inf<double>();
      int maxpos;
      double interc;

      for(j = 0; j < Q[i].Cuts.size(); j++)
       {
        interc = Q[i].Cuts[j].q + Q[i-1].rmin * Q[i].Cuts[j].m;

        if(interc >= max)
         {
           max = interc;
           maxpos = j;
         }
       }

      //cout<<"questo è il return della LineSearch "<<maxpos<<endl;

     return maxpos;
   }

/*--------------------------------------------------------------------------*/   

   
  void BenBound::Inizial()
  {
  	int i;
    double* limits;
    double lambda;
    double sqr;
    double barl; //\bar{l}_ij
    double SPcost; //lambda barl + f_ij rmin + lambda L /rmin
    int feas = 1;
    int in0opt;
    double zetam;
    double zetap;
    double zetaunico;
    double d;
    double tildez;
    double discrim; //se <= 0 quell'arco dà contributo nullo
    double ms; //appoggio per i conti, valore parziale della pendenza
    double counter = 0;
    double pointeps = 1e-9;

    DCRLagrangianSolver::LAGStatus lgstat;
    DCRLagrangianSolver::LinearCut c;
    c.m = 0;
    c.q = 0;

 
    limits = new double[2];

    limits = Limitrmin();//valori limite per cui è garantita l'ammissibilità di r_min
    //if(limits != NULL) {cout<<"valori limite per rmin: "<<limits[0]<<"-"<<limits[1]<<endl;}
    //if(limits == NULL) {cout<<"limit is null"<<endl;}
    
    mystep = 0.9; // before it was 0.9
    /*cout<<"sono qui"<<endl;*/

    if(BenStat == 0)
     {
      if(limits[0] == limits[1])
      {
        
        lagSol.updrmin(limits[1]);
      //  lambda = lagSol.getLambda();
        lagSol.Solve();
        counter_ite_Lag += lagSol.getNumIte();
        lgstat = lagSol.getStatus();
        lambda = lagSol.getLambda();

        SPLabels = lagSol.getSPLabels();
        ObjVal = lambda * Flow.burst / limits[0] - lambda * Flow.deadline + SPLabels[Flow.sinknode];

        rmin = limits[1];

        solvedflag = 1;

        //cout<<"risolto per limitrmin uguali"<<endl;
      }  
      else //caso generale in cui esista un intervallo non banale di ammissibilità
      {

        Q.resize(2);
        
        Q[0].solflag = 0;
        Q[1].solflag = 0;
        Q[1].branchedflag = 0;
        Q[1].infeasflag = 0;


        Q[0].rmin = limits[0];
        Q[1].rmin = limits[1];

       //cerchiamo per prima cosa un estremo sinistro vicino a quello trovato, ma garantito ammissibile, 
        //lo salviamo per ora in Q[1].inter.

         double temppoint; //restringiamo l'intervallo con l'ultimo valore non ammissibile come limite inferiore
         
         Q[1].inter = Q[0].rmin;

         //cerchiamo l'altro punto di ammissibilità

         while(feas != 0)
         {
           temppoint = Q[1].inter;
           Q[1].inter = myparam * Q[1].inter + (1 - myparam) * Q[1].rmin; //l'estremo sinistro varia, il destro resta fisso a Q[1].rmin

           lagSol.updrmin(Q[1].inter); 
           feas = lagSol.isFeasible();
           counter ++;
           //cout<<"cerco più volte?"<<endl;
         }

        Q[0].rmin = temppoint;

        //cout<<"temppoint = "<<Q[0].rmin<<endl;

       //////FASE PRELIMINARE DI RESTRIZIONE DELL'INTERVALLO.

        if(crit_capc < limits[1]) //caso in cui la criticità è nell'intervallo, partiamo da lì!
         {
          
          //cout<<"sono qui in crit_capc"<<endl;
         	lagSol.updrmin(crit_capc);

         	in0opt = lagSol.is0opt();

         	if(in0opt == 0) //in questo caso possiamo già tagliare tutto quello che è alla sua destra
         	 {  
            /*cout<<"sono dentro in0opt crit_cap"<<endl;*/
         	 	double tempd = crit_capc;
            /*cout<<"starting from"<<crit_capc<<endl;*/
         	 	double temps = Q[1].inter;
         	 	double cent = (tempd + temps) / 2;
               
               while(tempd - temps > pointeps * temps)//cerchiamo con una LS il primo punto per cui lambda = 0
               {
                  lagSol.updrmin(cent);

                  in0opt = lagSol.is0opt();

                  if(in0opt == 0) tempd = cent;
                  else temps = cent;

                  cent = (tempd + temps) / 2;
               }
              /*cout<<"Q[1].rmin = "<<tempd<<endl;*/
             Q[1].rmin = tempd;

             //cout<<"Q0rmin-Q1rmin = "<<Q[0].rmin<<"-"<<Q[1].rmin<<endl;

         	 }
         	else //anche se lambda != 0 vogliamo partire appena prima della criticità con l'iterazione.
         	 {
         	   mystep = 0.75; //mystep = 0.75;
             Q[1].inter = crit_capc - 2; //fermiamoci leggermente prima per evitare errori di calcolo

             //cout<<"Q[1].inter  = "<<Q[1].inter<<endl;
         	 }
        } //critcap < Q[1].rmin
       else //altrimenti proviamo comunque a restringere l'intervallo!
        { 
          //cout<<"sono di qua"<<endl;
          lagSol.updrmin(Q[1].rmin);

          in0opt = lagSol.is0opt();

         	if(in0opt == 0) //in questo caso possiamo già tagliare tutto quello che è alla sua destra
         	 {  
            /*cout<<"sono dentro in0opt"<<endl;*/

              //cout<<"entro nel taglio?"<<endl;
              /*cout<<"starting from"<<Q[1].rmin<<endl;*/
              double tempd = Q[1].rmin;
              double temps = Q[1].inter;
              double cent = (tempd + temps) / 2;
               
               while(tempd - temps > pointeps * temps)//cerchiamo con una LS il primo punto per cui lambda = 0
               {
                 lagSol.updrmin(cent);

                 in0opt = lagSol.is0opt();

                 if(in0opt == 0) tempd = cent;
                 else temps = cent;

                 cent = (tempd + temps) / 2;
               }
               /*cout<<"Q[1].rmin = "<<tempd<<endl;*/
 
               Q[1].rmin = tempd;

            } //se non ci riesce fa nulla

         }//critcap = Q[1].rmin
       
       /////CALCOLO DEL TAGLIO SINISTRO NELL'ESTREMO DESTRO DELL'INTERVALLO

    /*<in ogni caso dopo aver ristretto l'intervallo, iniziamo a popolarlo con un taglio estremale*/

        lagSol.updrmin(Q[1].rmin); //taglio sinistro per Q[1]
        lagSol.Solve();

        counter_ite_Lag = lagSol.getNumIte();
        lgstat = lagSol.getStatus();
        lambda = lagSol.getLambda();

        //cout<<"la soluzione del lagra era a lambda = "<<lambda<<endl;
        if(lgstat == 0) //in realtà dovrebbe essere garantito da Limitrmin.
        {

          SPLabels = lagSol.getSPLabels();

          ObjVal = lambda * Flow.burst / Q[1].rmin - lambda * Flow.deadline + SPLabels[Flow.sinknode];
          HeurVal = lagSol.getHeurVal();
        
          //cout<<ObjVal<<HeurVal<<endl;

          //before the next line was not commented
          if(abs(ObjVal - HeurVal) < eps*abs(ObjVal)) {solvedflag = 1;}
          
         //         cout<<"objval = "<<<<" Heurval ="<<lagSol.getHeurVal()<<endl;

        //taglio sinistro
        ///calcolo del valore tenendo conto del contributo di ciascun arco
        //anche le intercette vengono sommate di arco in arco, visto che alcune causano discontinuità.

                for(i = 0; i < numLinks; i++)
                 {
                  
                  if(Links[i].capacity > Q[0].rmin)
                  {
                  if(SPLabels[Links[i].endnode] > 1e200) SPLabels[Links[i].endnode] = 0;
                  if(SPLabels[Links[i].startnode] > 1e200) SPLabels[Links[i].startnode] = 0;

                  sqr = sqrt(lambda * MTU / Links[i].cost); //sqrt(lambda * L / f_ij)
                  d =  SPLabels[Links[i].startnode] - SPLabels[Links[i].endnode]; //d_i-d_j
                  //cout<<"d = "<<d<<endl;
                  if (d < 0)
                  {

                     if(lambda == 0) //primo caso semplificato
                      {
  
                      zetaunico = (-d) / Links[i].cost;  //(d_j - d_i)/f_ij

                      if(zetaunico > Q[0].rmin)
                      {
                       if(zetaunico <= Links[i].capacity)
                        { 
                         if(Q[1].rmin <= zetaunico)  ms = Links[i].cost; //caso già convesso
                         else  ms = (Links[i].cost * Q[0].rmin + d)/(Q[0].rmin - Q[1].rmin); //facciamo partire la retta approssimante da un punto noto garantito sopra a 0, ovvero Q[0].rmin

                         c.m += ms;
                         c.q += - ms * Q[1].rmin;
                        }
                       else
                       {
                        if(Q[1].rmin <= Links[i].capacity)  //caso rmin<c_ij<zeta (non dovrebbe succedere)
                          {ms = Links[i].cost; } //caso già convesso
                        else //caso con rmin > c_ij, e zeta > c_ij, qui c'è discontinuità! Aggiriamola
                           {ms = (Links[i].cost * Q[0].rmin + d)/(Q[0].rmin - Q[1].rmin);}
                        c.m += ms;
                        c.q += - ms * Q[1].rmin;
                          
                       } //zetaunico > c_ij
                       } //zetaunico > Q[0].rmin
                       //cout<<"for i ="<<i<<" c.q = "<<c.q<<"Q[0].rmin - Q[1].rmin = "<<Q[0].rmin - Q[1].rmin<<endl;
                       //if(i>50) exit(1);
                      }//lambda = 0
                     else //caso generale
                     { 
                      barl = MTU / Links[i].speed + Links[i].delay + Nodes[Links[i].startnode].delay;
                      discrim = pow(lambda * barl  + d ,2) - 4*Links[i].cost * lambda * MTU;
              
                     if(discrim > 0) //allora sono definiti zeta_\pm
                        {
                         zetap = (-lambda * barl - d + sqrt(discrim))/(2 * Links[i].cost); //radice più grande
        
                         if(zetap > sqr) //altrimenti parliamo delle radici del massimo
                         {
                          if(sqr < Links[i].capacity)
                           { 
                            if(sqr < Q[1].rmin)
                            {   
                            ms = (Links[i].cost * sqr + lambda * MTU / sqr + lambda * barl + d)/(sqr - Q[1].rmin); //la pendenza è quella della scorciatoia

                            c.m += ms;

                            if(ms >= Links[i].cost) //cerchiamo il punto di tangenza del nuovo taglio per definirne la pendenza
                             {

                               if(zetap < Links[i].capacity)
                               tildez = zetap;

                               else tildez = Links[i].capacity;

                             }
                            else  tildez = sqrt((lambda * MTU)/(Links[i].cost - ms));
                  
                            c.q += Links[i].cost * tildez + lambda * MTU / tildez + lambda * barl + d  - tildez * ms;
                            //if(ms < 0 ) cout<<"ms nel primo if = "<<ms<<" per d = "<<d<<endl;

                            }//sqr < r_min
                           }
                          else //quindi se sqrt >= c_ij ////FIXME: si potrebbe chiedere all'inizio che \phi(\bar{v})<0 altrimenti non fare nulla!
                           {

                            if(Links[i].capacity < Q[1].rmin)
                            {  
                             zetam = (-lambda * barl - d - sqrt(discrim))/(2 * Links[i].cost); //radice più piccola
                            
                             if(zetam < Links[i].capacity && zetam > Q[0].rmin)
                               { 
                                 ms = (Links[i].cost * Links[i].capacity + lambda * MTU / Links[i].capacity + lambda * barl + d)/(Links[i].capacity - Q[1].rmin);
                                 c.m += ms;
                                 c.q += - ms * Q[1].rmin;
                                 
                               } //se questo non succede gli archi avranno di nuovo contributo nullo

                               //if(ms < 0 ) cout<<"ms nel secondo if = "<<ms<<" per d = "<<d<<endl;;

                            } //c_ij <  r_min
                          }
                         }//se le radici si riferiscono al minimo

                        }//se sono definiti i punti in cui si annulla il contributo dell'arco
                         //altrimenti il contributo sarà nullo.
                      }//caso lambda > 0

                    //} //caso sqrt < \bar{r}_min

                   }//se d<0 altrimenti il contributo sarà nullo.
                  }
                  }//per tutti gli archi.

               //infine dobbiamo sommare i contributi della funzione g(z).
                
                ms = - lambda * Flow.burst / pow(Q[1].rmin,2);
                c.m += ms;
                c.q += lambda * Flow.burst / Q[1].rmin - lambda * Flow.deadline + SPLabels[Flow.sinknode] - Q[1].rmin * ms; //g(r_min) - r_min * ms

                // cout<<"Q[1].rmin = "<<Q[1].rmin<<" SpLabels[sinknode] = "<<SPLabels[Flow.sinknode]<<endl;
                // cout<<"ms ="<<ms<<" c.m ="<<c.m<<" c.q ="<<c.q<<endl;
                
        //fine calcolo taglio sinistro, salviamolo.

         Q[1].Cuts.push_back(c);

         Q[1].pCut.q = c.q; ///uno dei tagli che definisce l'attuale punto da cui parte
         Q[1].pCut.m = c.m; // la LS
       }
       else{ /*cout<<"Errore in inizial di BenBound"<<endl;*/ exit(1);}

       if(Q[1].pCut.m <= 0)//se la pendenza all'estremo è già negativa sono già sull'ottimo!
        {
          //cout<<"sono nel caso pcut.m<0?"<<endl;
          Q[1].solflag = 1;
          
          solvedflag = 0;

          //cout<<"risolto per pendenza all'estremo negativa"<<endl;

          ObjVal = lambda * Flow.burst / Q[1].rmin - lambda * Flow.deadline + SPLabels[Flow.sinknode];

          ApproxVal = ObjVal;

          rmin = Q[1].rmin;

          HeurVal = lagSol.getHeurVal();

          //BestLB e BestUB nel caso solvedflag = 1
          ////BestLB = ObjVal;
          ////BestUB = HeurVal;
        }
          //altrimenti inizializiamo Q[1].interval con un il suo valore e passiamo il tutto al solve
       
        Q[1].mCut.q = 0; 
        Q[1].mCut.m = 0; //diamo un'inizializzazione riconoscibile

       }//se esisteva un intervallo non banale di ammissibilità
     }//se Benstat = 0
     //cout<<"solved flag = "<<solvedflag<<endl;
     //exit(1);
  }



/*--------------------------------------------------------------------------*/


   double* BenBound::Limitrmin()
   {
    int i, j = 0; 
    int dx, sx, cent;
    int feas; //feasibility check
    //double rp, rm;
    
    sx = 0;
    dx = numLinks-1;
    cent = (dx + sx)/2;
    caps = new double[numLinks];

    for(i = 0; i < numLinks; i++)
      caps[i] = Links[i].capacity;
    
    capSort(sx, dx);

    ////DEBUG_DATA
    // double* red_caps = new double[numLinks];
    
    // red_caps[0] = caps[0];
    // j = 0;

    // for(i = 0; i < numLinks; i++){

    //   if(caps[i] != red_caps[j]){
    //     red_caps[j+1] = caps[i];
    //     j++;
    //   }
    // }
    // cout<<"types of caps:";
    // for(i = 0; i < j+1; i++) cout<<", "<<red_caps[i];    
    
    // cout<<";"<<endl;
    // j = 0;
    //END_DEBUG_DATA

    ///Non è vero che se per la minima capacità non c'è soluzione,
    //allora necessariamente il problema è vuoto (Waxman100_4).
    //quindi prendiamo come estremo sinistro il valore max{Flow.rate, Flow.burst/Flow.deadline}

    lagSol.LoadProblem(numNodes, numLinks, Flow, Links, Nodes, MTU, caps[0]);
    feas = lagSol.isFeasible();

/*
    if(feas == 1){
      lagSol.updrmin(std::max(Flow.rate, Flow.burst/Flow.deadline));
      feas = lagSol.isFeasible();
    }
*/

    while(feas != 0 && j < numLinks)
      { 
        j++;

        if(caps[j] != caps[j-1]) //stesso rmin implica stesso risultato
        {
          lagSol.updrmin(caps[j]);
          feas = lagSol.isFeasible();
        }
      }

    crit_capc = caps[j];

    if (feas == 0) 
      { 
        BenStat = OK;
        double* bounds = new double[2];
        
        if(Flow.rate < Flow.burst / Flow.deadline) bounds[0] = Flow.burst / Flow.deadline; //in ogni caso questo è il lower bound per r_min;
        else bounds[0] = Flow.rate;

        lagSol.updrmin(caps[dx]);
        feas = lagSol.isFeasible();   

        if (feas == 0) //se per l'ultima capacità va bene abbiamo il range;
         {
           bounds[1] = caps[dx];
           return bounds;
         }
        else //altrimenti ricerca binaria sul vettore delle capacità per trovare quella minima
         {
          sx = j; //a partire dal lower bound;

          //FIXME: si potrebbe migliorare, osservando che le capacità spesso si ripetono
          //e facendo un check su sx+1, visto che spesso esiste solo un valore ammissibile.
          while( dx - sx > 1 )
           { 
             cent = (sx + dx)/2;
             lagSol.updrmin(caps[cent]);
             feas = lagSol.isFeasible();

             if (feas == 0) sx = cent;
             else dx = cent;
           }

          lagSol.updrmin(caps[dx]);
          feas = lagSol.isFeasible();
          if(feas == 0) 
          {bounds[1] = caps[dx];}
        /*lagSol.updrmin(caps[dx]+1);
        feas = lagSol.isFeasible();
        if(feas==0)
        {
          cout<<"limitrmin"<<endl;
          exit(1);
        } }*/
 

          else {bounds[1] = caps[sx];}

        /* lagSol.updrmin(caps[sx]+1);
          feas = lagSol.isFeasible();
          if(feas==0)
          {
            cout<<"sx limitrmin"<<endl;
            exit(1);
          } }*///resta il dubbio che in questi casi ci potrebbe essere del margine

          /**FIXME: se SPT si sconnette per r_min uguale a una certa capacità,
          è chiaro che non si possa più salire, ma se il lagrangiano diventa semplicemente
          infeasible è altrettanto chiaro che, nel grande gap tra l'ultima capacità ammissibile
          e la prima inammissibile, non ci siano valori salvabili?*/

          return bounds;
          
         } 
      } //se il problema è ammissibile
    
    else
     {
       BenStat = Infeasible;
       return NULL;
     }

   }

/*--------------------------------------------------------------------------*/

   void BenBound::capSort(int sx, int dx) //Quicksort ricorsivo randomizzato.
   {
     int pivot, rango;
     srand((unsigned)time(NULL));

     if(sx < dx)
     { 
      pivot = rand()%(dx-sx)+sx;
 
      rango = Distrib(sx,pivot,dx);
      capSort(sx,rango-1);
      capSort(rango+1,dx);
     }
   }

/*--------------------------------------------------------------------------*/

   int BenBound::Distrib(int sx, int pv, int dx)  //metodo privato per quicksort 
   {                                    
    int i,j;

    if(pv != dx)   Swap(pv,dx);

    i = sx;
    j = dx-1;
 
    while(i <= j)
     {
      while(i<=j && caps[i]<=caps[dx])
        i++;

      while(i<=j && caps[j]>=caps[dx])
        j--;

      if(i<j)  Swap(i,j);
        
     }

    if(i != dx) Swap(i,dx);
   
   return i;
  }

/*--------------------------------------------------------------------------*/

  void BenBound::Swap(int a,int b)   //metodo privato per quicksort
  {                           
   double temp;
 
   temp = caps[a]; 
   caps[a] = caps[b];
   caps[b] = temp;
  }

/*--------------------------------------------------------------------------*/

  void BenBound::clean_up()
  {
    if(Links != 0) {delete[] Links; Links = 0;}
    if(Nodes != 0) {delete[] Nodes; Nodes = 0;}
    if(caps != 0) {delete[] caps; caps = 0;}

    for(int i = 0; i < Q.size(); i++)
      Q[i].Cuts.clear();
    Q.clear();

    xPlot.clear();

    SOLUTION.clear();
  }



/*--------------------------------------------------------------------------*/
/*-------------------------- End File BenBound.cpp -------------------------*/
/*--------------------------------------------------------------------------*/
