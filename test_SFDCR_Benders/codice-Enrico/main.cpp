#include "DCR.h"
#include "DCRGraph.h"
#include "DCRLagrangianSolver.h"
#include "SPT.h"
#include "CPX_LAG.h"
#include "BenBound.h" 
//#include "Benders.h"
#include "CPL_BEND.h" 
//#include "DCR_MFSP_SOCP_GRB.h"
//#include "DCR_SPT.h"
//#include "DCR_INDI.h" 
//#include "DCR_3PRONG.h"
#include <math.h>
#include <iostream>
#include <random>
#include <fstream>


void ReadData(char *filename, char format, int * nnodes, int *nlinks, int *nflows, double *MTU,
DCR::DCRFlow **flows, DCR::DCRLink **links, DCR::DCRNode **nodes, double perccap, double percdelay, int seed);

/*void MyRead(char *filename, char format, int * nnodes, int *nlinks, int *nflows, double *MTU,
DCR::DCRFlow **flows, DCR::DCRLink **links, DCR::DCRNode **nodes, double perccap, double percdelay, int seed);*/



int main (int argc, char ** argv)
{ 
  int i;
  char *filein;
  char format = 'm';  //formato Mnetgen
  filein=argv[1];
  int nnodes;
  int nlinks;
  int nflows;
  int attflow = 0; //flusso che stiamo servendo, ricordarsi di aggiornare
  double MTU;
  double ObjVal;
  double ApproxSol;
  double cpxObjVal;
  double lagObjVal;
  double contBound;
  double intBound;
  double imprcontBound;
  double lambda;
  double Heuristic;
  double rmin;
  double myrmin;
  double cplrmin;
  double beta;
  bool is = true;

  
  //ofstream fout("debuggw5.txt", ios::app);
  //fout<<"----------------------------------------------"<<endl;


  DCR::DCRFlow *flows = 0;
  DCR::DCRLink *links = 0;
  DCR::DCRNode *nodes = 0;

  BenBound::BndrStat stat;
  CPX_LAG::CPLStatus xstat;
  CPL_BEND::CPLBStatus contstat;
  CPL_BEND::CPLBStatus imprcontstat; 
  CPL_BEND::CPLBStatus intstat; 
  DCRLagrangianSolver::LAGStatus lagstat;
  //CPL_BEND::CPLBStatus BenStat;
  cout<<"reading data"<<endl;
  ReadData(filein, format, &(nnodes), &(nlinks), &(nflows), &(MTU), &(flows), &(links), &(nodes), 0, 0, 1);
  cout<<"done!"<<endl;
  
  BenBound BenSol;
  CPX_LAG rfissatoSOL;
  DCRLagrangianSolver lagSol;
  CPL_BEND contSol;
  CPL_BEND imprContSol;
  CPL_BEND intSol;

  double contdurat = 0;
  double intdurat = 0;
  double mydurat = 0;
  double imprContdurat = 0;

  int ntests = nflows;
  int j;
  default_random_engine generator;
  uniform_real_distribution<double> distribution(1,100);

  int integral = 0;
  int rlowrm = 0;
  int rbigc = 0;
  double delay = 0;
  double fo = 0;
  int counter = 0;


 /*lagSol.LoadProblem(nnodes, nlinks, flows[1], links, nodes, MTU,13402.8);
 lagSol.Solve();

 rfissatoSOL.LoadProblem(nnodes, nlinks, flows[1], links, nodes, MTU,13402.8);
 xstat = rfissatoSOL.Solve();

 double mine = lagSol.getObjVal();
 double cplex = rfissatoSOL.getObjVal();

 cout<<"nel punto incriminato il lagrangiano ha una precisione di "<<(mine-cplex)/cplex<<endl;*/

 /* ////PROFILING //aggiungere -pg nel makefile

  cout<<"main per profiling"<<endl;

 for(i=0; i<5; i++)
 {
  BenSol.LoadProblem(nnodes, nlinks, flows[i], links, nodes, MTU);
  BenSol.Solve();
  cout<<"servito il flusso "<<i<<endl;
 }*/

/*//TEST per Sperimentazione, crea un file in cui salviamo:
  //efficienza di Benders vs cplex
  //precisione nei casi non venga risolta all'ottimo
  //ovvero:
  //HeurVal - Benval su Benval
    double HB = 0;
  //HeurVal - intval su intval
    double HI = 0;
  //Intval - Benval su Benval
    double IB = 0;
  //Benval- contval1 su contval1
    double BCl = 0;
  //Benval- contval2 su contval2
    double BCu = 0;
  //Ben time | int time | conttime 2 | conttime 1
cout<<"finita l'inizializzazione, iniziamo i test"<<endl;
for(i = 11; i < 12; i ++)
  {
   intSol.DCRsetTime(is);
  // contSol.DCRsetTime(is);
   BenSol.DCRsetTime(is);
  // imprContSol.DCRsetTime(is);

   BenSol.LoadProblem(nnodes, nlinks, flows[i], links, nodes, MTU);
   cout<<"caricato il problema su Benders"<<endl;
   BenSol.DCRstartTime();
   BenSol.Solve();
   BenSol.DCRstopTime();
   mydurat = BenSol.getTime(); 
  // cout<<"ci abbiamo messo "<<mydurat<<" secondi"<<endl;
   myrmin = BenSol.getr_min();
   cout<<"inizia CPLEX"<<endl;

   Heuristic = BenSol.getHeurVal();
   ObjVal = BenSol.getObjVal();
   cout<<"Objval prima "<<ObjVal<<endl;

   contSol.LoadProblem(nnodes, nlinks, flows[i], links, nodes, MTU, 'C', 0);
   intSol.LoadProblem(nnodes, nlinks, flows[i], links, nodes, MTU, 'B', 0);
   rmin = intSol.getRmin();
  // imprContSol.LoadProblem(nnodes, nlinks, flows[i], links, nodes, MTU, 'B', 1);

   cout<<"caricato il problema su CPLEX"<<endl;

   /*imprContSol.DCRstartTime();
   imprcontstat = imprContSol.Solve();
   imprContSol.DCRstopTime();

   imprcontBound = imprContSol.getObjVal();*/

   //contSol.DCRstartTime();
 /*  contstat = contSol.Solve();
   //contSol.DCRstopTime();

   contBound = contSol.getObjVal();

 /* intSol.DCRstartTime();
   intstat = intSol.Solve();
   intSol.DCRstopTime();

   intBound = intSol.getObjVal();
   rmin = intSol.getRmin();

   contdurat = contSol.getTime();
   intdurat = intSol.getTime();
   
   imprContdurat = imprContSol.getTime();*/
  /* cout<<"differenza relativa tra gli rmin "<<rmin<<" contro il mio: "<<myrmin<<endl;
   cout<<"objval = "<<ObjVal<<" contro il valore di cplex: "<<intBound<<endl;   
   cout<<"l'euristica invece dice "<<Heuristic<<endl;
   cout<<"rilassamento continuo: "<<contBound;}//<< " mentre quello migliorato: "<<imprcontBound<<endl;
/*
   HB = (Heuristic - ObjVal)/Heuristic;
   HI = (Heuristic - intBound)/intBound;
   IB = (intBound - ObjVal)/ObjVal;
   //BCl = (intBound - contBound)/contBound;
   //cout<<"Objval dopo e cont2 "<<ObjVal<<" "<<imprcontBound<<endl;
   //BCu = (intBound - imprcontBound)/imprcontBound;
   cout<<" ";
   fout<<HB<<"\t"<<HI<<"\t"<<IB<<"\t"<</*BCl<<"\t"<<BCu<<"\t"<<*//*mydurat<<"\t"<<intdurat*//*<<"\t"<<imprContdurat<<"\t"<<contdurat<<"\t"*/;
/*   cout<<" "<<endl;
   if(IB < -5e-6) fout<<"!!!!!!";
   else{if(HI > 1e-5) fout<<"subopt";}
   if(IB < -1e-4) {cout<<"errore al flusso "<<i<<endl; exit(1);}
   fout<<endl;
   cout<<"servito il flusso "<<i<<endl<<endl;
  }

fout.close();*/

//////////////////////////////////////////////////////////////////////////////

 /*//////Plotting function, usare i metodi in BenBound e poi sfruttare matlab.
  ofstream fioutx("xFile3.txt");
  ofstream fiouty("yFile3.txt");
  int flowtoplot = 11;
  vector<double> xPlot;
  vector<double> yPlot;
  int npoints = 1000;

  BenSol.LoadProblem(nnodes, nlinks, flows[flowtoplot], links, nodes, MTU);
  cout<<"ho caricato tutto su Benders, inizio a cercare i punti"<<endl;
  xPlot = BenSol.getxPlotData();
  cout<<"inizio a trovare le altezze dei punti"<<endl;
  yPlot = BenSol.getyPlotData();

  for( i = 0; i < npoints; i++)
    fioutx<<xPlot[i]<<"\t";

  cout<<endl<<"---------------------------------------------"<<endl;

  for( i = 0; i < npoints; i++)
    fiouty<<yPlot[i]<<"\t";  

  cout<<endl<<"----------------------------------------"<<endl;
  fioutx.close();
  fiouty.close();
  

 /*
  BenSol.LoadProblem(nnodes, nlinks, flows[flowtoplot], links, nodes, MTU);
  BenSol.Solve();*/

  //////////////////////////////////////////////////////////////////////

  //CPL_BEND CplBen;

  //CplBen.LoadProblem(nnodes, nlinks, flows[0], links, nodes, MTU);
  //BenStat = CplBen.Solve();
  //ObjVal = CplBen.getObjVal();
  //CplBen.DCRWrite("myproblem.lp");
  //cout<<"solver status = "<<BenStat<<endl;
  //cout<<"ObjVal = "<<ObjVal<<endl;

/* for(i = 0; i< ntests; i++)
  {
   for(j = 0; j < nlinks; j++)///random costs
    flows[i].costs[j] = distribution(generator)*2.3;
 
  }*/
  /*BenSol.LoadProblem(nnodes, nlinks, flows[72], links, nodes, MTU);

  BenSol.Solve();

  ObjVal = BenSol.getObjVal();

  cout<<"ObjVal = "<<ObjVal<<endl;*/
  //cout<<"flow.burst = "<<flows[1].burst<<endl;
  //cout<<"flow.deadline = "<<flows[1].deadline<<endl;

  double* xsol, *rsol;
  //double rmin;
  double barl;
  vector<int> arcstoclose;
  vector<double> checksol;
  vector<int> checkposp;
  vector<int> checkposn;

  //arcstoclose = new int[1];

  xsol = new double[nlinks];
  rsol = new double[nlinks];

  ///////////////////////////////////////////////////////////////////////*/
///test correttezza Cplex con piccolo lagrangiano (comb convessa sol)
/* for(i = 11; i < 12; i++)
 {
    intSol.DCRsetTime(is);
    BenSol.DCRsetTime(is);
       BenSol.LoadProblem(nnodes, nlinks, flows[i], links, nodes, MTU);
    intSol.LoadProblem(nnodes, nlinks, flows[i], links, nodes, MTU, 'C',0);
    intstat = intSol.Solve();
    
     cout<<"risolto con cplex"<<endl;
 
    cout<<"caricato il problema su ben"<<endl;
    BenSol.Solve();

    cout<<"risolto con Benders"<<endl;
   // ObjVal = BenSol.getObjVal();
    //cout<<"risolv il benders"<<endl;
    //ApproxSol = BenSol.getApproxVal();

    cpxObjVal = intSol.getObjVal();
    rmin = intSol.getRmin();

    xsol = intSol.getXSol();
    rsol = intSol.getRSol();
    myrmin = BenSol.getr_min();
    //rmin = intSol.getRmin();
    cout<<"myrmin: "<<myrmin<<endl;

    lagSol.LoadProblem(nnodes, nlinks, flows[i], links, nodes, MTU, myrmin);

    lagSol.Solve();
    
    cout<<"risolto il lagrangiano"<<endl;
    checksol = lagSol.getCheckSol();
    checkposp = lagSol.getCheckSolPos();
    checkposn = lagSol.getCheckSolNeg();

    cout<<"soluzione nuova in r_ij = ";
    for(j=0; j<checksol.size(); j++)
      cout<< checksol[j] <<"-";

    cout<<endl;

    cout<<"posizioni positive : ";
    for(j=0; j<checkposp.size(); j++)
      cout<< checkposp[j] <<"-";

    cout<<endl;

    cout<<"posizioni negative : ";
    for(j=0; j<checkposn.size(); j++)
      cout<< checkposn[j] <<"-";

    cout<<endl;
    cout<<"soluzione CPLEX in r_ij = ";
    for(j=0; j<checksol.size(); j++)
      cout<< rsol[j] <<"-";

    cout<<endl;

    cout<<"soluzione CPLEX in x_ij = ";
    for(j=0; j<checksol.size(); j++)
      cout<< xsol[j] <<"-";

    cout<<endl;

 }*/
    /*lagSol.LoadProblem(nnodes, nlinks, flows[i], links, nodes, MTU, rmin);

    counter = 0;
    lagObjVal = 0;
    delay = 0;
    
    arcstoclose.clear();

    for(j = 0; j < nlinks; j++)
    {
      if(xsol[j] < 1e-9)
      { arcstoclose.push_back(j);
        counter++;}
    }
    //cout<<"sono qui 1"<<endl;
    //cout<<"gli archi da chiudere sono "<<counter<<" su un totale di "<<nlinks<<" archi"<<endl;

    lagSol.closeArcs(arcstoclose, counter);

    //cout<<"sono qui 2"<<endl;

    lagSol.Solve();

    lagstat = lagSol.getStatus();

    cout<<"lo status del lagrangiano è "<<lagstat;

    //cout<<"sono qui 3"<<endl;

    checksol.clear();
    checkpos.clear();

    checksol = lagSol.getCheckSol();
    checkpos = lagSol.getCheckSolPos();

    for(j = 0; j < checksol.size(); j++)
    {
      if(checksol[j] < rmin)
        {cout<<"viola il vincolo su rmin per "<<-(checksol[j]-rmin)/rmin<<endl;}

      if(checksol[j] > links[checkpos[j]].capacity)
        {cout<<"viola il vincolo su c_ij per "<<(checksol[j]-links[checkpos[j]].capacity)/links[checkpos[j]].capacity<<endl;}
      delay += MTU / checksol[j] + MTU / links[checkpos[j]].speed + links[checkpos[j]].delay + nodes[links[checkpos[j]].startnode].delay;
      lagObjVal += checksol[j]*flows[i].costs[checkpos[j]];
    }
   // cout<<"delay partial"<<delay<<endl;
     delay += flows[i].burst / rmin - flows[i].deadline;

    cout<<"il vincolo del ritardo dà un valore: "<<delay<<endl;
    cout<<"la differenza relativa tra il mini lagr e il lower bound è: "<<(lagObjVal- ApproxSol)/ApproxSol<<endl;
    cout<<"quella tra cplex e il mini lagr invece è: "<<(cpxObjVal - lagObjVal)/lagObjVal<<endl;

    cout<<"servito il flusso "<<i<<endl<<endl<<"------------------------"<<endl<<endl;

    if( i = 24) i =48;
          
 }

////////////////////////////////////////////////////////////////////////////////////
/////Test correttezza Cplex:
  /*for(i = 0; i < ntests + 1; i++)
  { 
    fo = 0;
    intSol.DCRsetTime(is);
    intSol.LoadProblem(nnodes, nlinks, flows[i], links, nodes, MTU, 'B');
    intstat = intSol.Solve();

    cpxObjVal = intSol.getObjVal();
    rmin = intSol.getRmin();

   // cout<<"rmin = "<<rmin<<endl;
   // cout<<"mtu = "<<MTU<<endl;

    xsol = intSol.getXSol();
    rsol = intSol.getRSol();
    delay = 0;
    for(j = 0; j < nlinks; j++)
    {
      if(xsol[j] > 1e-9)
      {
        if(1-xsol[j]!=0)
        cout<<"precisione sulle xsol: "<<1-xsol[j]<<endl;
        if((rsol[j]-rmin)<0)
        cout<<"precisione relativa su rmin: "<<(rsol[j]-rmin)/rmin<<endl;
        //if(rsol[j]-rmin < 0)
        //  {rsol[j] = rmin;}
        //cout<<"precisione relativa su c_ij: "<<(links[j].capacity - rsol[j])/rsol[j]<<endl;
        //cout<<"dove x_ij era: "<<xsol[j]<<endl;

        delay += MTU / rsol[j] + MTU / links[j].speed + links[j].delay + nodes[links[j].startnode].delay;
        fo += flows[i].costs[j] * rsol[j];
       // cout<<"partial delay "<<delay<<endl;
      }
    }
     //cout<<"partial delay  "<<delay<<endl;
    delay += flows[i].burst / rmin - flows[i].deadline;

    cout<<"precisione sul ritardo: "<<delay<<endl;
    //cout<<"precisione relativa sulla fo: "<<(fo- cpxObjVal)/cpxObjVal<<endl;
    //cout<<"infatti fo: "<<fo<<endl;
    //cout<<"cplex ObjVal: "<<cpxObjVal<<endl;
    cout<<"servito il flusso "<<i<<endl<<endl;
    
   cout<<"----------------------------------------------------"<<endl<<endl;
 

  }*/



////////////////////////////////////////////////////////////////////////////
cout<<"ntest = "<<ntests<<endl;
 /////////Test Vari
 for(i = 0; i < min(ntests,100); i++) 
  {

   cout<<"flow:"<<i<<endl;
   intSol.DCRsetTime(is);
   contSol.DCRsetTime(is);
   BenSol.DCRsetTime(is);
 // cout<<"valore base max tra "<<flows[i].rate<<" e "<<flows[i].burst/flows[i].deadline<<endl;

  //cout<<"il lagrangiano risolve"<<endl;
   BenSol.LoadProblem(nnodes, nlinks, flows[i], links, nodes, MTU);
   cout<<"problem loaded"<<endl;
  //lagSol.LoadProblem(nnodes, nlinks, flows[i], links, nodes, MTU,10000);

 //cout<<"status solver"<<xstat<<endl;

   //xSOL.LoadProblem(nnodes, nlinks, flows[i], links, nodes, MTU,40000);
 
 //xstat = xSOL.Solve();
 //cpxObjVal = xSOL.getObjVal();
 //cout<<"status solver"<<xstat<<endl;
 //cout<<"  cplxBound per 40000= "<<cpxObjVal<<endl;
   BenSol.DCRstartTime();
   BenSol.Solve();
   BenSol.DCRstopTime();

   Heuristic = BenSol.getHeurVal();
  //lagSol.Solve();
  

   ObjVal = BenSol.getObjVal();
   ApproxSol = BenSol.getApproxVal();

   myrmin = BenSol.getr_min();

/*
  //rfissatoSOL.LoadProblem(nnodes, nlinks, flows[i], links, nodes, MTU,rmin);
  contSol.LoadProblem(nnodes, nlinks, flows[i], links, nodes, MTU, 'C',0);
  intSol.LoadProblem(nnodes, nlinks, flows[i], links, nodes, MTU, 'B',0);

  //intSol.DCRWrite("myprob.lp");
  
  contSol.DCRstartTime();
  contstat = contSol.Solve();
  contSol.DCRstopTime();

  contBound = contSol.getObjVal();

  intSol.DCRstartTime();
  intstat = intSol.Solve();
  intSol.DCRstopTime();

  intBound = intSol.getObjVal();
  
  cplrmin = intSol.getRmin();

  lagSol.LoadProblem(nnodes, nlinks, flows[i], links, nodes, MTU, cplrmin);
  lagSol.Solve();
  lagstat = lagSol.getStatus();

  //lambda = lagSol.getLambda();
 // beta = lagSol.getOptBeta();

 // lagObjVal -= beta*lambda;

  cout<<"differenza relativa tra il lower bound di Benders e il valore di cplex: "<<-(ApproxSol - intBound)/intBound<<endl;
  //xstat = rfissatoSOL.Solve();
  //cpxObjVal = rfissatoSOL.getObjVal();

  //cout<<"status di CPLEX: "<<contstat<<endl; 
  
  //cout<<"il lagrangiano con rmin = "<<cplrmin<<" dà valore: "<<lagObjVal<<" relativo ad un lambda pari a "<<lambda<<endl;
  //cout<<"Bound di CPLEX continuo: "<<contBound<<endl;

  contdurat += contSol.getTime();
  intdurat += intSol.getTime();
  mydurat += BenSol.getTime();
 // 

 // cout<<"punti di minimo (Ben-cplex): "<<myrmin<<" - "<<cplrmin<<endl;

  //cout<<"differenza relativa tra i punti: "<<(myrmin - cplrmin)/cplrmin<<endl;
 // cout<<"differenza relativa tra i valori: "<<(ObjVal-intBound)/intBound<<endl;

 // cout<<"differenza relativa il valore della f.o. nel punto ed il miglior lower bound: "<<(lagObjVal - ApproxSol)/ApproxSol<<endl;
  //cout<<"con beta = "<<beta<<endl;
  
  lagObjVal = lagSol.getObjVal();
  lambda = lagSol.getLambda();

  if(abs(ObjVal-cpxObjVal)>100)
  cout<<"differenza relativa tra la f.o. nel punto e l'ottimo di cplex intero: "<<(lagObjVal-intBound)/lagObjVal<<endl;
  cout<<"valore approssimato per ObjVal = "<<ApproxSol<<endl;
  cout<<"ObjVal = "<<ObjVal<<" numero di iterazioni: "<<BenSol.getNumIterationBender()<<" valore di rmin: "<<BenSol.getr_min()<<endl;
  cout<<" lagBound = "<<lagObjVal<<", lambda = "<<lambda<<endl;
  cout<<"Bound di CPLEX intero: "<<intBound<<endl;
  cout<<"Bound di CPLEX continuo: "<<contBound<<endl;
  cout<<"Bound dell'euristica: "<<Heuristic<<endl;
  cout<<"servito il flusso "<<i<<endl;
  cout<<endl<<endl;

  if(lambda == 0 && lagstat == 0 && BenSol.getNumIterationBender()>0) {cout<<"Lagstat = "<<lagSol.getStatus()<<endl;exit(0);}*/

  cout<<"------------------------"<<endl<<endl;

}

 //cout<<"media tempi per CPLEX continuo: "<<contdurat/i<<endl;
 //cout<<"media tempi per BENdDCR: "<<mydurat/i<<endl;
 //cout<<"media tempi per CPLEX intero: "<<intdurat/i<<endl;

  return 0;

}


void ReadData(char *filename, char format, int * nnodes, int *nlinks, int *nflows, double *MTU,
DCR::DCRFlow **flows, DCR::DCRLink **links, DCR::DCRNode **nodes, double perccap, double percdelay, int seed)
{

	DCRGraph dcrgraph(filename,format,perccap,percdelay,seed);

	(*nnodes) = dcrgraph.getNumNodes();

	(*nlinks) = dcrgraph.getNumLinks();

	(*nflows) = dcrgraph.getNumFlows();

	(*MTU) =  dcrgraph.getMTU();

	(*nodes) = dcrgraph.getNodes();

	(*links) = dcrgraph.getLinks();

	(*flows) = dcrgraph.getFlows();

//cout<<"grafo aggiornato"<<endl;
}
