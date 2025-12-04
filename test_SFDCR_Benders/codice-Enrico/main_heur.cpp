#include "DCR.h"
#include "DCRGraph.h"
#include "DCRLagrangianSolver.h"
#include "SPT.h"
#include "CPX_LAG.h"
#include <math.h>
#include <iostream>
#include <random>
#include <fstream>
#include <vector>

using namespace std;

void ReadData(char *filename, char format, int * nnodes, int *nlinks, int *nflows, double *MTU,
DCR::DCRFlow **flows, DCR::DCRLink **links, DCR::DCRNode **nodes, double perccap, double percdelay, int seed);


int main (int argc, char ** argv)
{ 
  char *filein;
  string filename = argv[1];
  char format = 'm';  //formato Mnetgen
  filein=argv[1];
  int nnodes;
  int nlinks;
  int nflows;
  double MTU;
  double ObjVal;


  DCR::DCRFlow *flows = 0;
  DCR::DCRLink *links = 0;
  DCR::DCRNode *nodes = 0;

  CPX_LAG::CPLStatus Cstat;
  CPX_LAG::CPLStatus Bstat;
  DCRLagrangianSolver::LAGStatus lagstat;
  
  ReadData(filein, format, &(nnodes), &(nlinks), &(nflows), &(MTU), &(flows), &(links), &(nodes), 0, 0, 1);

  CPX_LAG cpxRelax;
  CPX_LAG cpxEsatto;
  DCRLagrangianSolver lagSol;

  //cout<<"nr flow ="<<nflows<<endl;

  //ofstream fout(filename+"_heur"+".csv");
  
  double* bounds = new double[2]; //TODO: renderlo indipendente dai flussi per fare le medie
  double step = 0;
  int ntest = min(nflows,1000);
  double rmin;

  //i seguenti sono vettori perché catturano il comportamento in funzione di alpha
  double* aritmetic_mean_feas_DLgap = new double[10]; //vediamo come si comporta l'euristica

  double* nr_feas = new double[10];


  for(int j = 0; j<10; j++){ //tutto a zero

    aritmetic_mean_feas_DLgap[j] = 0;
    nr_feas[j] = 0;
  }


  for(int j = 0; j<nlinks; j++){ //cerchiamo c_max

      if(links[j].capacity > bounds[1]) bounds[1] = links[j].capacity;
  } 

  lagSol.LoadProblem(nnodes, nlinks, flows[394], links, nodes, MTU, 1);
  cpxEsatto.LoadProblem(nnodes, nlinks, flows[394], links, nodes, MTU, 1, 'B');


  for(int count_feas = 8236.15; count_feas < 11000; count_feas += 100){
          lagSol.updrmin(count_feas); //aggiorniamo rmin
        cpxEsatto.updateRmin(count_feas);

        lagSol.Solve();
        lagstat = lagSol.getStatus();

        Bstat = cpxEsatto.Solve();

        cout<< lagstat<< "   "<<Bstat<<endl;}

  //cout<<filename<<" #nodes = "<<nnodes<<" #links = "<<nlinks<<" #nflows = "<<nflows<<endl;
  
  /*for(int att_flow = 0; att_flow < ntest; att_flow++){
    
    bounds[0] = flows[att_flow].rate; // \rho

    lagSol.LoadProblem(flows[att_flow]);

    for(int j = 0; j<10; j++){

        lagSol.DCRsetTime(true); //resettiamo i timer

        rmin = bounds[0]+((double) j/10)*(bounds[1]-bounds[0]);

        lagSol.updrmin(rmin); //aggiorniamo rmin

        lagSol.DCRstartTime();
        lagSol.Solve();
        lagSol.DCRstopTime();
        lagstat = lagSol.getStatus();

        if(lagstat == 0){ //caso feasible => DL e cont feasible!
            
            nr_feas[j]++;
            aritmetic_mean_feas_DLgap[j] += abs(lagSol.getHeurVal() - lagSol.getObjVal())/lagSol.getHeurVal();
        }
    }//per ogni rmin nel range
  }//per ogni flusso

  //calcoliamo le medie dei dati ottenuti

  double* perc_unfeasDL = new double[10];
  double* perc_unfeasCPX_R = new double[10];
  double* perc_feas = new double[10];

  fout<<"alpha"<<"\t"<<"(heur-DL)/DL"<<endl;

  for(int j=0; j<10; j++){

      if(nr_feas[j]){

        aritmetic_mean_feas_DLgap[j] /= nr_feas[j];
      }

      //infine scriviamo tutto in un csv
      fout<<(double) j/10<<"\t"<<aritmetic_mean_feas_DLgap[j]<<endl;

  }

  fout.close();*/

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