/*--------------------------------------------------------------------------*/
/*-------------------------- File SPT.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------------- IMPLEMENTATION -------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "SPT.h"
#include "OPTUtils.h"

#include <iostream>
#include <cstdlib>

using namespace std;

/*--------------------------------------------------------------------------*/
/*--------------------- IMPLEMENTATION OF SPT-----------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------------- PUBLIC METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*----------------------------- CONSTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

  SPT::SPT()
  {
    Links.resize(1);
    nhops = 0;
    s = 0;
    t = 0;
    stat = OK;
    Sol.resize(1); //serve?
    DualSol.resize(1);

  }         


/*--------------------------------------------------------------------------*/
/*----------------------------------LOAD DATA ------------------------------*/
/*--------------------------------------------------------------------------*/

  void SPT::LoadProblem (int nnodes, int nlinks, vector<SPTLink> links, int sourcenode, int sinknode)
  {

   clean_up();
   
   stat = OK;
   s = sourcenode;
   t = sinknode;
   numNodes = nnodes;
   numLinks = nlinks;
   copyDataArrays(links);

  } 

/*--------------------------------------------------------------------------*/

  void SPT::updCosts(vector<SPTLink> links)
  {
    int i;
    
    Sol.clear();
    DualSol.clear();

   //delete[] Sol;
   //delete[] DualSol;
   //Links = new SPT::SPTLink[numLinks];

   stat = OK;
   
    for(i = 0; i < numLinks; i++)
      Links[i].cost = links[i].cost;
  }

/*--------------------------------------------------------------------------*/
/*--------------------------- GET METHODS-----------------------------------*/
/*--------------------------------------------------------------------------*/

   vector<int> SPT::getPath(void)
   {
    return(Sol);
   }

/*--------------------------------------------------------------------------*/

   SPT::Status SPT::getStatus()
   {
    return(stat);
   }

/*--------------------------------------------------------------------------*/

   vector<double> SPT::getLabel()
   {
    return(DualSol);
   }

/*--------------------------------------------------------------------------*/

   int SPT::getNHops()
   {
    return(nhops);
   }
   
/*--------------------------------------------------------------------------*/
/*----------------------------------SOLVE-----------------------------------*/
/*--------------------------------------------------------------------------*/

   void SPT::Solve()
   {
    int i, j, h;
    //int maxIter = pow(numNodes,2);
    //int iterat=0; //se supera n^2, allora vi sarà almeno un ciclo di costo negativo (?)
    double *distance = new double[numNodes];
    int *previous = new int[numNodes];
    int *Q = new int[numNodes]; //FIFO queue
    double dist, inf = Inf<double>();
    int next;
    int HEAD; //FIFO head
    int TAIL; //FIFO tail

    //inizializzazione delle strutture
    for(i = 0; i < numNodes; i++)
    {
      distance[i] = inf; 
      previous[i] = -1;     
      Q[i] = -1;  
    }
    distance[s] = 0;
    HEAD = s;
    TAIL = s;

    //solve SHORTEST PATH
    //////////////////////////////////////
    while (HEAD != -1)
    {
      //extract node from queue Q
      h = HEAD;
      HEAD = Q[HEAD];
      Q[h] = -1;
      if(HEAD == -1) TAIL = -1;
  
      //consider all neighbours of h
      for(j = 0; j < numLinks; j++)
      {
        if(Links[j].startnode == h) 
        {
          next = Links[j].endnode;
          dist = Links[j].cost + distance[h];
      
          if(dist < distance[next])
          {
            distance[next] = dist;
            previous[next] = h;
            
            //insert node in queue Q
            if(HEAD == -1) { HEAD = next; TAIL = next; }
            else if((HEAD != next) && (TAIL != next) && (Q[next] == -1)) { Q[TAIL] = next; TAIL = next; }
                        
          }//if (distance label update)
        
        }//if (neighbour arc)
      
      }//for (all arcs)
    }//while (head != -1)
    //Sistema le variabili per i get methods.

    ////////va bene così? oppure new[nnodes] e ignorare quanto fatto in load pr, oppure delete e new??
    DualSol.resize(numNodes);
    //cout<<"distanza totale in SP = "<<distance[t]<<endl;

    for(j = 0; j < numNodes; j++) 
        DualSol[j] = distance[j]; //serve il for?

    i=t;
    int* indices = new int[numNodes];
    int index;
    nhops = 0; //serve per le chiamate successive

    while(i != s && stat == OK) 
     {
      index = getLink(previous[i],i); //se non trova un link conclude che il grafo è sconnesso
      indices[nhops] = index; //copia a ritroso gli indici del cammino ottimo
      i = previous[i];
      nhops++;  //salva il numero di questi indici
     }

    if(stat == OK)
     {
     Sol.resize(nhops);

     for(j = 0; j < nhops; j++)
      {
       Sol[nhops-j-1]=indices[j];  //prende il cammino ottimo rimettendolo in ordine s-t.
      }
     }
    
   delete [] indices;
   delete [] previous;
   delete [] distance;
   delete [] Q;
   }



/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

  SPT::~SPT()
   {
	  clean_up();
   }

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

   void SPT::copyDataArrays(vector<SPTLink> links)
   {
   	int i, j;

		Links.resize(numLinks);
	
  	for(j = 0; j < numLinks; j++)
	  {
		  Links[j].startnode = links[j].startnode;
		  Links[j].endnode = links[j].endnode;
		  Links[j].cost = links[j].cost;
      //inutile copiare Link.rstar a questo livello.
	  }	

   }

/*--------------------------------------------------------------------------*/

   int SPT::getLink(int from, int to)
   {
    int i;
  
     for(i = 0; i < numLinks; i++)
     {
       if(Links[i].startnode == from && Links[i].endnode == to)
         return i; 
     } 
    //cout<<"rmin troppo grande: si è sconnesso il grafo\n";
    stat = Error;
    return -1;
   }


/*--------------------------------------------------------------------------*/
 
   void SPT::clean_up()
   { 
    //cout<<"Links"<<endl;
	  Links.clear();
    //cout<<"Sol"<<endl;
	  Sol.clear();
    //cout<<"DualSol"<<endl;
    DualSol.clear();
   }  


/*--------------------------------------------------------------------------*/
/*------------------------------ End SPT.cpp -------------------------------*/
/*--------------------------------------------------------------------------*/