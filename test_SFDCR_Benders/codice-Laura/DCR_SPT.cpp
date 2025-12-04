/*--------------------------------------------------------------------------*/
/*---------------------------- File DCR_SPT.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------------- IMPLEMENTATION -------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "DCR_SPT.h"
#include <fstream>
#include <cstdlib>
#include <set>

/*--------------------------------------------------------------------------*/
/*----------------------- IMPLEMENTATION OF DCR_SPT ------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------------- PUBLIC METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*----------------------------- CONSTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

//#define DEBUG

DCR_SPT::DCR_SPT() : DCR()
{
	
	Nodes = 0;
	Links = 0;
	Flows = 0;
	
	Xsol = 0;
	Rsol = 0;	
	
	nhops = 0;
	
	//FIXME: is this okay ???
	heur = 1;
	
	//FIXME: is this okay ???
	objval = Inf<double>();
}

/*--------------------------------------------------------------------------*/
/*----------------------------------SOLVE-----------------------------------*/
/*--------------------------------------------------------------------------*/

DCR::DCRStatus DCR_SPT::DCRsolve( void )
{

	//select a heuristic method	
	if(heur == 1)	
	{
		DCRheurERAI();
		//debug
		//cout << "\nSPT solver option ERA-I" << endl;
	}
	if(heur == 2)
	{
		DCRheurERAH();
		//debug
		//cout << "\nSPT solver option ERA-H" << endl;
	}
	
	if(UBsol) return OK;
	
	return Infeasible;
}


/*--------------------------------------------------------------------------*/
/*----------------------------------GET RESULTS-----------------------------*/
/*--------------------------------------------------------------------------*/

double DCR_SPT::DCRgetObj( void )
{
	return objval;
}


/*--------------------------------------------------------------------------*/

int DCR_SPT::DCRgetUBSolNPaths( int k )
{
	if(k < 0 || k>= numFlows)
		throw( DCR::DCRException( "DCR_SPT::DCRgetUBSolNPaths(): wrong flow index." ) );
	if(UBsol) return 1;
	else
		throw( DCR::DCRException( "DCR_SPT::DCRgetUBSolNPaths(): no solution found." ) );	//FIXME: is this okay???
}

/*--------------------------------------------------------------------------*/

int DCR_SPT::DCRgetUBSolPath( int k, int p, int *X, double *R)
{
	int i;
	
	if(k < 0 || k>= numFlows)
		throw( DCR::DCRException( "DCR_SPT::DCRgetUBSolPath(): wrong flow index." ) );
	if(p != 0)
		throw( DCR::DCRException( "DCR_SPT::DCRgetPUBolPath(): wrong path index." ) );
	if(UBsol)
	{
	
		for(i = 0; i < nhops; i++)
		{
			X[i] = Xsol[i];
			R[i] = Rsol[i];
		}		
	}
	else
		throw( DCR::DCRException( "DCR_SPT::DCRgetUBSolPath(): no solution found." ) );//FIXME: is this okay ???		
	
	return nhops;
}

/*--------------------------------------------------------------------------*/

int DCR_SPT::DCRgetPSolNPaths( int k )
{
	throw( DCR::DCRException( "DCR_SPT::DCRgetPSolNPaths(): Continuous relaxation is not supported...sorry." ) );	
}

/*--------------------------------------------------------------------------*/

void DCR_SPT::DCRgetPSolPath(int k, int p, double *X, double *R)
{
	throw( DCR::DCRException( "DCR_SPT::DCRgetPSolPath(): Continuous relaxation is not supported...sorry." ) );
}


/*--------------------------------------------------------------------------*/
/*----------------------------------LOAD DATA ------------------------------*/
/*--------------------------------------------------------------------------*/

void DCR_SPT::DCRloadProblem(int nnodes, int nlinks, int nflows, 
	DCRFlow *flows, DCRLink *links, DCRNode *nodes, double mtu, DCRDelay deltype)
{
	int status;
	
	if(deltype != SRP)
		throw( DCR::DCRException( "DCR_SPT::DCRloadProblem(): this delay formula is not supported yet..sorry" ) );
	
	if(nflows > 1)
		throw( DCR::DCRException( "DCR_SPT::DCRloadProblem(): this class can only solve SFSP DCR problems..sorry" ) );
	
	///////////////////////////////
	//INIT
	/////////////////////////////
	clean_up();	
	
	UBsol = false;
	objval = Inf<double>();
	nhops = 0;
	////////////////////////////////
	
	numNodes = nnodes;
	numLinks = nlinks;
	numFlows = nflows;	
	MTU = mtu;
	dtype = deltype;	
	copyDataArrays(flows, links, nodes);
	Xsol = new int[nnodes-1];
	Rsol = new double[nnodes-1];
	
}

/*--------------------------------------------------------------------------*/
/*----------------------------------SETTERS---------------------------------*/
/*--------------------------------------------------------------------------*/


void DCR_SPT::DCRsetHeur(char h)
{
	if(h == '1') heur = 1;
	else if(h == '2') heur = 2;
	else throw( DCR::DCRException( "DCR_SPT::DCRsetHeur(): this heur. type is not supported yet..sorry" ) );
}

/*--------------------------------------------------------------------------*/
/*----------------------------------MODIFIERS-------------------------------*/
/*--------------------------------------------------------------------------*/

void DCR_SPT::DCRcloseArcs( int * whch, int na )
{
	int k, j, ai;
	
	for(k = 0; k < numFlows; k++)
	{
			for(j = 0; j < na; j++)
			{
				ai = whch[j];
				Flows[k].costs[ai] = Inf<double>();
			}		
	} 
}

/*--------------------------------------------------------------------------*/

void DCR_SPT::DCRcloseArcs( int k , int * whch, int na )
{
	int j, ai;
	
	for(j = 0; j < na; j++)
	{
		ai = whch[j];
		Flows[k].costs[ai] = Inf<double>();
	}		
}

/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

DCR_SPT::~DCR_SPT()
{

	clean_up();
}

/*--------------------------------------------------------------------------*/
/*------------------------------ LOCAL (PRIVATE) METHODS -------------------*/
/*--------------------------------------------------------------------------*/

/* This private method implements ERA-I heuristic algorithm
 *  for SFSP DCR problems, see Orda paper.*/
void DCR_SPT::DCRheurERAI()
{
	//Dijkstra
	int s, t, i, j, h, next;
	double tcost;
	int lindex, nh;
	double inf, dist, rmin, cost;
	double *distance = new double[numNodes];
	int *previous = new int[numNodes];
	int *Q = new int[numNodes]; //FIFO queue
	int HEAD; //FIFO head
	int TAIL; //FIFO tail
	std::set<double> myset;
	std::set<double>::iterator it;
		
	inf = Inf<double>();//used to find the minimum
	cost = -1;//init cost
	
	s = Flows[0].sourcenode;
	t = Flows[0].sinknode;

	
	for(i = 0; i < numLinks; i++)
		myset.insert(Flows[0].caps[i]);//NB always use individual capacity
	
	//consider REDUCED GRAPH	
	////////////////////////////////////
    for (it=myset.begin(); it!=myset.end(); ++it)
    {
		rmin = *it;
				
		if(rmin < Flows[0].rate) continue;
	
		//init SHORTEST PATH DATA STRUCTURES
		//////////////////////////////////////
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
				if(Flows[0].caps[j] < rmin) continue;//NB always use individual capacity
								
				if(Links[j].startnode == h) 
				{
					next = Links[j].endnode;
					dist = MTU/Flows[0].caps[j] + MTU/Links[j].speed + Nodes[h].delay + Links[j].delay + distance[h];
			
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
		
		//CHECK SOLUTION 
		//////////////////////////////////////
		//check feasibility
		if((distance[t] + Flows[0].burst/rmin) <= Flows[0].deadline)
		{
			//calculate cost
			tcost = 0;
			i = t;
			while( i != s )
			{
				lindex = DCRgetLink(previous[i], i);
				tcost += Flows[0].costs[lindex]*Flows[0].caps[lindex];
				i = previous[i];
			}	
			
			//maybe update best solution
			if((tcost < cost) || (cost < 0)) 
			{ 
				cost = tcost; 
		
				i = t;
				j = 0;	
				while( i != s )
				{
					lindex = DCRgetLink(previous[i], i);
					Xsol[j] = lindex;
					Rsol[j] = Flows[0].caps[lindex];
					j++;
					i = previous[i];
				}//while
					
				nhops = j;
				 				
			}//if (sol update)		
		}//if (feasible sol)
		
	}//for (all capacities values -> reduced graph)

	delete [] distance;
	delete [] previous;
	delete [] Q;			

	if( cost  < 0 ) { UBsol = false; return; }

	UBsol = true;
	objval = cost;
	
	/*cout << "\n ERA-I cost is " << cost << " with #hops " << nhops << endl;*/
		
}

/*--------------------------------------------------------------------------*/

/* This private method implements ERA-H  heuristic algorithm
 * for SFSP DCR problems, see Orda paper.*/
void DCR_SPT::DCRheurERAH()
{
	int s, t, i, j, h, next, lindex, nh;
	double tcost, tr0, dl, inf, dist, rmin, cost, r0;
	double *distance = new double[numNodes];
	int *previous = new int[numNodes];
	int *Q = new int[numNodes]; //FIFO queue
	int HEAD; //FIFO head
	int TAIL; //FIFO tail
	std::set<double> myset;
	std::set<double>::iterator it;
	bool found;
	
	inf = Inf<double>();
	cost = Inf<double>();
	found = false;

	s = Flows[0].sourcenode;
	t = Flows[0].sinknode;
	
	for(i = 0; i < numLinks; i++)
		myset.insert(Flows[0].caps[i]);
	
	//consider REDUCED GRAPH	
	////////////////////////////////////
    for (it=myset.begin(); it!=myset.end(); ++it)
    {
		rmin = *it;
		
		if(rmin < Flows[0].rate) continue;
	
		//init SHORTEST PATH DATA STRUCTURES
		//////////////////////////////////////
		for(i = 0; i < numNodes; i++)
		{
			distance[i] = inf; 
			previous[i] = s;			
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
			
			//CHECK SOLUTION 
			//////////////////////////////////////		
			if(h == t)
			{
						
				if((distance[t] + Flows[0].burst/rmin) < Flows[0].deadline)
				{
											
					//calculate r0					
					dl = 0;
					nh = 0;
					i = t;
					while( i != s )
					{
						lindex = DCRgetLink(previous[i], i);
						dl += MTU/Links[lindex].speed + Nodes[previous[i]].delay + Links[lindex].delay;
						nh ++;
						i = previous[i];
					}				
					tr0 = (Flows[0].burst + nh*MTU)/(Flows[0].deadline - dl);
					if (tr0 < Flows[0].rate)	tr0 = Flows[0].rate;
							
					//calculate cost for r0		
					tcost = 0;										
					i = t;			
					while( i != s )
					{
						lindex = DCRgetLink(previous[i], i);
						tcost += Flows[0].costs[lindex] * tr0;
						i = previous[i];
					}	
										
					//maybe update solution 
					if(tcost < cost) 
					{ 
						cost = tcost; 
						r0 = tr0; 
						found = true; 

						j = 0;
						i = t;
						while( i != s )
						{
							lindex = DCRgetLink(previous[i], i);
							Xsol[j] = lindex;
							Rsol[j] = r0;
							j++;
							i = previous[i];
						}//while
						
						nhops = j;
					}//if (update sol)		
				
				}//if (feasible sol)			
			
			}//if (sink is extracted -> check solution)
	
			//consider all neighbours of h 
			for(j = 0; j < numLinks; j++)
			{
				if(Flows[0].caps[j] < rmin) continue;
		
				if(Links[j].startnode == h) 
				{
					next = Links[j].endnode;
					dist = MTU/rmin + MTU/Links[j].speed + Nodes[h].delay + Links[j].delay + distance[h];
			
					if(dist < distance[next])
					{
						distance[next] = dist;
						previous[next] = h;
						
						//insert next node in queue Q if not already present
						if(HEAD == -1) { HEAD = next; TAIL = next; }
						else if((HEAD != next) && (TAIL != next) && (Q[next] == -1)) { Q[TAIL] = next; TAIL = next; }
												
					}//if (distance label update)
				
				}//if (neighbour arc)
			
			}//for (all arcs)
		
		}//while (head != -1)
		
	}//for (all capacities values -> reduced graph)
	
	delete [] distance;
	delete [] previous;
	delete [] Q;	
		

	if( !found ) { UBsol = false; return; }

	UBsol = true;
	objval = cost;
	
	/*cout << "\n ERA-H cost is " << cost << " with #hops " << nhops << endl;*/

}

/*--------------------------------------------------------------------------*/

/* This private method returns the index of the link 
 * having the given startnode (from) and endnode (to). */
int DCR_SPT::DCRgetLink(int from, int to)
{
	int i;
	
	for(i = 0; i < numLinks; i++)
	{
		if(Links[i].startnode == from && Links[i].endnode == to)
			return i;	
	}	
}

/*--------------------------------------------------------------------------*/

/* This private method deallocates dinamically allocated memory. */
void DCR_SPT::clean_up()
{
		
	if (Nodes) {delete [] Nodes; Nodes = 0;}
	if (Links) {delete [] Links; Links = 0;}
	if (Flows) {
		for(int i = 0; i < numFlows; i++) {
			delete [] Flows[i].caps;
			delete [] Flows[i].costs;
		}
		delete [] Flows; 
		Flows = 0;
	}
	
	if (Xsol) {delete [] Xsol; Xsol = 0;}
	if (Rsol) {delete [] Rsol; Rsol = 0;}
	
}

/*--------------------------------------------------------------------------*/

/* This private method is used to copy array data (flows, links and nodes) from the load() function 
 * into the class as class members.*/
void DCR_SPT::copyDataArrays(DCRFlow *flows, DCRLink *links, DCRNode *nodes)
{
	int i, j, k;
	
	try
	{
		Nodes = new DCR::DCRNode[numNodes];
	} 
	catch(exception& e)
	{ cout << "Standard exception: " << e.what() << endl; }
	
	for(i = 0; i < numNodes; i++)
	{
		Nodes[i].delay = nodes[i].delay; 
	}
	
	
	try
	{
		Links = new DCR::DCRLink[numLinks];
	}
	catch(exception& e)
	{cout << "Standard exception: " << e.what() << endl;}
	
	for(j = 0; j < numLinks; j++)
	{
		Links[j].speed = links[j].speed; 
		Links[j].capacity = links[j].capacity; 
		Links[j].delay = links[j].delay; 
		Links[j].startnode = links[j].startnode;
		Links[j].endnode = links[j].endnode;
	}	
	
	
	try
	{
		Flows = new DCR::DCRFlow[numFlows];
	}
	catch(exception& e)
	{ cout << "Standard exception: " << e.what() << endl; }
	
	for(k = 0; k < numFlows; k++)
	{
		Flows[k].sourcenode = flows[k].sourcenode; 
		Flows[k].sinknode = flows[k].sinknode; 
		Flows[k].burst =  flows[k].burst;  
		Flows[k].rate = flows[k].rate; 
		Flows[k].deadline = flows[k].deadline; 
		 
		try
		{
		Flows[k].costs = new double[numLinks];
		}
		catch(exception& e)
		{ cout << "Standard exception: " << e.what() << endl; }
		
		for(j = 0; j < numLinks; j++)
			Flows[k].costs[j] = flows[k].costs[j];
		
		try
		{
		Flows[k].caps = new double[numLinks];
		}
		catch(exception& e)
		{ cout << "Standard exception: " << e.what() << endl; }
		
		for(j = 0; j < numLinks; j++)
			Flows[k].caps[j] = flows[k].caps[j];
	}
}


/*--------------------------------------------------------------------------*/
/*-------------------------- End File DCR_SPT.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
