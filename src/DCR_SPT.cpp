/*--------------------------------------------------------------------------*/
/*---------------------------- File DCR_SPT.cpp -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the DCR_SPT class.
 *
* \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Laura Galli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Luca Mencarelli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Enrico Sorbera \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 * 
 * \copyright &copy; by Antonio Frangioni
 */ 

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

// gives default values to all data members: no instance is loaded yet
// (Nodes/Links/Flows are null), the default heuristic is ERA-I (heur == 1)
// and no solution has been found yet (objval == +Inf)

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

/// runs the heuristic selected by DCRsetHeur() (ERA-I by default)
/** Dispatches to DCRheurERAI() or DCRheurERAH() according to the current
 * value of heur, then reports OK if a feasible (mixed-)integer solution was
 * found (UBsol == true) or Infeasible otherwise. */

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

// returns the cost of the best feasible path found by DCRsolve()

double DCR_SPT::DCRgetObj( void )
{
	return objval;
}


/*--------------------------------------------------------------------------*/

// always returns 1 (a single path) if a solution was found, since this
// class only solves *single-path* SFSP instances

int DCR_SPT::DCRgetUBSolNPaths( int k )
{
	if(k < 0 || k>= numFlows)
		throw( DCR::DCRException( "DCR_SPT::DCRgetUBSolNPaths(): wrong flow index." ) );
	if(UBsol) return 1;
	else
		throw( DCR::DCRException( "DCR_SPT::DCRgetUBSolNPaths(): no solution found." ) );	//FIXME: is this okay???
}

/*--------------------------------------------------------------------------*/

// copies the best path found (arcs Xsol[] and corresponding rates Rsol[])
// into the caller-provided arrays X and R, and returns the number of hops

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

// loads a new SFSP DCR instance: only the SRP delay formula and single-flow
// (nflows == 1) instances are supported (an exception is thrown otherwise);
// deep-copies the topology/link/node/flow data and allocates Xsol/Rsol with
// enough room for a Hamiltonian path (at most nnodes - 1 hops)

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

// selects which heuristic DCRsolve() will run: '1' for ERA-I, '2' for ERA-H

void DCR_SPT::DCRsetHeur(char h)
{
	if(h == '1') heur = 1;
	else if(h == '2') heur = 2;
	else throw( DCR::DCRException( "DCR_SPT::DCRsetHeur(): this heur. type is not supported yet..sorry" ) );
}

/*--------------------------------------------------------------------------*/
/*----------------------------------MODIFIERS-------------------------------*/
/*--------------------------------------------------------------------------*/

// closes the arcs listed in whch[] (na of them) for every flow, by setting
// their cost to +Infinity, effectively excluding them from any future path

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

// closes the arcs listed in whch[] (na of them) for flow k only

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

// releases all dynamically allocated memory

DCR_SPT::~DCR_SPT()
{

	clean_up();
}

/*--------------------------------------------------------------------------*/
/*------------------------------ LOCAL (PRIVATE) METHODS -------------------*/
/*--------------------------------------------------------------------------*/

/* This private method implements ERA-I heuristic algorithm
 *  for SFSP DCR problems, see Orda paper.*/
/** ERA-I ("Extended Routing Algorithm - Individual rates") heuristic.
 *
 * The candidate values for the reserved rate r_min are the distinct
 * capacities of the flow's links (Flows[0].caps[]); for each candidate
 * rmin (in increasing order, skipped if below the flow's minimum rate),
 * the "reduced graph" containing only the arcs with capacity >= rmin is
 * built implicitly (by simply ignoring the other arcs while relaxing),
 * and a Shortest Path Tree from the source s is computed on it via a
 * FIFO-queue label-correcting algorithm (a Bellman-Ford variant: a node is
 * (re-)enqueued whenever its distance label is improved, and is not
 * re-inserted while already in the queue) using, as the length of each
 * arc, its contribution MTU/rmin + MTU/speed + node_delay + link_delay to
 * the end-to-end transmission delay for that specific candidate rmin (note
 * arc costs are always nonnegative here, so a simpler Dijkstra could also
 * be used, but the FIFO/Bellman-Ford scheme handles the general case
 * uniformly, consistently with SPT::Solve()).
 *
 * For each candidate rmin, once the tree is computed, the resulting s-t
 * path is checked for feasibility of the deadline constraint (transmission
 * delay along the path, plus the burst term MTU*burst/rmin, must not
 * exceed Flow.deadline); if feasible, its routing cost is computed by
 * pricing every arc of the path at its *own* individual capacity (i.e.,
 * Flows[0].costs[a] * Flows[0].caps[a] for every arc a on the path -- this
 * is the "Individual rates" of ERA-I), and, if cheaper than the best
 * solution found so far over all candidate values of rmin, it replaces it
 * (Xsol/Rsol/nhops/objval are updated accordingly). */
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
/** ERA-H ("Extended Routing Algorithm - Homogeneous/Hop-optimized rates")
 * heuristic.
 *
 * Like DCRheurERAI(), for each candidate rmin (again the distinct link
 * capacities of the flow, in increasing order) a FIFO-queue label-
 * correcting Shortest Path computation is run from the source s over the
 * (implicit) reduced graph of arcs with capacity >= rmin, using the same
 * per-arc transmission-delay length MTU/rmin + MTU/speed + node_delay +
 * link_delay. The key difference with ERA-I is *when* and *how* a
 * candidate path is evaluated: here the check is performed every time the
 * sink t is extracted from the FIFO queue (i.e., as soon as its distance
 * label distance[t] is finalized for the current partial tree, which may
 * happen several times as the label-correcting algorithm improves it),
 * rather than only once at the end.
 *
 * Whenever the current s-t path (backtracked via previous[]) satisfies the
 * deadline constraint using rmin, ERA-H does not price the path's arcs at
 * their individual capacities as ERA-I does; instead it computes, in
 * closed form, the *smallest common rate* r0 that all arcs on the path
 * could share while still meeting the deadline, given the path's fixed
 * (non-transmission) delay dl and number of hops nh:
 * r0 = max( Flow.rate , ( burst + nh * MTU ) / ( deadline - dl ) ). Pricing
 * every arc of the path at this common rate r0 (Flows[0].costs[a] * r0)
 * typically yields a cheaper, still-feasible solution than using each
 * arc's own capacity; as usual, the best (cheapest) such path found over
 * all candidate values of rmin replaces the incumbent (Xsol/Rsol/nhops/
 * objval, and the local r0/cost). */
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
