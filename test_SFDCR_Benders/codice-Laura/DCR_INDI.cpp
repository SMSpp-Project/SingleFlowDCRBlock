/*--------------------------------------------------------------------------*/
/*---------------------------- File DCR_INDI.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------------- IMPLEMENTATION -------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "DCR_INDI.h"
#include "DCR_MFSP_SOCP_CPX.h"//FIXME: and GUROBI ???
#include <fstream>
#include <cstdlib>
#include <set>
#include <queue>

/*--------------------------------------------------------------------------*/
/*----------------------- IMPLEMENTATION OF DCR_INDI -----------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------------- PUBLIC METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*----------------------------- CONSTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

//#define DEBUG

DCR_INDI::DCR_INDI() : DCR()
{
	
	Nodes = 0;
	Links = 0;
	Flows = 0;
	
	Xsol = 0;
	Rsol = 0;	
	
	nhops = -1;
	
	heur = 1;
	
	objval = Inf<double>();
	
	socp = 0;
}

/*--------------------------------------------------------------------------*/
/*----------------------------------SOLVE-----------------------------------*/
/*--------------------------------------------------------------------------*/

DCR::DCRStatus DCR_INDI::DCRsolve( void )
{
	if(heur == 1)	
	{
		DCRheurWS();//widest-shortest
		//debug
		//cout << "\nINDI solver option Widest-Shortest" << endl;
	}
	if(heur == 2)
	{
		DCRheurSW();//shortest-widest
		//debug
		//cout << "\nINDI solver option Shortest-Widest" << endl;
	}
	
	if(UBsol) return OK;
	
	return Infeasible;
}


/*--------------------------------------------------------------------------*/
/*----------------------------------GET RESULTS-----------------------------*/
/*--------------------------------------------------------------------------*/

double DCR_INDI::DCRgetObj( void )
{
	return objval;
}

/*--------------------------------------------------------------------------*/

int DCR_INDI::DCRgetUBSolNPaths( int k )
{
	if(k < 0 || k>= numFlows)
		throw( DCR::DCRException( "DCR_INDI::DCRgetUBSolNPaths(): wrong flow index." ) );
	if(UBsol) return 1;
	throw( DCR::DCRException( "DCR_INDI::DCRgetUBSolNPaths(): no solution found." ) );	//FIXME: is this okay???
}

/*--------------------------------------------------------------------------*/

int DCR_INDI::DCRgetUBSolPath( int k, int p, int *X, double *R)
{
	int i;
	
	if(k < 0 || k>= numFlows)
		throw( DCR::DCRException( "DCR_INDI::DCRgetUBSolPath(): wrong flow index." ) );
	if(p != 0)
		throw( DCR::DCRException( "DCR_INDI::DCRgetUBSolPath(): wrong path index." ) );
	if(UBsol)
	{
	
		for(i = 0; i < nhops; i++)
		{
			X[i] = Xsol[i];
			R[i] = Rsol[i];
		}	
		return nhops;	
	}
	throw( DCR::DCRException( "DCR_INDI::DCRgetUBSolPath(): no solution found." ) );//FIXME: is this okay ???		
	
}

/*--------------------------------------------------------------------------*/

int DCR_INDI::DCRgetPSolNPaths( int k )
{
	throw( DCR::DCRException( "DCR_INDI::DCRgetPSolNPaths(): Continuous relaxation is not supported...sorry." ) );	
}

/*--------------------------------------------------------------------------*/

void DCR_INDI::DCRgetPSolPath(int k, int p, double *X, double *R)
{
	throw( DCR::DCRException( "DCR_INDI::DCRgetPSolPath(): Continuous relaxation is not supported...sorry." ) );
}

/*--------------------------------------------------------------------------*/
/*----------------------------------LOAD DATA ------------------------------*/
/*--------------------------------------------------------------------------*/

void DCR_INDI::DCRloadProblem(int nnodes, int nlinks, int nflows, 
	DCRFlow *flows, DCRLink *links, DCRNode *nodes, double mtu, DCRDelay deltype)
{
	int status;
	
	if(deltype != SRP)
		throw( DCR::DCRException( "DCR_INDI::DCRloadProblem(): this delay formula is not supported yet..sorry" ) );
	
	if(nflows > 1)
		throw( DCR::DCRException( "DCR_INDI::DCRloadProblem(): this class can only solve SFSP DCR problems..sorry" ) );
	
	///////////////////////////////
	//INIT
	/////////////////////////////
	clean_up();	
	
	UBsol = false;
	objval = Inf<double>();
	nhops = -1;
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

void DCR_INDI::DCRsetHeur(char h)
{
	if(h == '1') heur = 1;
	else if(h == '2') heur = 2;
	else throw( DCR::DCRException( "DCR_INDI::DCRsetHeur(): this heur. type is not supported yet..sorry" ) );
}

/*--------------------------------------------------------------------------*/
/*----------------------------------MODIFIERS-------------------------------*/
/*--------------------------------------------------------------------------*/

void DCR_INDI::DCRcloseArcs( int * whch, int na )
{
	int k, j;
	
	for(k = 0; k < numFlows; k++)
	{
			for(j = 0; j < na; j++)
				Flows[k].costs[whch[j]] = Inf<double>();
	} 
}

/*--------------------------------------------------------------------------*/

void DCR_INDI::DCRcloseArcs( int k , int * whch, int na )
{
	int j;
	
	for(j = 0; j < na; j++)
		Flows[k].costs[whch[j]] = Inf<double>();
}

/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

DCR_INDI::~DCR_INDI()
{

	clean_up();
}

/*--------------------------------------------------------------------------*/
/*------------------------------ LOCAL (PRIVATE) METHODS -------------------*/
/*--------------------------------------------------------------------------*/

/*Widest-shortest based heuristic*/
void DCR_INDI::DCRheurWS()
{
	
	INDIWSPhase1();
	
	if(nhops < 0) {UBsol = false; return;}
		
	INDIPhase2();
	
}

/*--------------------------------------------------------------------------*/

/*Shortest-widest based heuristic*/
void DCR_INDI::DCRheurSW()
{
	INDISWPhase1();
	
	if(nhops <  0) {UBsol = false; return;}
	
	INDIPhase2();
}

/*--------------------------------------------------------------------------*/

int DCR_INDI::WidestShortest(int maxhops, int s, int t, int *previous)
{
	int j, h, next;
	double rmin;	
	bool found = false;	
	int hopcount = -1;
	int *d = new int[numNodes];
	std::queue<int> myqueue;
	std::set<double> myset;
	std::set<double>::reverse_iterator rit;	

	/****************************************************************/	
	//create a set with all (individual) capacity values
	for(j = 0; j < numLinks; j++)
		myset.insert(Flows[0].caps[j]);
		
	/****************************************************************/	
	//consider capacity values of the set in *decreasing* order	
    for (rit=myset.rbegin(); rit!=myset.rend(); ++rit)
    {
		rmin = *rit;
		
		//check if min capacity >= flow rate
		if(rmin < Flows[0].rate) continue;
	
		//init array previous 	
		for(j = 0; j < numNodes; j++)
			previous[j] = -1;		
		
		//empty queue
		while(! myqueue.empty() )
			myqueue.pop();
					
		//insert source node in queue
		myqueue.push(s);
		
		d[s] = 0;
	
		while (! myqueue.empty())
		{			
			//extract front node from queue 
			h = myqueue.front();
			myqueue.pop();
			
			if(h == t) 
			{
				if(d[t] <= maxhops) { found = true; break; }
			}
			
			//consider FS(h)
			for(j = 0; j < numLinks; j++)
			{								
				if(Links[j].startnode == h) 
				{
					if(Flows[0].caps[j] < rmin) continue;
					
					//get end-node
					next = Links[j].endnode;				
					
					//if not visited, set previous for end-node
					if(previous[next] == -1) 
					{	
						previous[next] = h;
						d[next] = d[h] + 1;
						myqueue.push(next);						
					}//if
					
				}//if
				
			}//for FS(h) 
		
		}//while
		
		if(found) break;
		
	}//for (rit)s
	
	/****************************************************************/
	//returns -1 if no path was found, otherwise returns #hops in the path
	if (found) hopcount = d[t];
	
	delete [] d;

	return hopcount;
}

/*--------------------------------------------------------------------------*/

int DCR_INDI::HopShortest(int s, int t, int *previous)
{
	int j, h, next;
	double rmin;	
	int *d = new int[numNodes];
	bool found = false;
	int hopcount = -1;
	std::queue<int> myqueue;
	std::set<double> myset;
	std::set<double>::reverse_iterator rit;	

	/****************************************************************/
	//init array previous 	
	for(j = 0; j < numNodes; j++)
		previous[j] = -1;	
			
	//insert source node in queue
	myqueue.push(s);

	d[s] = 0;

	while (! myqueue.empty())
	{			
		//extract front node from queue 
		h = myqueue.front();
		myqueue.pop();
			
		//asa sink is extracted from the queue, the hop-shortest path has been found 
		if(h == t) {found = true; break;}
			
		//consider FS(h)
		for(j = 0; j < numLinks; j++)
		{								
			if(Links[j].startnode == h) 
			{
				if(Flows[0].caps[j] < Flows[0].rate) continue;
					
				//get end-node
				next = Links[j].endnode;				
					
				//if not visited, set previous for end-node
				if(previous[next] == -1) 
				{	
					previous[next] = h;
					d[next] = d[h] + 1;
					myqueue.push(next);						
				}//if
					
			}//if
				
		}//for FS(h) 
		
	}//while
	
	/****************************************************************/
	//returns -1 if no path was found, otherwise returns #hops in the path
	if (found) hopcount = d[t];
	
	delete [] d;

	return hopcount;	
}

/*--------------------------------------------------------------------------*/

void DCR_INDI::INDIWSPhase1()
{
	int i, j, h;	
	int *previous = new int[numNodes];
	
	nhops = WidestShortest(numNodes, Flows[0].sourcenode, Flows[0].sinknode, previous);
	
	if(nhops < 0) { delete [] previous; return; }
	
	i = Flows[0].sinknode;
	j = 0;
	while( i != Flows[0].sourcenode )
	{
		h = DCRgetLink(previous[i], i);
		Xsol[j] = h;
		i = previous[i];
		j++;
	}	
	
	delete [] previous;	
}

/*--------------------------------------------------------------------------*/

void DCR_INDI::INDISWPhase1()
{
	int i, j, h;

	int *previous = new int[numNodes];
	
	nhops = HopShortest(Flows[0].sourcenode, Flows[0].sinknode, previous);
	
	if(nhops < 0) { delete [] previous; return; }
	
	WidestShortest(nhops, Flows[0].sourcenode, Flows[0].sinknode, previous);
	
	i = Flows[0].sinknode;
	j = 0;
	while( i != Flows[0].sourcenode )
	{
		h = DCRgetLink(previous[i], i);
		Xsol[j] = h;
		i = previous[i];
		j++;		
	}	
		
	delete [] previous;	
}

/*--------------------------------------------------------------------------*/

/*Compute reservable rates with fixed path*/
void DCR_INDI::INDIPhase2()
{
	int i, c, snh;
	double *r, cost;
	int *toclose, *x, nc;
	bool *app;
	int status;
	
	//create an array with the arcs that must be closed	
	app = new bool[numLinks];
	for(i = 0; i < numLinks; i++)
		app[i] = true;
	for(i = 0; i < nhops; i++)
		app[Xsol[i]] = false;
	nc = numLinks - nhops;
	toclose = new int[nc];
	c = 0;
	for(i = 0; i < numLinks ; i++)
	{
		if(app[i]) 
		{
			toclose[c]=i; 
			c++;
		}
	}
	
	x = new int[numNodes-1];
	r = new double[numNodes -1];
	socp = new DCR_MFSP_SOCP_CPX;
	
	socp->DCRloadProblem(numNodes, numLinks, numFlows, Flows, Links, Nodes, MTU, dtype);
	
	socp->DCRsetFsbEps(1e-06);		
	socp->DCRsetOptEps(1e-04);
	socp->DCRcloseArcs(toclose, nc);	
	
	status = socp->DCRsolve();			
	if (status)
	{
		UBsol = false;
		
		delete [] app;
		delete [] x;
		delete [] r;
		delete socp;
		
		return;		
	}
		
	snh = socp->DCRgetUBSolPath(0,0,x,r);	

	if(snh != nhops)
		throw( DCR::DCRException( "DCR_INDI::INDI_Phase2(): internal error: wrong num hops" ) );
		
	cost = 0;
	for(i = 0; i < nhops; i++)
	{
			Xsol[i] = x[i];
			Rsol[i] = r[i];	
			cost += Flows[0].costs[Xsol[i]]*Rsol[i];	
	}
		
	UBsol = true;
	
	//get objval (if feasible)
	objval = cost;
	
	delete [] app;
	delete [] x;
	delete [] r;
	delete socp;	
}

/*--------------------------------------------------------------------------*/

/* This private method returns the index of the link 
 * having the given startnode (from) and endnode (to). */
int DCR_INDI::DCRgetLink(int from, int to)
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
void DCR_INDI::clean_up()
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
void DCR_INDI::copyDataArrays(DCRFlow *flows, DCRLink *links, DCRNode *nodes)
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
/*-------------------------- End File DCR_INDI.cpp -------------------------*/
/*--------------------------------------------------------------------------*/
