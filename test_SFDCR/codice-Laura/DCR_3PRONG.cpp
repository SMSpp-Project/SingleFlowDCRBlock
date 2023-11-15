/*--------------------------------------------------------------------------*/
/*---------------------------- File DCR_3PRONG.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------------- IMPLEMENTATION -------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "DCR_3PRONG.h"
#include <fstream>
#include <cstdlib>
#include <set>
#include <queue>

/*--------------------------------------------------------------------------*/
/*----------------------- IMPLEMENTATION OF DCR_3PRONG ---------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------------- PUBLIC METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*----------------------------- CONSTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

//#define DEBUG

DCR_3PRONG::DCR_3PRONG() : DCR()
{
	
	Nodes = 0;
	Links = 0;
	Flows = 0;
	
	Xsol = 0;
	Rsol = 0;	
	
	nhops = -1;
	
	objval = Inf<double>();
	
	socp = 0;
	eraI = 0;
	eraH = 0;
}

/*--------------------------------------------------------------------------*/
/*----------------------------------SOLVE-----------------------------------*/
/*--------------------------------------------------------------------------*/

//TODO
DCR::DCRStatus DCR_3PRONG::DCRsolve( void )
{
	int status;
	
	UBsol = false;

	//try ERA_I first
	eraI = new DCR_SPT;
	eraI->DCRloadProblem(numNodes, numLinks, numFlows, Flows, Links, Nodes, MTU, dtype);
	eraI->DCRsetHeur('1');
	status = eraI->DCRsolve();

	//if ERA_I infeasible stop
	if(status)	
	{
		delete eraI;
		return Infeasible;		 
	}
	
	delete eraI;
	
	//else try ERA_H
	eraH = new DCR_SPT;
	eraH->DCRloadProblem(numNodes, numLinks, numFlows, Flows, Links, Nodes, MTU, dtype);
	eraH->DCRsetHeur('2');
	status = eraH->DCRsolve(); 

	//if eraH is feasible stop
	if(!status)
	{
		objval = eraH->DCRgetObj();
		nhops = eraH->DCRgetUBSolPath(0,0,Xsol,Rsol);
		UBsol = true;
		delete eraH;
		return OK;						
	}
		
	delete eraH;
			
	//else try SOCP
	socp = new DCR_MFSP_SOCP_CPX;
	socp->DCRloadProblem(numNodes, numLinks, numFlows, Flows, Links, Nodes, MTU, dtype);
	socp->DCRsetFsbEps(1e-06);		
	socp->DCRsetOptEps(1e-04);	
	status = socp->DCRsolve();
			
	//if socp is feasible stop
	if(!status)
	{
		objval = socp->DCRgetObj();
		nhops = socp->DCRgetUBSolPath(0,0,Xsol,Rsol);
		UBsol = true;
		delete socp;
		return OK;						
	}
		
	delete socp;
	
	return Infeasible;

}	

/*--------------------------------------------------------------------------*/
/*----------------------------------GET RESULTS-----------------------------*/
/*--------------------------------------------------------------------------*/

double DCR_3PRONG::DCRgetObj( void )
{
	return objval;
}

/*--------------------------------------------------------------------------*/

int DCR_3PRONG::DCRgetUBSolNPaths( int k )
{
	if(k < 0 || k>= numFlows)
		throw( DCR::DCRException( "DCR_INDI::DCRgetUBSolNPaths(): wrong flow index." ) );
	if(UBsol) return 1;
	throw( DCR::DCRException( "DCR_INDI::DCRgetUBSolNPaths(): no solution found." ) );	//FIXME: is this okay???
}

/*--------------------------------------------------------------------------*/

int DCR_3PRONG::DCRgetUBSolPath( int k, int p, int *X, double *R)
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

int DCR_3PRONG::DCRgetPSolNPaths( int k )
{
	throw( DCR::DCRException( "DCR_INDI::DCRgetPSolNPaths(): Continuous relaxation is not supported...sorry." ) );	
}

/*--------------------------------------------------------------------------*/

void DCR_3PRONG::DCRgetPSolPath(int k, int p, double *X, double *R)
{
	throw( DCR::DCRException( "DCR_INDI::DCRgetPSolPath(): Continuous relaxation is not supported...sorry." ) );
}

/*--------------------------------------------------------------------------*/
/*----------------------------------LOAD DATA ------------------------------*/
/*--------------------------------------------------------------------------*/

void DCR_3PRONG::DCRloadProblem(int nnodes, int nlinks, int nflows, 
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
/*----------------------------------MODIFIERS-------------------------------*/
/*--------------------------------------------------------------------------*/

void DCR_3PRONG::DCRcloseArcs( int * whch, int na )
{
	int k, j;
	
	for(k = 0; k < numFlows; k++)
	{
			for(j = 0; j < na; j++)
				Flows[k].costs[whch[j]] = Inf<double>();
	} 
}

/*--------------------------------------------------------------------------*/

void DCR_3PRONG::DCRcloseArcs( int k , int * whch, int na )
{
	int j;
	
	for(j = 0; j < na; j++)
		Flows[k].costs[whch[j]] = Inf<double>();
}

/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

DCR_3PRONG::~DCR_3PRONG()
{

	clean_up();
}

/*--------------------------------------------------------------------------*/
/*------------------------------ LOCAL (PRIVATE) METHODS -------------------*/
/*--------------------------------------------------------------------------*/


/* This private method returns the index of the link 
 * having the given startnode (from) and endnode (to). */
int DCR_3PRONG::DCRgetLink(int from, int to)
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
void DCR_3PRONG::clean_up()
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
void DCR_3PRONG::copyDataArrays(DCRFlow *flows, DCRLink *links, DCRNode *nodes)
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
/*-------------------------- End File DCR_3PRONG.cpp -------------------------*/
/*--------------------------------------------------------------------------*/
