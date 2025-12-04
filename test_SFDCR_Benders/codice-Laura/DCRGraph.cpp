/*--------------------------------------------------------------------------*/
/*---------------------------- File DCRGraph.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------------- IMPLEMENTATION -------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "DCRGraph.h"

#include <string.h>
#include <fstream>


/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

//FIXME: what is this for ??? do i need it ???

/*#if( OPT_USE_NAMESPACES )
 using namespace MMCFGraph_di_unipi_it;
#else
 using namespace std;
#endif*/

/*--------------------------------------------------------------------------*/
/*----------------------- IMPLEMENTATION OF DCRGraph  ----------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------- PUBLIC METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

/*----------------------------- CONSTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

DCRGraph::DCRGraph( const char *const FN , char FT , double myperccap, double mypercdelay, int myseed) : Graph(FN, FT)
{
	//randgen= new OPTrand();
	perc_cap = myperccap;
	perc_delay = mypercdelay;
	seed = myseed;
	//randgen->srand(seed);

	//init arrays to 0
	NodeDelays = LinkDelays = FlowBursts = FlowDeadlines = 0;
	
	//read DCR part: .dcr file
 	int l = strlen( FN );
	char *Name = new char[ l + 5 ];  // temporary string containing the constant
	strcpy( Name , FN );             // part of the pathname + space for `.dcr'
	strcpy( Name + l , ".dcr" );

	ifstream inFile( Name );
	if( ! inFile.is_open() )
		throw( DCR::DCRException( "DCRGraph::DCRGraph(): cannot open .dcr file" ) );
	
	//allocate memory	
	int i;
	try
	{
		NodeDelays = new double[NNodes];
		LinkDelays = new double[NArcs];
		FlowBursts = new double[NComm];
		FlowDeadlines = new double[NComm];
	}
	catch(exception& e)
	{
		cout << "Standard exception: " << e.what() << endl;
	}
	
	//read node delays
	for(i = 0; i < NNodes; i++)
		inFile >> NodeDelays[i];
		
	//read link delays
	for(i = 0; i < NArcs; i++)
		inFile >> LinkDelays[i];
		
	//read flow burst and deadlines
	for(i = 0; i < NComm; i++)
	{
		inFile >> FlowBursts[i];
		inFile >> FlowDeadlines[i];
	}
	
	//read MTU;
	inFile >> MTU;

	inFile.close();
	
	delete [] Name;

}  // end( DCRGraph( char* , char ) )

/*--------------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/

double DCRGraph::getSupply(int Kflow)
{
	int j;
	double def, supply = 0;
	int single = 0;
	
	for(j = 0; j < NNodes; j++)
	{
		def = DeficitKJ( Kflow , j );
			
		if( def > 0 ) 
		{ supply = def; single ++; }	
	}
	
	if( supply == 0 ) 
	{ throw( DCR::DCRException( "DCRGraph::getSupply(): 0 supply" ) ); }
	
	//NB. we want O-D instances, so **only one positive supply** must exist
	if( single > 1 ) throw( DCR::DCRException( "DCRGraph::getSupply(): multiple positive supplies found" ) );
	
	
	return supply;	
}

/*--------------------------------------------------------------------------*/

int DCRGraph::getSource(int Kflow)
{
	int j, source = -1;
	double def;
	int single = 0;
	
	for(j = 0; j < NNodes; j++)
	{
		def = DeficitKJ( Kflow, j );
		
		//NB. source has a negative deficit, i.e., a positive supply 
		if( def < 0 ) 
		{ source = j; single++; }
	}
	
	if( source == -1 ) throw( DCR::DCRException( "DCRGraph::getSource(): no source found" ) );
	
	//NB. we want O-D instances, so **only one source** must exist
	if( single > 1 ) throw( DCR::DCRException( "DCRGraph::getSource(): multiple sources found" ) );
	
	return source;
}

/*--------------------------------------------------------------------------*/

int DCRGraph::getSink(int Kflow)
{
	int j, sink = -1;
	double def;
	int single = 0;
	
	for(j = 0; j < NNodes; j++)
	{
		def = DeficitKJ( Kflow , j );
		
		//NB. destination has a positive deficit, i.e., a negative supply
		if ( def  > 0 ) 
		{ sink = j; single++; }
	}
	
	if( sink == -1 ) throw( DCR::DCRException( "DCRGraph::getSink(): no sink found" ) );
	
	//NB. we want O-D instances, so **only one sink** must exist
	if( single > 1 ) throw( DCR::DCRException( "DCRGraph::getSink(): multiple sinks found" ) );
	
	return sink;
}	

/*--------------------------------------------------------------------------*/
 
DCR::DCRFlow * DCRGraph::getFlows()
{ 
	int k,j;
	DCR::DCRFlow *Flows = 0;

	try
	{
		Flows = new DCR::DCRFlow[NComm];
	}
	catch(exception& e)
	{ cout << "Standard exception: " << e.what() << endl; }
	
	for(k = 0; k < NComm; k++)
	{
		Flows[k].sourcenode = getSource(k); 
		Flows[k].sinknode = getSink(k); 
		Flows[k].burst =  FlowBursts[k];  
		Flows[k].rate = getSupply(k); 
		Flows[k].deadline = FlowDeadlines[k]; //* ( 1 + perc_delay * randgen->rand() ); 
		 
		try
		{
		Flows[k].costs = new double[NArcs];
		}
		catch(exception& e)
		{ cout << "Standard exception: " << e.what() << endl; }
		
		for(j = 0; j < NArcs; j++)
			Flows[k].costs[j] = CostKJ(k,j);
		
		//FIXME: how to deal with null capacities ???
		try
		{
		Flows[k].caps = new double[NArcs];
		}
		catch(exception& e)
		{ cout << "Standard exception: " << e.what() << endl; }
		
		for(j = 0; j < NArcs; j++)
			Flows[k].caps[j] = CapacityKJ(k,j); //* ( 1 - perc_cap * randgen->rand() );
		
	}
	
	return Flows;
	
 }
 
/*--------------------------------------------------------------------------*/

DCR::DCRLink * DCRGraph::getLinks()
{
	int j;
	DCR::DCRLink *Links = 0;
	
	try
	{
		Links = new DCR::DCRLink[NArcs];
	}
	catch(exception& e)
	{cout << "Standard exception: " << e.what() << endl;}
	
	for(j = 0; j < NArcs; j++)
	{
		Links[j].speed =  Links[j].capacity = TotalCapacityJ( j ); 
		Links[j].delay = LinkDelays[j]; 
		Links[j].startnode = (int) Startn[j] -1;
		Links[j].endnode = (int) Endn[j] -1;
	}	
	
	return Links;	
}

/*--------------------------------------------------------------------------*/

DCR::DCRNode * DCRGraph::getNodes()
{ 
	int i;
	DCR::DCRNode *Nodes = 0;
	
	try
	{
		Nodes = new DCR::DCRNode[NNodes];
	} 
	catch(exception& e)
	{ cout << "Standard exception: " << e.what() << endl; }
	
	for(i = 0; i < NNodes; i++)
	{
		Nodes[i].delay = NodeDelays[i]; 
	}
	
	return Nodes;
}

/*--------------------------------------------------------------------------*/

int DCRGraph::getNumNodes()
{ return NNodes; }

/*--------------------------------------------------------------------------*/

int DCRGraph::getNumLinks()
{ return NArcs; }

/*--------------------------------------------------------------------------*/

int DCRGraph::getNumFlows()
{ return NComm; }

/*--------------------------------------------------------------------------*/

double DCRGraph::getMTU()
{ return MTU; }

/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

DCRGraph::~DCRGraph()
{
	delete [] NodeDelays;
	delete [] LinkDelays;
	delete [] FlowBursts;
	delete [] FlowDeadlines;
}  


/*--------------------------------------------------------------------------*/
/*-------------------------- End File DCRGraph.cpp -------------------------*/
/*--------------------------------------------------------------------------*/
