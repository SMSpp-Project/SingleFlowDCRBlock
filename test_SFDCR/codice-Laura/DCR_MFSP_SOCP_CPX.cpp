/*--------------------------------------------------------------------------*/
/*---------------------------- File DCR_MFSP_SOCP_CPX.cpp ------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------------- IMPLEMENTATION -------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "DCR_MFSP_SOCP_CPX.h"
#include <fstream>
#include <cstdlib>
#include <math.h>

/*--------------------------------------------------------------------------*/
/*----------------------- IMPLEMENTATION OF DCR_MFSP_SOCP_CPX --------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------------- PUBLIC METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*----------------------------- CONSTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

//#define DEBUG
#define VARNAMES

DCR_MFSP_SOCP_CPX::DCR_MFSP_SOCP_CPX() : DCR()
{

	//init all class members
	Nodes = 0;
	Links = 0;
	Flows = 0;
		
	env = 0;
	lp = 0;
	
	//FIXME: is this okay ???
	flagBigM = 0;
	
	r_offset = 0;
	rmin_offset = 0;
	s_offset = 0;
	t_offset = 0;	
	rp_offset = 0;
	theta_offset = 0;
	/* FOFFSET=#vars for a (any) *single flow* of the given network.
	 * NB. The number of variables only depends on the network topology and 
	 * it is independent of the flow considered.
	 * This offset is used to distinguish variables of one flow from the other
	 * in the multi-flow case, where variables are ordered one flow after the other.
	 * NB. There are no variables in common to all flows.*/
	FOFFSET = 0;
	
	Xsol = 0;
	Rsol = 0;	
	
	//FIXME: is this okay ???
	objval = Inf<double>();
}

/*--------------------------------------------------------------------------*/
/*----------------------------------WRITE LP--------------------------------*/
/*--------------------------------------------------------------------------*/

void DCR_MFSP_SOCP_CPX::DCRwrite( char *filename )
{
	int status;
	
	status = CPXwriteprob (env, lp, filename, NULL);
	if ( status ) 
			throw( DCR::DCRException( "DCR_MFSP_SOCP_CPX::DCRwrite(): failed to write CPLEX lp file" ) );
}

/*--------------------------------------------------------------------------*/
/*----------------------------------SOLVE-----------------------------------*/
/*--------------------------------------------------------------------------*/

DCR::DCRStatus DCR_MFSP_SOCP_CPX::DCRsolve( void )
{
	int status;
	int probtype;
	int lpstat;
	
	//debug
	//cout << "\nSOCP CPX solver " << endl;
		
	probtype = CPXPROB_MIQCP;
	
	status = CPXsetdblparam (env, CPX_PARAM_EPAGAP, 1e-04);
	//if ( status ) 
	//	throw( DCR::DCRException( "DCR_MFSP_SOCP_CPX::DCRsetOptEps(): failed to set CPLEX absolute MIP gap tolerance" ) );
	

	//Mixed-Integer case
	//////////////////////////////////////
	if(probtype == CPXPROB_MIQCP)
	{	
		status = CPXmipopt (env, lp);
		if( status )
		{
			//throw( DCR::DCRException( "DCR_MFSP_SOCP_CPX::DCRsolve(): failed to solve problem" ) );
			UBsol = false;
			return Error;
		}

		lpstat = CPXgetstat(env, lp);

		//FIXME: more CPX codes...
		if( (lpstat == CPXMIP_OPTIMAL) || 
		    (lpstat == CPXMIP_OPTIMAL_TOL)
		  ) { getSol(); UBsol=true; return OK; }
		  
		UBsol=false;
		
		if(lpstat == CPXMIP_INFEASIBLE) return Infeasible;
		
		return Error;
	}
	
	//Continuous case
	///////////////////////////////////
	if(probtype == CPXPROB_QCP)
	{
		status = CPXbaropt(env, lp);
		if( status )
			throw( DCR::DCRException( "DCR_MFSP_SOCP_CPX::DCRsolve(): failed to solve problem LP relax." ) );
		
		lpstat = CPXgetstat(env, lp);
		
		//FIXME: more CPX codes...
		if(lpstat == CPX_STAT_OPTIMAL) 
		{ getSol(); Psol=true; return OK; }
		
		UBsol = false;
		
		if(lpstat == CPX_STAT_INFEASIBLE) return Infeasible;
		
		return Error;
	}
	
	//wrong type
	//////////////////
	throw( DCR::DCRException( "DCR_MFSP_SOCP_CPX::DCRsolve(): cannot solve problem, CPX probtype is wrong" ) );
	
}

/*--------------------------------------------------------------------------*/
/*----------------------------------READ RESULTS-----------------------------*/
/*--------------------------------------------------------------------------*/

double DCR_MFSP_SOCP_CPX::DCRgetObj( void )
{	
	return objval;
}

/*--------------------------------------------------------------------------*/

int DCR_MFSP_SOCP_CPX::DCRgetNodes( void )
{
	return CPXgetnodecnt(env,lp);
}

/*--------------------------------------------------------------------------*/

int DCR_MFSP_SOCP_CPX::DCRgetUBSolNPaths( int k )
{
	if(k < 0 || k>= numFlows)
		throw( DCR::DCRException( "DCR_MFSP_SOCP_CPX::DCRgetUBSolNPaths(): wrong flow index." ) );
	if(UBsol) return 1;
	else
		throw( DCR::DCRException( "DCR_MFSP_SOCP_CPX::DCRgetUBSolNPaths(): no solution found." ) );	//FIXME: is this okay???
}

/*--------------------------------------------------------------------------*/

int DCR_MFSP_SOCP_CPX::DCRgetUBSolPath( int k, int p, int *X, double *R)
{
	int i, nhops;
	
	
	if(k < 0 || k>= numFlows)
		throw( DCR::DCRException( "DCR_MFSP_SOCP_CPX::DCRgetUBSolPath(): wrong flow index." ) );
	if(p != 0)
		throw( DCR::DCRException( "DCR_MFSP_SOCP_CPX::DCRgetPUBolPath(): wrong path index." ) );
	if(UBsol)
	{
		nhops = 0;
		for(i = 0; i < numLinks; i++)
		{
			if(Xsol[i + k*numLinks] > 0.5)
			{
				X[nhops] = i;
				R[nhops] = Rsol[i + k*numLinks];
				nhops++;
			}
		}		
	}
	else
		throw( DCR::DCRException( "DCR_MFSP_SOCP_CPX::DCRgetUBSolPath(): no solution found." ) );//FIXME: is this okay ???		
	
	return nhops;
}

/*--------------------------------------------------------------------------*/

int DCR_MFSP_SOCP_CPX::DCRgetPSolNPaths( int k )
{
	if(k < 0 || k>= numFlows)
		throw( DCR::DCRException( "DCR_MFSP_SOCP_CPX::DCRgetPSolNPaths(): wrong flow index." ) );
	if(Psol) return 1;
	else
		throw( DCR::DCRException( "DCR_MFSP_SOCP_CPX::DCRgetPSolNPaths(): no solution found." ) );	//FIXME: is this okay???
}

/*--------------------------------------------------------------------------*/

void  DCR_MFSP_SOCP_CPX::DCRgetPSolPath(int k, int p, double *X, double *R)
{
	int i;
	
	if(k < 0 || k>= numFlows)
		throw( DCR::DCRException( "DCR_MFSP_SOCP_CPX::DCRgetPSolPath(): wrong flow index." ) );
	if(p != 0)
		throw( DCR::DCRException( "DCR_MFSP_SOCP_CPX::DCRgetPSolPath(): wrong path index." ) );
	if(Psol)
	{
		for(i = 0; i < numLinks; i++)
		{
			X[i] = Xsol[i + k*numLinks];
			R[i] = Rsol[i + k*numLinks];
		}		
	}
	else
		throw( DCR::DCRException( "DCR_MFSP_SOCP_CPX::DCRgetPSolPath(): no solution found." ) );//FIXME: is this okay ???	
}

/*--------------------------------------------------------------------------*/
/*----------------------------------LOAD DATA ------------------------------*/
/*--------------------------------------------------------------------------*/

void DCR_MFSP_SOCP_CPX::DCRloadProblem(int nnodes, int nlinks, int nflows, 
	DCRFlow *flows, DCRLink *links, DCRNode *nodes, double mtu, DCRDelay deltype)
{
	int k;
	int status;
		
	if(deltype != SRP)
		throw( DCR::DCRException( "DCR_MFSP_SOCP_CPX::DCRloadProblem(): this delay formula is not supported yet..sorry" ) );
	
	
	///////////////////////////////////
	//INIT
	////////////////////////////////// 
	clean_up();	
	
	Psol = false;
	UBsol = false;
	objval = Inf<double>();
	//////////////////////////////////////////////////////////
	
	numNodes = nnodes;
	numLinks = nlinks;
	numFlows = nflows;	
	MTU = mtu;
	dtype = deltype;	
	copyDataArrays(flows, links, nodes);
	Xsol = new double[numLinks*numFlows];
	Rsol = new double[numLinks*numFlows];		


	//open CPLEX env
	env = CPXopenCPLEX(&status);
	if( env == 0 ) 
		throw( DCR::DCRException( "DCR_MFSP_SOCP_CPX::DCRloadProblem(): failed to open CPLEX env" ) );
	
#ifdef DEBUG
	//turn on CPLEX data checking
	status = CPXsetintparam (env, CPX_PARAM_DATACHECK, CPX_ON);
		if ( status ) 
			throw( DCR::DCRException( "DCR_MFSP_SOCP_CPX::DCRloadProblem(): failed to turn on CPLEX data check" ) );

	//turn on CPLEX output
	status = CPXsetintparam (env, CPX_PARAM_SCRIND, CPX_ON);
		if ( status ) 
			throw( DCR::DCRException( "DCR_MFSP_SOCP_CPX::DCRloadProblem(): failed to turn on CPLEX output" ) );
#endif

	//create CPLEX problem
	lp = CPXcreateprob(env, &status, "DCR");
	if ( lp == 0 )
			throw( DCR::DCRException( "DCR_MFSP_SOCP_CPX::DCRloadProblem(): failed to create CPLEX prob" ) );

	//set CPLEX obj sen
	CPXchgobjsen(env, lp, CPX_MIN);

	FOFFSET = 3*numLinks +2;
	if(flagBigM) FOFFSET += numLinks*2;
	
	for(k = 0; k < nflows; k++)
	{
		//set CPLEX cols
		addCols(k);
	
		//set CPLEX rows	
		addFlowRows(k);
		addRateBounds(k);	
		if(flagBigM) addDelayRowsBigM(k);
		else addDelayRows(k);		
	}

	addCapacityRows();
}


/*--------------------------------------------------------------------------*/
/*----------------------------------SETTERS---------------------------------*/
/*--------------------------------------------------------------------------*/

void DCR_MFSP_SOCP_CPX::DCRsetModel(char method)
{
	if(method=='b') flagBigM = 1;	
	else if(method=='p') flagBigM = 0;
	else throw( DCR::DCRException( "DCR_MFSP_SOCP_CPX::DCRsetModel(): failed to set model" ) );
}

void DCR_MFSP_SOCP_CPX::DCRsetOptEps( double OE )
{
	int status;
	
	optEps = OE;
		
	//FIXME: how about this : set absolute MIP gap tolerance ???
	status = CPXsetdblparam (env, CPX_PARAM_EPAGAP, optEps);
	if ( status ) 
		throw( DCR::DCRException( "DCR_MFSP_SOCP_CPX::DCRsetOptEps(): failed to set CPLEX absolute MIP gap tolerance" ) );
	
	
	status = CPXsetdblparam(env,CPX_PARAM_EPOPT,optEps);           // simplex tolerance for optimality
	if ( status ) 
		throw( DCR::DCRException( "DCR_MFSP_SOCP_CPX::DCRsetOptEps(): failed to set CPLEX optimality tolerance" ) );
}

/*--------------------------------------------------------------------------*/

void DCR_MFSP_SOCP_CPX::DCRsetFsbEps( double OE )
{
	int status;
	
	fsbEps = OE;
		
	status = CPXsetdblparam(env,CPX_PARAM_EPRHS,fsbEps);           // simplex tolerance for feasibility
	if ( status ) 
		throw( DCR::DCRException( "DCR_MFSP_SOCP_CPX::DCRsetFsbEps(): failed to set CPLEX feasibility tolerance" ) );
}

/*--------------------------------------------------------------------------*/

void DCR_MFSP_SOCP_CPX::DCRsetCuts( int setval )
{
	int status;
	//setval =-1 -> switch off cplex cuts:
	
	status = CPXsetintparam (env, CPX_PARAM_CLIQUES, setval);
	if ( status ) 
		throw( DCR::DCRException( "DCR_MFSP_SOCP_CPX::DCRsetNoCuts(): failed to switch off CPLEX cuts" ) );
	status = CPXsetintparam (env, CPX_PARAM_COVERS, setval);
	if ( status ) 
		throw( DCR::DCRException( "DCR_MFSP_SOCP_CPX::DCRsetNoCuts(): failed to switch off CPLEX cuts" ) );
	status = CPXsetintparam (env, CPX_PARAM_DISJCUTS, setval);
	if ( status ) 
		throw( DCR::DCRException( "DCR_MFSP_SOCP_CPX::DCRsetNoCuts(): failed to switch off CPLEX cuts" ) );
	status = CPXsetintparam (env, CPX_PARAM_FLOWCOVERS, setval);
	if ( status ) 
		throw( DCR::DCRException( "DCR_MFSP_SOCP_CPX::DCRsetNoCuts(): failed to switch off CPLEX cuts" ) );
	status = CPXsetintparam (env, CPX_PARAM_FLOWPATHS, setval);
	if ( status ) 
		throw( DCR::DCRException( "DCR_MFSP_SOCP_CPX::DCRsetNoCuts(): failed to switch off CPLEX cuts" ) );
	status = CPXsetintparam (env, CPX_PARAM_FRACCUTS, setval);
	if ( status ) 
		throw( DCR::DCRException( "DCR_MFSP_SOCP_CPX::DCRsetNoCuts(): failed to switch off CPLEX cuts" ) );
	status = CPXsetintparam (env, CPX_PARAM_GUBCOVERS, setval);
	if ( status ) 
		throw( DCR::DCRException( "DCR_MFSP_SOCP_CPX::DCRsetNoCuts(): failed to switch off CPLEX cuts" ) );
	status = CPXsetintparam (env, CPX_PARAM_IMPLBD, setval);
	if ( status ) 
		throw( DCR::DCRException( "DCR_MFSP_SOCP_CPX::DCRsetNoCuts(): failed to switch off CPLEX cuts" ) );
	status = CPXsetintparam (env, CPX_PARAM_MCFCUTS, setval);
	if ( status ) 
		throw( DCR::DCRException( "DCR_MFSP_SOCP_CPX::DCRsetNoCuts(): failed to switch off CPLEX cuts" ) );
	status = CPXsetintparam (env, CPX_PARAM_MIRCUTS, setval);
	if ( status ) 
		throw( DCR::DCRException( "DCR_MFSP_SOCP_CPX::DCRsetNoCuts(): failed to switch off CPLEX cuts" ) );
	status = CPXsetintparam (env, CPX_PARAM_ZEROHALFCUTS, setval);
	if ( status ) 
		throw( DCR::DCRException( "DCR_MFSP_SOCP_CPX::DCRsetNoCuts(): failed to switch off CPLEX cuts" ) );	
}

/*--------------------------------------------------------------------------*/

void DCR_MFSP_SOCP_CPX::DCRsetNodeLimit( int nodelimit )
{
	int status;
	
	//nodelimit = 0 -> root node
	status = CPXsetlongparam (env, CPX_PARAM_NODELIM, nodelimit);
	if ( status ) 
		throw( DCR::DCRException( "DCR_MFSP_SOCP_CPX::DCRsetNodeLimit(): failed to set CPLEX  node limit" ) );
}

/*--------------------------------------------------------------------------*/

void DCR_MFSP_SOCP_CPX::DCRsetHeur( int setval )
{
	int status;
	
	//setval = -1 -> switch off heuristics	
	status = CPXsetlongparam(env, CPX_PARAM_HEURFREQ, setval);
	if ( status ) 
		throw( DCR::DCRException( "DCR_MFSP_SOCP_CPX::DCRsetHeur(): failed to set CPLEX heurs " ) );
}

/*--------------------------------------------------------------------------*/

void DCR_MFSP_SOCP_CPX::DCRsetNodeFile(int setval)
{
	int status;
	
	//setval = 3 -> use compressed tree storage on disk
	status = CPXsetintparam(env, CPX_PARAM_NODEFILEIND, setval);
	if ( status ) 
		throw( DCR::DCRException( "DCR_MFSP_SOCP_CPX::DCRsetNoHeur(): failed to set CPLEX node file " ) );
}

/*--------------------------------------------------------------------------*/

void DCR_MFSP_SOCP_CPX::DCRsetQCPType( void )
{
	int i, k, c, status;
	int *indices;
	char *ctype;
	
	try{
		indices = new int[numLinks*numFlows];
		ctype = new char[numLinks*numFlows];
	}catch(exception& e)
	{
		cout << "Standard exception: " << e.what() << endl;
	}
	
	c = 0;
	for(k = 0; k < numFlows; k++)
	{
		for(i = 0; i < numLinks; i++)
		{
			indices[c] = i + k*FOFFSET;
			ctype[c] = 'C';
			c++;
		}
	}
	
	status = CPXchgctype (env, lp, numLinks*numFlows, indices, ctype);
	if( status )
			throw( DCR::DCRException( "DCR_MFSP_SOCP_CPX::DCRsetQCPType(): failed to set CPLEX vars type" ) );
			
	status = CPXchgprobtype (env, lp, CPXPROB_QCP);
	if( status )
			throw( DCR::DCRException( "DCR_MFSP_SOCP_CPX::DCRsetQCPType(): failed to set CPLEX prob type" ) );
	

	delete [] indices;
	delete [] ctype;	
	
		
}
	
/*--------------------------------------------------------------------------*/
void DCR_MFSP_SOCP_CPX::DCRsetMIQCPType( void )
{
	int i, k, c, status;
	int *indices;
	char *ctype;
		
	try{
		indices = new int[numLinks*numFlows];
		ctype = new char[numLinks*numFlows];
	}catch(exception& e)
	{
		cout << "Standard exception: " << e.what() << endl;
	}
	
	c = 0;
	for(k = 0; k < numFlows; k++)
	{
		for(i = 0; i < numLinks; i++)
		{
			indices[c] = i + k*FOFFSET;
			ctype[c] = 'B';
			c++;
		}
	}
	
	status = CPXchgctype (env, lp, numLinks*numFlows, indices, ctype);
	if( status )
			throw( DCR::DCRException( "DCR_MFSP_SOCP_CPX::DCRsetMIQCPType(): failed to set CPLEX vars type" ) );
			
	status = CPXchgprobtype (env, lp, CPXPROB_MIQCP);
	if( status )
			throw( DCR::DCRException( "DCR_MFSP_SOCP_CPX::DCRsetMIQCPType(): failed to set CPLEX prob type" ) );
	
	
	delete [] indices;
	delete [] ctype;	
	
}

/*--------------------------------------------------------------------------*/

void DCR_MFSP_SOCP_CPX::DCRsetLog(char *filename)
{
	int status;
	
	//CPXFILEptr logfileptr;
    //logfileptr = CPXfopen (filename, "w");
	
	status = CPXsetlogfilename (env, filename, "w");
	if( status )
			throw( DCR::DCRException( "DCR_MFSP_SOCP_CPX::DCRsetLog(): failed to set CPLEX log" ) );

}

/*--------------------------------------------------------------------------*/

void DCR_MFSP_SOCP_CPX::DCRsetMIQCPStrat(int setval)
{
	int status;
	
	//0 = auto(let CPLEX decide), 1 = QCP relax., 2 = LP relax.
    status =  CPXsetintparam(env, CPX_PARAM_MIQCPSTRAT, setval);
	if( status )
			throw( DCR::DCRException( "DCR_MFSP_SOCP_CPX::DCRsetMIQCPStrat(): failed to set CPLEX MIQCP strategy" ) );	
}

/*--------------------------------------------------------------------------*/

void DCR_MFSP_SOCP_CPX::DCRsetNodeSelection(int setval)
{
	int status;
	
	//0 = Depth-first search 
	status = CPXsetintparam(env, CPX_PARAM_NODESEL, setval);
	if( status )
		throw( DCR::DCRException( "DCR_MFSP_SOCP_CPX::DCRsetNodeSelection(): failed to set CPLEX MIQCP node selection strategy" ) );	
}
/*--------------------------------------------------------------------------*/

void DCR_MFSP_SOCP_CPX::DCRsetTimeLimit( double setval )
{
	int status;
	
	//time limit in seconds
	status = CPXsetdblparam(env, CPX_PARAM_TILIM, setval);
	if( status )
	throw( DCR::DCRException( "DCR_MFSP_SOCP_CPX::DCRsetTimeLimit(): failed to set CPLEX time limit" ) );	
}

/*--------------------------------------------------------------------------*/
/*----------------------------------MODIFIERS-------------------------------*/
/*--------------------------------------------------------------------------*/

void DCR_MFSP_SOCP_CPX::DCRcloseArcs( int * whch, int na)
{
	int status;
	int j, k, c;
	int cnt;
	int *indices;
	char *lu;
	double *bd;
	
	// close all arcs with indices in whch for all flows
	// this is done using CPLEX function CPXchgbds
	// and setting the upper bound of the corresponding
	// variables to 0. 

	//number of variables whose bound must be changed
	cnt = na * numFlows;
	
	//array containing variable indices
	indices = new int[cnt];
	
	//array containing info on how to change the bound
	//'L' = lower bound must be changed
	//'U' = upper bound must be changed
	//'B' = both (lower and upper) bounds must be changed
	lu = new char[cnt];
	
	//array containing the new bound values
	bd = new double[cnt];
	
	c = 0;
	for(k = 0; k < numFlows; k++)
	{
		for(j = 0; j < na; j++)
		{
			 indices[c] = whch[j] + FOFFSET*k;
			 lu[c] = 'U';
			 bd[c] = 0;
			 c++;
		}		
	}
	
	status = CPXchgbds (env, lp, cnt, indices, lu, bd);
	if( status )
	throw( DCR::DCRException( "DCR_MFSP_SOCP_CPX::DCRcloseArcs(): failed to change CPLEX bounds" ) );	
	
	delete [] indices;
	delete [] lu;
	delete [] bd;	
	
}

void DCR_MFSP_SOCP_CPX::DCRcloseArcs( int k , int * whch, int na )
{
	int status;
	int j, c;
	int cnt;
	int *indices;
	char *lu;
	double *bd;
	
	// close all arcs with indices in whch for flow k
	// this is done using CPLEX function CPXchgbds
	// and setting the upper bound of the corresponding
	// variables to 0. 

	//number of variables whose bound must be changed
	cnt = na;
	
	//array containing variable indices
	indices = new int[cnt];
	
	//array containing info on how to change the bound
	//'L' = lower bound must be changed
	//'U' = upper bound must be changed
	//'B' = both (lower and upper) bounds must be changed
	lu = new char[cnt];
	
	//array containing the new bound values
	bd = new double[cnt];
	
	c = 0;	
	for(j = 0; j < na; j++)
	{
		indices[c] = whch[j] + FOFFSET*k;
		lu[c] = 'U';
		bd[c] = 0;
		c++;
	}
	
	status = CPXchgbds (env, lp, cnt, indices, lu, bd);
	if( status )
	throw( DCR::DCRException( "DCR_MFSP_SOCP_CPX::DCRcloseArcs(): failed to change CPLEX bounds" ) );	
	
	delete [] indices;
	delete [] lu;
	delete [] bd;	
}

/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

DCR_MFSP_SOCP_CPX::~DCR_MFSP_SOCP_CPX()
{
	
	clean_up();

}

/*--------------------------------------------------------------------------*/
/*------------------------------ LOCAL (PRIVATE) METHODS -------------------*/
/*--------------------------------------------------------------------------*/

/* This private method adds columns*/
 /* Number of vars is: 
 * *****************************
 * In case of perspective SOCP model:
 * *****************************
 * The order of the vars is:
 * x_link1, ... ,x_linkm, r_link1, ... ,r_linkm, r_min, s_link1, ..., s_linkm, t
 * In case of bigM formulation we should add the following variables:
 * (rp_ij), (theta_ij)
 * */
void DCR_MFSP_SOCP_CPX::addCols( int flowidx )
{
	
	int i, c;
	int status, ccnt;
	double *obj, *lb, *ub;
	char **vnames = 0;
	char *ctype;

	ccnt = 3*numLinks +2;

if(flagBigM)
	ccnt += numLinks*2;


	try
	{
		obj = new double[ccnt]; 
		lb = new double[ccnt];
		ub = new double[ccnt];
		ctype = new char[ccnt];
#ifdef VARNAMES
		vnames = new char*[ccnt];
		for(i = 0; i < ccnt; i++)
			vnames[i] = new char[100];
#endif
	}catch(exception& e)
	{
		cout << "Standard exception: " << e.what() << endl;
	}

	//lb, ub, obj for x_ij
	for(i = 0; i < numLinks; i++) 
	{
		lb[i] = 0;	
		ub[i] = 1;
		obj[i] = 0;
		ctype[i] = 'B';
#ifdef VARNAMES
			sprintf(vnames[i],"x_L%d_F%d",i,flowidx);
#endif
	} 
	
	r_offset = i;
	
	//lb, ub, obj for r_ij 
	for(i; i < numLinks*2; i++)
	{
		lb[i] = 0;
		ub[i] = CPX_INFBOUND;
		obj[i] = Flows[flowidx].costs[i-numLinks];
		ctype[i] = 'C';
#ifdef VARNAMES
			sprintf(vnames[i],"r_L%d_F%d",i-r_offset,flowidx);
#endif
	}
	
	rmin_offset = i;
	
	//r_min
	lb[i] = Flows[flowidx].rate;
	ub[i] = CPX_INFBOUND;
	obj[i] = 0;
	ctype[i] = 'C';
#ifdef VARNAMES
			sprintf(vnames[i],"r_min_F%d",flowidx);
#endif
	i++;
	
	s_offset = i;
	
	//s_ij
	for(i; i < numLinks*3+1; i++)
	{
		lb[i] = 0;
		ub[i] = CPX_INFBOUND;
		obj[i] = 0;
		ctype[i] = 'C';
#ifdef VARNAMES
			sprintf(vnames[i],"s_L%d_F%d",i-s_offset,flowidx);
#endif
	}
	
	t_offset = i;
	
	//t
	lb[i] = 0;
	ub[i] = CPX_INFBOUND;
	obj[i] = 0;
	ctype[i] = 'C';
#ifdef VARNAMES
			sprintf(vnames[i],"t_F%d",flowidx);
#endif	
	i++;

if(flagBigM)
{
	c = i;
	rp_offset = c;
	
	//lb, ub, obj for (rp_ij) 
	for(i = 0; i < numLinks; i++)
	{
		lb[c] = 0;
		ub[c] = CPX_INFBOUND;
		obj[c] = 0;
		ctype[c] = 'C';
#ifdef VARNAMES
			sprintf(vnames[c],"rp_L%d_F%d",i,flowidx);
#endif
		c++;
	}
	
	theta_offset = c;
	
	//lb, ub, obj for (theta_ij) 
	for(i = 0; i < numLinks; i++)
	{
		lb[c] = 0;
		ub[c] = CPX_INFBOUND;
		obj[c] = 0;
		ctype[c] = 'C';
#ifdef VARNAMES
			sprintf(vnames[c],"theta_L%d_F%d",i,flowidx);
#endif
		c++;
	}

}	
	//add CPLEX columns
	status = CPXnewcols(env, lp, ccnt, obj, lb, ub, ctype, vnames);
	if(status)
		throw( DCR::DCRException( "DCR_MFSP_SOCP_CPX::addCols(): failed to create CPLEX cols" ) );
		
	delete [] lb;
	delete [] ub;
	delete [] obj;
	delete [] ctype;
	
#ifdef VARNAMES
for(i = 0; i < ccnt; i++)
	delete [] vnames[i];
delete vnames;
#endif
	
}

/*--------------------------------------------------------------------------*/

/* This private method adds flow-conservation constraints
 * \sum_{ji \in BS(i)} x_ji^k  - \sum_{ij \in FS(i)} x_ij^k = b(i)^k \forall i \in N 
 * */
void DCR_MFSP_SOCP_CPX::addFlowRows( int flowidx ) 
{
	
	int r, j, c, rcnt, nzcnt, status;
	int *rmatbeg, *rmatind;
	double *rmatval, *rhs;
	char *sense;

	//add one row at a time
	rcnt = 1;
	try
	{
		rmatbeg = new int[1];
		rhs = new double[1];
		sense = new char[1];
	}catch(exception& e)
	{
		cout << "Standard exception: " << e.what() << endl;
	}
	
	rmatbeg[0] = 0;
	sense[0] = 'E';

	for(r = 0; r < numNodes; r++)
	{
		c = 0;
		nzcnt = getNumInLinks(r) + getNumOutLinks(r);
	
		try
		{
			rmatind = new int[nzcnt];
			rmatval = new double[nzcnt];
		}catch(exception& e)
		{
			cout << "Standard exception: " << e.what() << endl;
		}
		
		for(j = 0; j < numLinks; j++)
		{
			if(isLinkInNode(r, j))
			{
				rmatind[c] = j + FOFFSET*flowidx;
				rmatval[c] = 1;
				c++;
			}	
			else if(isLinkOutNode(r, j))
			{
				rmatind[c] = j + FOFFSET*flowidx;
				rmatval[c] = -1;
				c++;
			}			
		}
			
		if( isFlowSource(r, flowidx) )
			rhs[0] = -1;
		else if ( isFlowSink(r, flowidx) )
			rhs[0] = 1;
		else
			rhs[0] = 0;			
		
		status = CPXaddrows(env, lp, 0, rcnt, nzcnt, rhs, sense, rmatbeg, rmatind, rmatval, NULL, NULL);
		if(status)
		throw( DCR::DCRException( "DCR_SFSP_SOCP_CPX::addFlowRows(): failed to create CPLEX rows" ) );
			
		delete [] rmatval;
		delete [] rmatind;		
	}
	
	delete [] rhs;
	delete [] sense;
	delete [] rmatbeg;
	
}

/*--------------------------------------------------------------------------*/

/* This private method adds rate-bound constraints*/
void DCR_MFSP_SOCP_CPX::addRateBounds( int flowidx )
{

	int r, c, rcnt, nzcnt, status;
	int *rmatbeg, *rmatind;
	double *rmatval, *rhs;
	char *sense;
	
	//r_ij - c_ij x_ij <= 0 \forall (ij) \in A
	rcnt = numLinks;
	nzcnt = 2*rcnt;
	try
	{
		rmatbeg = new int[rcnt];
		rhs = new double[rcnt];
		sense = new char[rcnt];
		rmatind = new int[nzcnt];
		rmatval = new double[nzcnt];
	}catch(exception& e)
	{
		cout << "Standard exception: " << e.what() << endl;
	}
	
	c = 0;
	for(r = 0; r < numLinks; r++)
	{
		rmatbeg[r] = c;
		rhs[r] = 0;
		sense[r] = 'L';
		
		//r_ij
		rmatind[c] = r + numLinks + FOFFSET*flowidx;
		rmatval[c] = 1;
		c++;
		
		//c_ij x_ij
		rmatind[c] = r + FOFFSET*flowidx;
		//rmatval[c] = - Links[r].mutualcap;
		rmatval[c] = - Flows[flowidx].caps[r];
		c++; 
	}
	
	status = CPXaddrows(env, lp, 0, rcnt, nzcnt, rhs, sense, rmatbeg, rmatind, rmatval, NULL, NULL);
	if(status)
		throw( DCR::DCRException( "DCR_SFSP_SOCP_CPX::addRateBouds(): failed to create CPLEX rows" ) );
	
	delete [] rmatbeg;
	delete [] rhs;
	delete [] sense;
	delete [] rmatind;
	delete [] rmatval;
	
	//r_min - r_ij + c_max x_ij <= c_max
	double c_max;
	
	//find c_max for flow flowidx
	c_max = 0;
	for(r = 0; r < numLinks; r++)
	{
		if (Flows[flowidx].caps[r] > c_max)
				c_max = Flows[flowidx].caps[r];
	}

	rcnt = numLinks;
	nzcnt = rcnt*3;
	try
	{
		rmatbeg = new int[rcnt];
		rhs = new double[rcnt];
		sense = new char[rcnt];
		rmatind = new int[nzcnt];
		rmatval = new double[nzcnt];
	}catch(exception& e)
	{
		cout << "Standard exception: " << e.what() << endl;
	}

	c = 0;
	for(r = 0; r < numLinks; r++)
	{
		rmatbeg[r] = c;
		//rhs[r] = Links[r].mutualcap;
		//rhs[r] = Flows[flowidx].caps[r];
		rhs[r] = c_max;
		sense[r] = 'L';
		
		// r_min	
		rmatind[c] = 2*numLinks + FOFFSET*flowidx;
		rmatval[c] = 1;
		c++;
		
		// - r_ij
		rmatind[c] = r + numLinks + FOFFSET*flowidx;
		rmatval[c] = -1;
		c++;
		
		// c_ij x_ij
		rmatind[c] = r + FOFFSET*flowidx;
		//rmatval[c] = Links[r].mutualcap;
		//rmatval[c] = Flows[flowidx].caps[r];
		rmatval[c] = c_max;
		c++;
	}
	
	status = CPXaddrows(env, lp, 0, rcnt, nzcnt, rhs, sense, rmatbeg, rmatind, rmatval, NULL, NULL);
	if(status)
		throw( DCR::DCRException( "DCR_SFSP_SOCP_CPX::addRateBounds(): failed to create CPLEX rows" ) );
	
	delete [] rmatbeg;
	delete [] rhs;
	delete [] sense;
	delete [] rmatind;
	delete [] rmatval;
	
}

/*--------------------------------------------------------------------------*/
/* This private method adds max-delay constraint according
 * to SRP (Strictly-Rate-Proportional) formula 
 * using a SOCP-perspective formulation*/
void DCR_MFSP_SOCP_CPX::addDelayRows( int flowidx )
{
	
	int r, c, rcnt, nzcnt, status;
	int *rmatbeg, *rmatind;
	double *rmatval, *rhs;
	char *sense;
	
	//max delay constraint
	rcnt = 1;
	nzcnt = 1 + 2*numLinks;
	try
	{
		rmatbeg = new int[rcnt];
		rhs = new double[rcnt];
		sense = new char[rcnt];
		rmatind = new int[nzcnt];
		rmatval = new double[nzcnt];
	}catch(exception& e)
	{
		cout << "Standard exception: " << e.what() << endl;
	}

	
	rmatbeg[0] = 0;
	rhs[0] = Flows[flowidx].deadline;
	sense[0] = 'L';
	
	c = 0;
	//t
	rmatind[c] = 3*numLinks + 1 + FOFFSET*flowidx;
	rmatval[c] = 1;
	c++;
	
	//s_ij
	for(r = 0; r < numLinks; r++)
	{
		rmatind[c] = r + numLinks*2 +1 + FOFFSET*flowidx;
		rmatval[c] = 1;  
		c++;
	}
	
	//x_ij
	for(r = 0; r < numLinks; r++)
	{
		rmatind[c] = r + FOFFSET*flowidx;
		rmatval[c] = MTU/Links[r].speed + Links[r].delay + Nodes[Links[r].startnode].delay;  
		c++;
	}
	
	status = CPXaddrows(env, lp, 0, rcnt, nzcnt, rhs, sense, rmatbeg, rmatind, rmatval, NULL, NULL);
	if(status)
		throw( DCR::DCRException( "DCR_SFSP_SOCP_CPX::addDelayRows(): failed to create CPLEX rows" ) );
		
	delete [] rmatbeg;
	delete [] rhs;
	delete [] sense;
	delete [] rmatind;
	delete [] rmatval;
	
	
	//SOCP constraints: s_ij r_ij >= L x_ij^2
	int quadnzcnt, *quadrow, *quadcol;
	double *quadval;
	
    quadnzcnt = 2;
 
	try
	{
		quadrow = new int[quadnzcnt];
		quadcol = new int[quadnzcnt];
		quadval = new double[quadnzcnt];
	}catch(exception& e)
	{
		cout << "Standard exception: " << e.what() << endl;
	}
    
    /*prova 1*/
    /*int linnzcnt;
    int *linind;
    double *linval;
    
    linnzcnt = 1;
    linind = new int[linnzcnt];
    linval = new double[linnzcnt];*/
    /*fine prova 1*/
    
   
    for(r = 0; r < numLinks; r++)
    {
		//s_ij
		quadrow[0] = numLinks*2 +1 + r + FOFFSET*flowidx;
		//r_ij
		quadcol[0] = numLinks + r + FOFFSET*flowidx;
		quadval[0] = 1;
		
		//r_ij
		//quadrow[1] = numLinks + r;
		//s_ij
		//quadcol[1] = numLinks*2 +1 + r;
		//quadval[1] = 0.5;

		
		//- L x_ij^2
		quadrow[1] = r + FOFFSET*flowidx;
		quadcol[1] = r + FOFFSET*flowidx;
		quadval[1] = -MTU;
		
		/*prova 1*/
		// - L x_ij
		/*linind[0] = r;
		linval[0] = - MTU;*/
		//fine prova 1

		status = CPXaddqconstr(env, lp, 0, quadnzcnt, 0, 'G', NULL, NULL, quadrow, quadcol, quadval, NULL);
		if(status)
		throw( DCR::DCRException( "DCR_SFSP_SOCP_CPX::addDelayRows(): failed to create CPLEX quad rows" ) );
	}
	
	
	delete [] quadrow;
	delete [] quadcol;
	delete [] quadval;

	//prova 1
	/*delete [] linind;
	delete [] linval;*/
	//fine prova 1
	
	//SOCP constraints: t r_min >= \sigma	
	quadnzcnt = 1;
 
	try
	{
		quadrow = new int[quadnzcnt];
		quadcol = new int[quadnzcnt];
		quadval = new double[quadnzcnt];
    }catch(exception& e)
	{
		cout << "Standard exception: " << e.what() << endl;
	}
	
    quadrow[0] = numLinks*3 +1 + FOFFSET*flowidx;
	quadcol[0] = numLinks*2 + FOFFSET*flowidx;
	quadval[0] = 1;
		
	status = CPXaddqconstr(env, lp, 0, quadnzcnt, Flows[flowidx].burst, 'G', NULL, NULL, quadrow, quadcol, quadval, NULL);
	if(status)
		throw( DCR::DCRException( "DCR_SFSP_SOCP_CPX::addDelayRows(): failed to create CPLEX quad row" ) );
		
	delete [] quadrow;
	delete [] quadcol;
	delete [] quadval;
}

/*--------------------------------------------------------------------------*/

/* This private method adds max-delay constraint according
 * to SRP (Strictly-Rate-Proportional) formula 
 * using a SOCP-bigM formulation*/
void DCR_MFSP_SOCP_CPX::addDelayRowsBigM( int flowidx )
{
	
	//NB. bigM must be large enough
	//double bigM = 100*Flows[flowidx].deadline; 
	double a, b;
	a = sqrt(MTU);
	b = MTU/Flows[flowidx].rate;
	double bigM = fmax(a, b) + 1;
	
	
	/******************************************************************/
	//max delay constraint:
	//t + sum_{ij} \theta_{ij} + x_{ij}(L/w_{ij} + l_{ij} + n_i) <= delta
	/******************************************************************/ 
	
	int r, c, rcnt, nzcnt, status;
	int *rmatbeg, *rmatind;
	double *rmatval, *rhs;
	char *sense;

	rcnt = 1;
	nzcnt = 1 + 2*numLinks;
	try
	{
		rmatbeg = new int[rcnt];
		rhs = new double[rcnt];
		sense = new char[rcnt];
		rmatind = new int[nzcnt];
		rmatval = new double[nzcnt];
	}catch(exception& e)
	{
		cout << "Standard exception: " << e.what() << endl;
	}
	
	rmatbeg[0] = 0;
	rhs[0] = Flows[flowidx].deadline;
	sense[0] = 'L';
	
	c = 0;
		
	//t
	rmatind[c] = t_offset + FOFFSET*flowidx;
	rmatval[c] = 1;
	c++;
	
	//theta_ij 		
	for(r = 0; r < numLinks; r++)
	{
		rmatind[c] = theta_offset + r + FOFFSET*flowidx;
		rmatval[c] = 1;  
		c++;
	}
	
	//x_ij
	for(r = 0; r < numLinks; r++)
	{
		rmatind[c] = r + FOFFSET*flowidx;
		rmatval[c] = MTU/Links[r].speed + Links[r].delay + Nodes[Links[r].startnode].delay;  
		c++;
	}
	
	status = CPXaddrows(env, lp, 0, rcnt, nzcnt, rhs, sense, rmatbeg, rmatind, rmatval, NULL, NULL);
	if(status)
		throw( DCR::DCRException( "DCR_MFSP_SOCP_CPX::addDelayRows(): failed to create CPLEX rows" ) );
		
	delete [] rmatbeg;
	delete [] rhs;
	delete [] sense;
	delete [] rmatind;
	delete [] rmatval;		
	
	/*******************************************************************/
	//SOCP constraints: s_ij rp_ij >= L 
	/*******************************************************************/
	
	int quadnzcnt, *quadrow, *quadcol;
	double *quadval;
 
    quadnzcnt = 1;
 
	try
	{
		quadrow = new int[quadnzcnt];
		quadcol = new int[quadnzcnt];
		quadval = new double[quadnzcnt];
	}catch(exception& e)
	{
		cout << "Standard exception: " << e.what() << endl;
	}
    
 	for(r = 0; r < numLinks; r++)
	{
		//s_ij
		quadrow[0] = s_offset + r + FOFFSET*flowidx;
		//rp_ij
		quadcol[0] = rp_offset + r + FOFFSET*flowidx;
		quadval[0] = 1;		
			
		status = CPXaddqconstr(env, lp, 0, quadnzcnt, MTU, 'G', NULL, NULL, quadrow, quadcol, quadval, NULL);
		if(status)
		throw( DCR::DCRException( "DCR_MFSP_SOCP_CPX::addDelayRows(): failed to create CPLEX quad rows" ) );
	
	}
	
	delete [] quadrow;
	delete [] quadcol;
	delete [] quadval;
	
	/*******************************************************************/
	//SOCP constraints: t r_min >= \sigma	
	/*******************************************************************/
	
	quadnzcnt = 1;
 
	try
	{
		quadrow = new int[quadnzcnt];
		quadcol = new int[quadnzcnt];
		quadval = new double[quadnzcnt];
    }catch(exception& e)
	{
		cout << "Standard exception: " << e.what() << endl;
	}
	
	//t
	quadrow[0] = t_offset + FOFFSET*flowidx;
	
	//r_min
	quadcol[0] = rmin_offset + FOFFSET*flowidx;
	quadval[0] = 1;
		
	status = CPXaddqconstr(env, lp, 0, quadnzcnt, Flows[flowidx].burst, 'G', NULL, NULL, quadrow, quadcol, quadval, NULL);
	if(status)
		throw( DCR::DCRException( "DCR_MFSP_SOCP_CPX::addDelayRows(): failed to create CPLEX quad row" ) );
	
	delete [] quadrow;
	delete [] quadcol;
	delete [] quadval;
	
	/******************************************************************/
	//theta_ij - M x_ij <= 0
	/******************************************************************/

	rcnt = numLinks;
	nzcnt = 2*rcnt;
		
	try
	{
		rmatbeg = new int[rcnt];
		rhs = new double[rcnt];
		sense = new char[rcnt];
		rmatind = new int[nzcnt];
		rmatval = new double[nzcnt];
	}catch(exception& e)
	{
		cout << "Standard exception: " << e.what() << endl;
	}
	
	c = 0;
	for(r = 0; r < numLinks; r++)
	{
		rmatbeg[r] = c;
		rhs[r] = 0;
		sense[r] = 'L';
		
		//theta_ij
		rmatind[c] = theta_offset + r + FOFFSET*flowidx;
		rmatval[c] = 1;
		c++;
		
		//- M x_ij
		rmatind[c] = r + FOFFSET*flowidx;
		rmatval[c] = - bigM; //M is set equal to the flow deadline
		c++; 
	}
	
	status = CPXaddrows(env, lp, 0, rcnt, nzcnt, rhs, sense, rmatbeg, rmatind, rmatval, NULL, NULL);
	if(status)
		throw( DCR::DCRException( "DCR_MFSP_SOCP_CPX::addDelayRows(): failed to create CPLEX rows" ) );
	
	delete [] rmatbeg;
	delete [] rhs;
	delete [] sense;
	delete [] rmatind;
	delete [] rmatval;
	
	/******************************************************************/
	//theta_ij - s_ij + - M x_ij >= - M
	/******************************************************************/
	
	rcnt = numLinks;
	nzcnt = 3*rcnt;
		
	try
	{
		rmatbeg = new int[rcnt];
		rhs = new double[rcnt];
		sense = new char[rcnt];
		rmatind = new int[nzcnt];
		rmatval = new double[nzcnt];
	}catch(exception& e)
	{
		cout << "Standard exception: " << e.what() << endl;
	}
	
	c = 0;
	for(r = 0; r < numLinks; r++)
	{
		rmatbeg[r] = c;
		rhs[r] = - bigM;
		sense[r] = 'G';
		
		//theta_ij
		rmatind[c] = theta_offset + r + FOFFSET*flowidx;
		rmatval[c] = 1;
		c++;
		
		//- M x_ij
		rmatind[c] = r + FOFFSET*flowidx;
		rmatval[c] = - bigM; //M is set equal to the flow deadline
		c++; 
			
		//- s_ij
		rmatind[c] = s_offset + r + FOFFSET*flowidx;
		rmatval[c] = -1;
		c++;
	}
	
	status = CPXaddrows(env, lp, 0, rcnt, nzcnt, rhs, sense, rmatbeg, rmatind, rmatval, NULL, NULL);
	if(status)
		throw( DCR::DCRException( "DCR_MFSP_SOCP_CPX::addDelayRows(): failed to create CPLEX rows" ) );
	
	delete [] rmatbeg;
	delete [] rhs;
	delete [] sense;
	delete [] rmatind;
	delete [] rmatval;
	
	/******************************************************************/
	//rp_ij^k - r_ij^k + M x_ij^k  <= M
	/******************************************************************/
	
	rcnt = numLinks;
	nzcnt = 3*rcnt;
		
	try
	{
		rmatbeg = new int[rcnt];
		rhs = new double[rcnt];
		sense = new char[rcnt];
		rmatind = new int[nzcnt];
		rmatval = new double[nzcnt];
	}catch(exception& e)
	{
		cout << "Standard exception: " << e.what() << endl;
	}
	
	c = 0;
	for(r = 0; r < numLinks; r++)
	{
		rmatbeg[r] = c;
		rhs[r] = bigM;
		sense[r] = 'L';
		
		//rp_ij
		rmatind[c] = rp_offset + r + FOFFSET*flowidx;
		rmatval[c] = 1;
		c++;
		
		// M x_ij
		rmatind[c] = r + FOFFSET*flowidx;
		rmatval[c] =  bigM; //M is set equal to the flow deadline
		c++; 
			
		//- r_ij
		rmatind[c] = r_offset + r + FOFFSET*flowidx;
		rmatval[c] = -1;
		c++;
	}
	
	status = CPXaddrows(env, lp, 0, rcnt, nzcnt, rhs, sense, rmatbeg, rmatind, rmatval, NULL, NULL);
	if(status)
		throw( DCR::DCRException( "DCR_MFSP_SOCP_CPX::addDelayRows(): failed to create CPLEX rows" ) );
	
	delete [] rmatbeg;
	delete [] rhs;
	delete [] sense;
	delete [] rmatind;
	delete [] rmatval;

	
}
/*--------------------------------------------------------------------------*/

/* This private method adds mutual capacity constraints in the SOCP model*/
/* \sum_{k \in FLows} r_ij^k <= c_ij \forall (i,j) \in Links*/
void DCR_MFSP_SOCP_CPX::addCapacityRows( void )
{	
	
	int r, c, rcnt, nzcnt, status, k;
	int *rmatbeg, *rmatind;
	double *rmatval, *rhs;
	char *sense;
	
	rcnt = numLinks;
	nzcnt = rcnt * numFlows;
	
	rmatbeg = new int[rcnt];
	rhs = new double[rcnt];
	sense = new char[rcnt];
	rmatind = new int[nzcnt];
	rmatval = new double[nzcnt];
	
	c = 0;
	
	for (r = 0; r < numLinks; r++)
	{
		rmatbeg[r] = c;
		rhs[r] = Links[r].capacity;
		sense[r] = 'L'; 
		
		for(k = 0; k < numFlows; k++)
		{
			rmatind[c] = FOFFSET*k + r_offset + r;
			rmatval[c] = 1;
			c++;
		}
	}
	
	status = CPXaddrows(env, lp, 0, rcnt, nzcnt, rhs, sense, rmatbeg, rmatind, rmatval, NULL, NULL);
	if(status)
		throw( DCR::DCRException( "DCR_MFSP_SOCP_CPX::addCapacityRows(): failed to create CPLEX rows" ) );
		
	delete [] rmatbeg;
	delete [] rhs;
	delete [] sense;
	delete [] rmatind;
	delete [] rmatval;	
	
}

/*--------------------------------------------------------------------------*/

/* This private method returns the number of entering arcs into node*/
int DCR_MFSP_SOCP_CPX::getNumInLinks(int node)
{
	int j, nin;
	
	nin = 0;
	for(j = 0; j < numLinks; j++)
	{
		if (Links[j].endnode == node) nin++;
	}
	
	return nin;
}

/*--------------------------------------------------------------------------*/

/* This private method returns the number of exiting arcs from node*/
int DCR_MFSP_SOCP_CPX::getNumOutLinks(int node)
{
	int j, nout;
	
	nout = 0;
	for(j = 0; j < numLinks; j++)
	{
		if(Links[j].startnode == node) nout++;
	}
	
	return nout;
}

/*--------------------------------------------------------------------------*/

/* This private method checks if arc is entering node*/
bool DCR_MFSP_SOCP_CPX::isLinkInNode(int node, int arc)
{
	if(Links[arc].endnode == node) return true;
	
	return false;
}

/*--------------------------------------------------------------------------*/

/* This private method checks if arc is exiting from node*/
bool DCR_MFSP_SOCP_CPX::isLinkOutNode(int node, int arc)
{
	if(Links[arc].startnode == node)return true;
	
	return false;
}

/*--------------------------------------------------------------------------*/

/* This private method checks if node is the source of flow*/
bool DCR_MFSP_SOCP_CPX::isFlowSource(int node, int flow)
{
	if(Flows[flow].sourcenode == node) return true;
	
	return false;
}

/*--------------------------------------------------------------------------*/

/* This private method checks if node is the sink of flow*/
bool DCR_MFSP_SOCP_CPX::isFlowSink(int node, int flow)
{
	if(Flows[flow].sinknode == node ) return true;
	
	return false;
}

/*--------------------------------------------------------------------------*/

/* This private method is used to copy array data (flows, links and nodes) from the load() function 
 * into the class as class members.*/
void DCR_MFSP_SOCP_CPX::copyDataArrays(DCRFlow *flows, DCRLink *links, DCRNode *nodes)
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

/* This private method is used to get the objval and the solution vector when
 * DCRsolve() function returns OK.
 * It can be either a continuous solution or a mixed-integer solution.
 * In both cases it stores the solution splitting it into
 * Xsol array (for path variables) and Rsol array (for rate variables).
 * */
void DCR_MFSP_SOCP_CPX::getSol()
{
	int status;
	int i, k, c, nc;
	double *sol;
	
	/////////////////////////////////////////
	//GET OBJVAL
	//////////////////////////////////////////
	status = CPXgetobjval (env, lp, &objval);
	if( status )
		throw( DCR::DCRException( "DCR_MFSP_SOCP_CPX::getSol(): failed to get CPLEX objval" ) );	
	
	///////////////////////////////////////
	//GET SOLUTION
	////////////////////////////////////////
	
	nc = CPXgetnumcols(env, lp);
	
	try
	{
		sol = new double[nc];
	}catch(exception& e)
	{
		cout << "Standard exception: " << e.what() << endl;
	}
	
	//get the whole solution vector from cplex
	status = CPXgetx (env, lp, sol, 0, nc-1);
	if( status )
		throw( DCR::DCRException( "DCR_MFSP_SOCP_CPX::getSol(): failed to get CPLEX solution" ) );
	
	c = 0;
	
	//store path variables in XSol
	for(k = 0; k < numFlows; k++)
	{
		for(i = 0; i < numLinks; i++)
		{
			Xsol[c] = sol[k*FOFFSET + i];
			c++;
		}
	}
	
	c = 0;

	//store rate variables in RSol	
	for(k = 0; k < numFlows; k++)
	{
		//cout << "\n rmin of flow " << k << " : " << sol[k*FOFFSET + rmin_offset] << endl; 
		for(i = 0; i < numLinks; i++)
		{
			Rsol[c] = sol[k*FOFFSET + r_offset + i];
			c++;
		}
	}
		
	delete [] sol;
	
}

/*--------------------------------------------------------------------------*/

/* This private method deallocates dinamically allocated method */

void DCR_MFSP_SOCP_CPX::clean_up()
{
	int status = 0;
		
	if ( lp != 0 ) 
	{
		status = CPXfreeprob (env, &lp);
		if ( status ) 
			throw( DCR::DCRException( "DCR_MFSP_SOCP_CPX::destructor(): failed to free CPLEX prob." ) );
	}
	
	if ( env != 0 ) 
	{
		status = CPXcloseCPLEX (&env);
		if ( status ) 
			throw( DCR::DCRException( "DCR_MFSP_SOCP_CPX::destructor(): failed to free CPLEX env." ) );
	}
		
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
/*-------------------------- End File DCR_MFSP_SOCP_CPX.cpp-----------------*/
/*--------------------------------------------------------------------------*/
