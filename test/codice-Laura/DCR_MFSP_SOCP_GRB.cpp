/*--------------------------------------------------------------------------*/
/*---------------------------- File DCR_MFSP_SOCP_GRB.cpp ----------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------------- IMPLEMENTATION -------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

//#define DEBUG_GRB
//#define VARNAMES

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "DCR_MFSP_SOCP_GRB.h"
#include <fstream>
#include <cstdlib>
#include <math.h>

/*--------------------------------------------------------------------------*/
/*----------------------- IMPLEMENTATION OF DCR_MFSP_SOCP_GRB---------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------- PUBLIC METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

/*----------------------------- CONSTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

DCR_MFSP_SOCP_GRB::DCR_MFSP_SOCP_GRB() : DCR()
{
	//init all class members
	Nodes = 0;
	Links = 0;
	Flows = 0;
		
	env = 0;
	model = 0;
	
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

void DCR_MFSP_SOCP_GRB::DCRwrite( char *filename )
{	
	int error;
	
	error = GRBwrite(model, filename);
	if (error) 
	{
		printf("%s\n", GRBgeterrormsg(env));
		throw( DCR::DCRException( "DCR_MFSP_SOCP_GRB::DCRwrite(): failed to write GRB problem" ) );
	}
}

/*--------------------------------------------------------------------------*/
/*----------------------------------SOLVE-----------------------------------*/
/*--------------------------------------------------------------------------*/

DCR::DCRStatus DCR_MFSP_SOCP_GRB::DCRsolve( void )
{
	int error;
	int solstatus;
	
	//debug
	//cout << "\nSOCP GRB solver " << endl;
	
	error = GRBoptimize(model);
	if (error) 
	{
		printf("%s\n", GRBgeterrormsg(env));
		throw( DCR::DCRException( "DCR_MFSP_SOCP_GRB::DCRsolve(): failed to solve GRB problem" ) );
	}	
	
	error = GRBgetintattr(model, "Status", &solstatus);
	if (error)
	{		printf("%s\n", GRBgeterrormsg(env));
			throw( DCR::DCRException( "DCR_MFSP_SOCP_GRB::DCRsolve: failed to get GUROBI status" ) );
	}
	
	if(solstatus == 2) 
	{
		int ismip;
		
		error = GRBgetintattr(model, "IsMIP", &ismip);
		if (error)
		{	printf("%s\n", GRBgeterrormsg(env));
			throw( DCR::DCRException( "DCR_MFSP_SOCP_GRB::DCRsolve: failed to get GUROBI model type" ) );
		}
		
		if(ismip) UBsol = true; 
		else      Psol = true; 	
		getSol();				
		return OK;
	}
	
	UBsol = false;
	
	if(solstatus == 3) return Infeasible;
	
	//FIXME: add more solver status codes...
	
	return Error;
	
}

/*--------------------------------------------------------------------------*/
/*----------------------------------READ RESULTS-----------------------------*/
/*--------------------------------------------------------------------------*/


double DCR_MFSP_SOCP_GRB::DCRgetObj( void )
{	
	return objval;
}


/*--------------------------------------------------------------------------*/

double DCR_MFSP_SOCP_GRB::DCRgetNodes( void )
{
	int error;
	double nodecnt;
		
	error = GRBgetdblattr(model, "NodeCount", &(nodecnt));
	if (error) 
	{
		printf("%s\n", GRBgeterrormsg(env));
		throw( DCR::DCRException( "DCR_MFSP_SOCP_GRB::DCRgetNodes(): failed to get GRB nodes" ) );
	}
	
	return nodecnt;
}

/*--------------------------------------------------------------------------*/

double DCR_MFSP_SOCP_GRB::DCRgetRuntime( void )
{
	int error;
	double runtime;
		
	error = GRBgetdblattr(model, "Runtime", &(runtime));
	if (error) 
	{
		printf("%s\n", GRBgeterrormsg(env));
		throw( DCR::DCRException( "DCR_MFSP_SOCP_GRB::DCRgetRuntime(): failed to get GRB runtime" ) );
	}
	
	return runtime;
}

/*--------------------------------------------------------------------------*/

int DCR_MFSP_SOCP_GRB::DCRgetUBSolNPaths( int k )
{
	if(k < 0 || k>= numFlows)
		throw( DCR::DCRException( "DCR_MFSP_SOCP_GRB::DCRgetUBSolNPaths(): wrong flow index." ) );
	if(UBsol) return 1;
	else
		throw( DCR::DCRException( "DCR_MFSP_SOCP_GRB::DCRgetUBSolNPaths(): no solution found." ) );	//FIXME: is this okay???
}

/*--------------------------------------------------------------------------*/

int DCR_MFSP_SOCP_GRB::DCRgetUBSolPath( int k, int p, int *X, double *R)
{
	int i, nhops;
	
	
	if(k < 0 || k>= numFlows)
		throw( DCR::DCRException( "DCR_MFSP_SOCP_GRB::DCRgetUBSolPath(): wrong flow index." ) );
	if(p != 0)
		throw( DCR::DCRException( "DCR_MFSP_SOCP_GRB::DCRgetPUBolPath(): wrong path index." ) );
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
		throw( DCR::DCRException( "DCR_MFSP_SOCP_GRB::DCRgetUBSolPath(): no solution found." ) );//FIXME: is this okay ???		
	
	return nhops;
}

/*--------------------------------------------------------------------------*/

int DCR_MFSP_SOCP_GRB::DCRgetPSolNPaths( int k )
{
	if(k < 0 || k>= numFlows)
		throw( DCR::DCRException( "DCR_MFSP_SOCP_GRB::DCRgetPSolNPaths(): wrong flow index." ) );
	if(Psol) return 1;
	else
		throw( DCR::DCRException( "DCR_MFSP_SOCP_GRB::DCRgetPSolNPaths(): no solution found." ) );	//FIXME: is this okay???
}

/*--------------------------------------------------------------------------*/

void  DCR_MFSP_SOCP_GRB::DCRgetPSolPath(int k, int p, double *X, double *R)
{
	int i;
	
	if(k < 0 || k>= numFlows)
		throw( DCR::DCRException( "DCR_MFSP_SOCP_GRB::DCRgetPSolPath(): wrong flow index." ) );
	if(p != 0)
		throw( DCR::DCRException( "DCR_MFSP_SOCP_GRB::DCRgetPSolPath(): wrong path index." ) );
	if(Psol)
	{
		for(i = 0; i < numLinks; i++)
		{
			X[i] = Xsol[i + k*numLinks];
			R[i] = Rsol[i + k*numLinks];
		}		
	}
	else
		throw( DCR::DCRException( "DCR_MFSP_SOCP_CPX::DCRgetPSolPath(): no solution found." ) );	//FIXME: is this okay ???
}


/*--------------------------------------------------------------------------*/
/*----------------------------------LOAD DATA ------------------------------*/
/*--------------------------------------------------------------------------*/

void DCR_MFSP_SOCP_GRB::DCRloadProblem(int nnodes, int nlinks, int nflows, 
	DCRFlow *flows, DCRLink *links, DCRNode *nodes, double mtu, DCRDelay deltype)
{
	int k;
	int error;
	
	if(deltype != SRP)
		throw( DCR::DCRException( "DCR_MFSP_SOCP_CPX::DCRloadProblem(): this delay formula is not supported yet..sorry" ) );
	
	///////////////////////////////////
	//INIT
	////////////////////////////////// 
	clean_up();	
	
	Psol = false;
	UBsol = false;
	objval = Inf<double>();	
	///////////////////////////////////////////////////
	
	numNodes = nnodes;
	numLinks = nlinks;
	numFlows = nflows;	
	MTU = mtu;
	dtype = deltype;	
	copyDataArrays(flows, links, nodes);
	Xsol = new double[numLinks*numFlows];
	Rsol = new double[numLinks*numFlows];		
	
	//open GUROBI env
	 error = GRBloadenv(&env, "dcr.log");
	if (error) 
	{
		printf("%s\n", GRBgeterrormsg(env));
		throw( DCR::DCRException( "DCR_MFSP_SOCP_GRB::DCRloadProblem(): failed to create GUROBI env" ) );
	} 		
			
	//create GUROBI problem
	error = GRBnewmodel(env, &model, "dcr", 0, NULL, NULL, NULL, NULL, NULL);
	if (error) 
	{
		printf("%s\n", GRBgeterrormsg(env));
		throw( DCR::DCRException( "DCR_MFSP_SOCP_GRB::DCRloadProblem(): failed to create GUROBI model" ) );
	} 			

	//set GUROBI obj sen (1 = minimize)
	error = GRBsetintattr(model, "ModelSense", 1);
	if (error) 
	{
		printf("%s\n", GRBgeterrormsg(env));
		throw( DCR::DCRException( "DCR_MFSP_SOCP_GRB::DCRloadProblem(): failed to set GUROBI obj sense" ) );
	} 	
	
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

void DCR_MFSP_SOCP_GRB::DCRsetModel(char method)
{
	if(method=='b') flagBigM = 1;	
	else if(method=='p') flagBigM = 0;
	else throw( DCR::DCRException( "DCR_MFSP_SOCP_GRB::DCRsetModel(): failed to set model" ) );
}

/*--------------------------------------------------------------------------*/

void DCR_MFSP_SOCP_GRB::DCRsetOptEps( double OE )
{
	int error;
	
	optEps = OE;
		
	//FIXME: is this okay -- see also CPLEX version ???
	error = GRBsetdblparam (env, "IntFeasTol", optEps);
	if ( error ) 
		throw( DCR::DCRException( "DCR_MFSP_SOCP_GRB::DCRsetOptEps(): failed to set GUROBI integer feasibility tolerance" ) );
	
	
	error = GRBsetdblparam(env,"OptimalityTol",optEps);           //simplex tolerance for optimality
	if ( error ) 
		throw( DCR::DCRException( "DCR_MFSP_SOCP_GRB::DCRsetOptEps(): failed to set GUROBI optimality tolerance" ) );
}

/*--------------------------------------------------------------------------*/

void DCR_MFSP_SOCP_GRB::DCRsetFsbEps( double OE )
{
	int error;
	
	fsbEps = OE;
		
	error = GRBsetdblparam(env,"FeasibilityTol",fsbEps);           //simplex tolerance for feasibility
	if ( error ) 
		throw( DCR::DCRException( "DCR_MFSP_SOCP_GRB::DCRsetFsbEps(): failed to set GUROBI feasibility tolerance" ) );
}

/*--------------------------------------------------------------------------*/

void DCR_MFSP_SOCP_GRB::DCRsetQCPType( void )
{
	int i, error, k;
	char *ctype;
	
	try{
		ctype = new char[numLinks];
	}catch(exception& e)
	{
		cout << "Standard exception: " << e.what() << endl;
	}
	
	for(i = 0; i < numLinks; i++)
	{
		ctype[i] = 'C';
	}

	for(k = 0; k < numFlows; k++)
	{		
		error = GRBsetcharattrarray (model, "VType", k*FOFFSET, numLinks, ctype);
		if( error )
				throw( DCR::DCRException( "DCR_MFSP_SOCP_GRB::DCRsetQCPType(): failed to set GUROBI vars type" ) );
	}
	
	delete [] ctype;		
}
	
/*--------------------------------------------------------------------------*/

void DCR_MFSP_SOCP_GRB::DCRsetMIQCPType( void )
{
int i, error, k;
	char *ctype;
	
	try{
		ctype = new char[numLinks];
	}catch(exception& e)
	{
		cout << "Standard exception: " << e.what() << endl;
	}
	
	for(i = 0; i < numLinks; i++)
	{
		ctype[i] = 'B';
	}
	
	for(k = 0; k < numFlows; k++)
	{
		error = GRBsetcharattrarray (model, "VType", k*FOFFSET , numLinks, ctype);
		if( error )
			throw( DCR::DCRException( "DCR_MFSP_SOCP_GRB::DCRsetQCPType(): failed to set GUROBI vars type" ) );
	}
	
	delete [] ctype;		
}

/*--------------------------------------------------------------------------*/

void DCR_MFSP_SOCP_GRB::DCRsetCuts( int setval )
{
	int error;
	//setval = 0 -> switch off grb cuts:
	
	error = GRBsetintparam (env, "Cuts", setval);
	if ( error ) 
		throw( DCR::DCRException( "DCR_MFSP_SOCP_GRB::DCRsetCuts(): failed to set GUROBI cuts" ) );
}

/*--------------------------------------------------------------------------*/

void DCR_MFSP_SOCP_GRB::DCRsetNodeLimit( double nodelimit )
{
	int error;
	
	error = GRBsetdblparam (env, "NodeLimit", nodelimit);
	if ( error ) 
		throw( DCR::DCRException( "DCR_MFSP_SOCP_GRB::DCRsetNodeLimit(): failed to set GUROBI  node limit" ) );
}

/*--------------------------------------------------------------------------*/

void DCR_MFSP_SOCP_GRB::DCRsetHeur( double setval )
{
	int error;
	
	error = GRBsetdblparam(env, "Heuristics", setval);
	if ( error ) 
		throw( DCR::DCRException( "DCR_MFSP_SOCP_GRB::DCRsetHeur(): failed to set GUROBI heurs " ) );
}

/*--------------------------------------------------------------------------*/

void DCR_MFSP_SOCP_GRB::DCRsetNodeFile(double setval)
{
	int error;
	
	error = GRBsetdblparam(env, "NodefileStart", setval);
	if ( error ) 
		throw( DCR::DCRException( "DCR_MFSP_SOCP_GRB::DCRsetNoHeur(): failed to set GUROBI node file " ) );
}

/*--------------------------------------------------------------------------*/

void DCR_MFSP_SOCP_GRB::DCRsetLog(char *filename)
{
	int error;	
		
	error = GRBsetstrparam(env,"LogFile",filename);
	if( error )
			throw( DCR::DCRException( "DCR_MFSP_SOCP_GRB::DCRsetLog(): failed to set GUROBI log" ) );

}

/*--------------------------------------------------------------------------*/

void DCR_MFSP_SOCP_GRB::DCRsetMIQCPStrat(int setval)
{
	int error;
	
	error = GRBsetintparam(env, "MIQCPMethod", setval);
	if( error )
			throw( DCR::DCRException( "DCR_MFSP_SOCP_GRB::DCRsetMIQCPStrat(): failed to set GUROBI MIQCP strategy" ) );
}

/*--------------------------------------------------------------------------*/

void DCR_MFSP_SOCP_GRB::DCRsetTimeLimit( double setval )
{
	int error;
	
	//time limit in seconds 
	error = GRBsetdblparam(env, "TimeLimit", setval);
	if( error )
			throw( DCR::DCRException( "DCR_MFSP_SOCP_GRB::DCRsetTimeLimit(): failed to set GUROBI time limit" ) );	
}

/*--------------------------------------------------------------------------*/
/*----------------------------------MODIFIERS-------------------------------*/
/*--------------------------------------------------------------------------*/

void DCR_MFSP_SOCP_GRB::DCRcloseArcs( int * whch, int na )
{
	int error;
	int j, k;
	
	// close all arcs with indices in whch for all flows
	// this is done using GUROBI function GRBsetdblattrelement
	// and setting the upper bound of the corresponding
	// variables to 0. 

	for(k = 0; k < numFlows; k++)
	{
		for(j = 0; j < na; j++)
		{
				error = GRBsetdblattrelement(model, "UB", whch[j] + FOFFSET*k, 0.0);
				if( error )
				throw( DCR::DCRException( "DCR_MFSP_SOCP_GRB::DCRsetdblattrelement(): failed to set GUROBI var UB" ) );	
		}		
	}
		
}

void DCR_MFSP_SOCP_GRB::DCRcloseArcs( int k , int * whch, int na )
{
	int error;
	int j;
	
	// close all arcs with indices in whch for flow k
	// this is done using GUROBI function GRBsetdblattrelement
	// and setting the upper bound of the corresponding
	// variables to 0. 


	for(j = 0; j < na; j++)
	{
		
		error = GRBsetdblattrelement(model, "UB", whch[j] + FOFFSET*k, 0.0);
		if( error )
			throw( DCR::DCRException( "DCR_MFSP_SOCP_GRB::DCRsetdblattrelement(): failed to set GUROBI var UB" ) );	
	}		
	
}

/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

DCR_MFSP_SOCP_GRB::~DCR_MFSP_SOCP_GRB()
{

	clean_up();	
	
	
}

/*--------------------------------------------------------------------------*/
/*------------------------------ LOCAL (PRIVATE) METHODS -------------------*/
/*--------------------------------------------------------------------------*/


void DCR_MFSP_SOCP_GRB::addCols( int flowidx )
/* Number of vars is: 
 * *****************************
 * In case of perspective SOCP model:
 * *****************************
 * #x_ij + #r_ij + r_min
 * If m is the number of links, we have 2*m vars
 * The order of the vars is:
 * x_link1, ... ,x_linkm, r_link1, ... ,r_linkm, r_min, s_link1, ..., s_linkm, t
 * */
{
	int i, error, ncols, c;
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
		ctype[i] = GRB_BINARY;
#ifdef VARNAMES
			sprintf(vnames[i],"x_L%d_F%d",i,flowidx);
#endif
	}
	
	r_offset = i; 
	
	//lb, ub, obj for r_ij 
	for(i; i < numLinks*2; i++)
	{
		lb[i] = 0;
		ub[i] = GRB_INFINITY;
		obj[i] = Flows[flowidx].costs[i-numLinks];
		ctype[i] = GRB_CONTINUOUS;
#ifdef VARNAMES
			sprintf(vnames[i],"r_L%d_F%d",i-r_offset,flowidx);
#endif
	}
	
	rmin_offset = i;
	
	//r_min
	lb[i] = Flows[flowidx].rate;
	ub[i] = GRB_INFINITY;
	obj[i] = 0;
	ctype[i] = GRB_CONTINUOUS;
#ifdef VARNAMES
			sprintf(vnames[i],"r_min_F%d",flowidx);
#endif
	i++;
	
	s_offset = i;
	
	//s_ij
	for(i; i < numLinks*3+1; i++)
	{
		lb[i] = 0;
		ub[i] = GRB_INFINITY;
		obj[i] = 0;
		ctype[i] = GRB_CONTINUOUS;
#ifdef VARNAMES
			sprintf(vnames[i],"s_L%d_F%d",i-s_offset,flowidx);
#endif
	}
	
	t_offset = i;
	
	//t
	lb[i] = 0;
	ub[i] = GRB_INFINITY;
	obj[i] = 0;
	ctype[i] = GRB_CONTINUOUS;
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
		ub[c] = GRB_INFINITY;
		obj[c] = 0;
		ctype[c] = GRB_CONTINUOUS;
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
		ub[c] = GRB_INFINITY;
		obj[c] = 0;
		ctype[c] = GRB_CONTINUOUS;
#ifdef VARNAMES
			sprintf(vnames[c],"theta_L%d_F%d",i,flowidx);
#endif
		c++;
	}
	
}
	
	//add GRB columns
	error = GRBaddvars(model, ccnt, 0, NULL, NULL, NULL, obj, lb, ub, ctype, NULL);
	if( error ) 
	{ 	printf("%s\n", GRBgeterrormsg(env));
		throw( DCR::DCRException( "DCR_SFSP_SOCP_GRB::addCols(): failed to create GUROBI cols" ) );
	}
		
	error = GRBupdatemodel(model);
	if (error)
	{		printf("%s\n", GRBgeterrormsg(env));
			throw( DCR::DCRException( "DCR_SFSP_SOCP_GRB::addCols(): failed to update GUROBI cols" ) );
	}
	
	#ifdef DEBUG_GRB
	error = GRBgetintattr(model, "NumVars", &ncols);
	
	if( !error )
		cout << "\n***OK: Added " << ncols << " GRB cols ***\n"; 
	#endif
	
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
void DCR_MFSP_SOCP_GRB::addFlowRows( int flowidx )
{

	int r, j, c, rcnt, nzcnt, error, nrows;
	int *rmatind;
	double *rmatval, rhs;
	char sense;

	//add one row at a time
	sense = GRB_EQUAL;
	for(r = 0; r < numNodes; r++)
	{
		c = 0;
		nzcnt = getNumInLinks(r) + getNumOutLinks(r);
	
	try{
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
			rhs = -1;
		else if ( isFlowSink(r, flowidx) )
			rhs = 1;
		else
			rhs = 0;			
		
		error = GRBaddconstr(model, nzcnt, rmatind, rmatval, sense, rhs, NULL);
		if (error) 
		{
			printf("%s\n", GRBgeterrormsg(env));	
			throw( DCR::DCRException( "DCR_SFSP_SOCP_GRB::addFlowRows(): failed to create GUROBI rows" ) );
		}
		
		error = GRBupdatemodel(model);
		if (error)
		{
			printf("%s\n", GRBgeterrormsg(env));
			throw( DCR::DCRException( "DCR_SFSP_SOCP_GRB::addFlowRows(): failed to update GUROBI rows" ) );
		}
		
			
		delete [] rmatval;
		delete [] rmatind;		
	}
	
	
	
	#ifdef DEBUG_GRB
	error = GRBgetintattr(model, "NumConstrs", &nrows);
	
	if( !error )
		cout << "\n***OK: Added " << nrows << " GRB rows (+ flow rows)***\n"; 
	#endif
	
}

/*--------------------------------------------------------------------------*/

/* This private method adds rate-bound constraints*/
void DCR_MFSP_SOCP_GRB::addRateBounds( int flowidx )
{
	int r, c, rcnt, nzcnt, error, nrows;
	int *rmatbeg, *rmatind;
	double *rmatval, *rhs;
	char *sense;
	
	//r_ij - c_ij x_ij <= 0 \forall (ij) \in A
	rcnt = numLinks;
	nzcnt = 2*rcnt;
	
	try{
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
		sense[r] = GRB_LESS_EQUAL;
		
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
	
		
	error = GRBaddconstrs(model, rcnt, nzcnt, rmatbeg, rmatind, rmatval, sense, rhs, NULL);
	if (error)
	{
		printf("%s\n", GRBgeterrormsg(env));
		throw( DCR::DCRException( "DCR_SFSP_SOCP_GRB::addRateBounds(): failed to create GUROBI rows" ) );
	}
	
	error = GRBupdatemodel(model);
	if (error)
	{
		printf("%s\n", GRBgeterrormsg(env));
		throw( DCR::DCRException( "DCR_SFSP_SOCP_GRB::addRateBounds(): failed to update GUROBI rows" ) );
	}
	
	
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
	
	try{
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
		sense[r] = GRB_LESS_EQUAL;
		
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
	
	error = GRBaddconstrs(model, rcnt, nzcnt, rmatbeg, rmatind, rmatval, sense, rhs, NULL);
	if (error)
	{
		printf("%s\n", GRBgeterrormsg(env));
		throw( DCR::DCRException( "DCR_SFSP_SOCP_GRB::addRateBounds() cont.: failed to create GUROBI rows" ) );
	}
	
	error = GRBupdatemodel(model);
	if (error)
	{
		printf("%s\n", GRBgeterrormsg(env));
		throw( DCR::DCRException( "DCR_SFSP_SOCP_GRB::addRateBounds(): failed to update GUROBI rows" ) );
	}
		
	delete [] rmatbeg;
	delete [] rhs;
	delete [] sense;
	delete [] rmatind;
	delete [] rmatval;	
	

	
	#ifdef DEBUG_GRB
	error = GRBgetintattr(model, "NumConstrs", &nrows);
	
	if( !error )
		cout << "\n***OK: Added " << nrows << " GRB rows (+ rate bounds)***\n"; 
	#endif
}

/*--------------------------------------------------------------------------*/

/* This private method adds max-delay constraint according
 * to SRP (Strictly-Rate-Proportional) formula 
 * using a SOCP-perspective formulation*/
void DCR_MFSP_SOCP_GRB::addDelayRows( int flowidx )
{
	int r, c, rcnt, nzcnt, error, nrows, nqrows;
	int *rmatind;
	double *rmatval, rhs, sense;
	
	//max delay constraint
	rcnt = 1;
	nzcnt = 1 + 2*numLinks;
	
	try{
	rmatind = new int[nzcnt];
	rmatval = new double[nzcnt];
	}catch(exception& e)
	{
		cout << "Standard exception: " << e.what() << endl;
	}
	
	rhs = Flows[flowidx].deadline;
	sense = GRB_LESS_EQUAL;
	
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
	
	error = GRBaddconstr(model, nzcnt, rmatind, rmatval, sense, rhs, NULL);
	if (error) 
	{
		printf("%s\n", GRBgeterrormsg(env));	
		throw( DCR::DCRException( "DCR_SFSP_SOCP_GRB::addDelayRows(): failed to create GUROBI rows" ) );
	}
	
	error = GRBupdatemodel(model);
	if (error)
	{
		printf("%s\n", GRBgeterrormsg(env));
		throw( DCR::DCRException( "DCR_SFSP_SOCP_GRB::addDelayRows(): failed to update GUROBI SOCP rows" ) );
	}
		
		
	delete [] rmatind;
	delete [] rmatval;
	
	
	//SOCP constraints: s_ij r_ij >= L x_ij^2
	int quadnzcnt, *quadrow, *quadcol;
	double *quadval;
	
    quadnzcnt = 2;
 
	try{
    quadrow = new int[quadnzcnt];
    quadcol = new int[quadnzcnt];
    quadval = new double[quadnzcnt];
	}catch(exception& e)
	{
		cout << "Standard exception: " << e.what() << endl;
	}
    
    rhs = 0;
	sense = GRB_GREATER_EQUAL;
    
     
   
    for(r = 0; r < numLinks; r++)
    {
		//s_ij
		quadrow[0] = numLinks*2 +1 + r + FOFFSET*flowidx;
		//r_ij
		quadcol[0] = numLinks + r + FOFFSET*flowidx;
		quadval[0] = 1;
		
			
		//- L x_ij^2
		quadrow[1] = r + FOFFSET*flowidx;
		quadcol[1] = r + FOFFSET*flowidx;
		quadval[1] = -MTU;
		
			
		error = GRBaddqconstr(model, 0, NULL, NULL, quadnzcnt, quadrow, quadcol, quadval, sense, rhs, NULL);
		if (error) 
		{
			printf("%s\n", GRBgeterrormsg(env));	
			throw( DCR::DCRException( "DCR_SFSP_SOCP_GRB::addDelayrows(): failed to create GUROBI SOCP rows" ) );
		}
		
		error = GRBupdatemodel(model);
		if (error)
		{
			printf("%s\n", GRBgeterrormsg(env));
			throw( DCR::DCRException( "DCR_SFSP_SOCP_GRB::addDelayRows(): failed to update GUROBI SOCP rows" ) );
		}
				
	}
		
	delete [] quadrow;
	delete [] quadcol;
	delete [] quadval;

	
	//SOCP constraints: t r_min >= \sigma	
	rhs = Flows[flowidx].burst;
	sense = GRB_GREATER_EQUAL;
	
	quadnzcnt = 1;
 
	try{
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
		
	error = GRBaddqconstr(model, 0, NULL, NULL, quadnzcnt, quadrow, quadcol, quadval, sense, rhs, NULL);
	if (error) 
	{
		printf("%s\n", GRBgeterrormsg(env));	
		throw( DCR::DCRException( "DCR_SFSP_SOCP_GRB::addDelayrows() cont.: failed to create GUROBI SOCP rows" ) );
	}
	
	error = GRBupdatemodel(model);
	if (error)
	{
		printf("%s\n", GRBgeterrormsg(env));
		throw( DCR::DCRException( "DCR_SFSP_SOCP_GRB::addDelayRows(): failed to update GUROBI SOCP rows" ) );
	}
		
	delete [] quadrow;
	delete [] quadcol;
	delete [] quadval;
	
	
	#ifdef DEBUG_GRB
	error = GRBgetintattr(model, "NumConstrs", &nrows);
	error = GRBgetintattr(model, "NumQConstrs", &nqrows);
	
	if( !error )
		cout << "\n***OK: Added " << nrows << " GRB rows (+ delay rows)***\n"; 
		cout << "\n***OK: Added " << nqrows << " GRB quadratic rows (+ delay rows)***\n"; 
	#endif
	
}

/*--------------------------------------------------------------------------*/

/* This private method adds max-delay constraint according
 * to SRP (Strictly-Rate-Proportional) formula 
 * using a SOCP-bigM formulation*/
void DCR_MFSP_SOCP_GRB::addDelayRowsBigM( int flowidx )
{
	int error;
	
	//NB. bigM must be large enough
	//double bigM = 100*Flows[flowidx].deadline;
	double a, b;
	a = sqrt(MTU);
	b = MTU/Flows[flowidx].rate;
	double bigM = fmax(a, b) +1;
		
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
		rmatind = new int[nzcnt];
		rmatval = new double[nzcnt];
	}catch(exception& e)
	{
		cout << "Standard exception: " << e.what() << endl;
	}
	
		
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
	
	error = GRBaddconstr(model, nzcnt, rmatind, rmatval, GRB_LESS_EQUAL, Flows[flowidx].deadline, NULL);
	if (error) 
	{
		printf("%s\n", GRBgeterrormsg(env));	
		throw( DCR::DCRException( "DCR_SFSP_SOCP_GRB::addDelayRows(): failed to create GUROBI rows" ) );
	}
	
	error = GRBupdatemodel(model);
	if (error)
	{
		printf("%s\n", GRBgeterrormsg(env));
		throw( DCR::DCRException( "DCR_SFSP_SOCP_GRB::addDelayRows(): failed to update GUROBI SOCP rows" ) );
	}
		
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
			
		error = GRBaddqconstr(model, 0, NULL, NULL, quadnzcnt, quadrow, quadcol, quadval, GRB_GREATER_EQUAL, MTU, NULL);
		if (error) 
		{
			printf("%s\n", GRBgeterrormsg(env));	
			throw( DCR::DCRException( "DCR_SFSP_SOCP_GRB::addDelayrows() cont.: failed to create GUROBI SOCP rows" ) );
		}
	
		error = GRBupdatemodel(model);
		if (error)
		{
			printf("%s\n", GRBgeterrormsg(env));
			throw( DCR::DCRException( "DCR_SFSP_SOCP_GRB::addDelayRows(): failed to update GUROBI SOCP rows" ) );
		}
	
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
		
	error = GRBaddqconstr(model, 0, NULL, NULL, quadnzcnt, quadrow, quadcol, quadval, GRB_GREATER_EQUAL, Flows[flowidx].burst, NULL);
	if (error) 
	{
		printf("%s\n", GRBgeterrormsg(env));	
		throw( DCR::DCRException( "DCR_SFSP_SOCP_GRB::addDelayrows() cont.: failed to create GUROBI SOCP rows" ) );
	}
	
	error = GRBupdatemodel(model);
	if (error)
	{
		printf("%s\n", GRBgeterrormsg(env));
		throw( DCR::DCRException( "DCR_SFSP_SOCP_GRB::addDelayRows(): failed to update GUROBI SOCP rows" ) );
	}
	
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
	
	error = GRBaddconstrs(model, rcnt, nzcnt, rmatbeg, rmatind, rmatval, sense, rhs, NULL);
	if (error) 
	{
		printf("%s\n", GRBgeterrormsg(env));	
		throw( DCR::DCRException( "DCR_SFSP_SOCP_GRB::addDelayRows(): failed to create GUROBI rows" ) );
	}
	
	error = GRBupdatemodel(model);
	if (error)
	{
		printf("%s\n", GRBgeterrormsg(env));
		throw( DCR::DCRException( "DCR_SFSP_SOCP_GRB::addDelayRows(): failed to update GUROBI SOCP rows" ) );
	}
	
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
	
	error = GRBaddconstrs(model, rcnt, nzcnt, rmatbeg, rmatind, rmatval, sense, rhs, NULL);
	if (error) 
	{
		printf("%s\n", GRBgeterrormsg(env));	
		throw( DCR::DCRException( "DCR_SFSP_SOCP_GRB::addDelayRows(): failed to create GUROBI rows" ) );
	}
	
	error = GRBupdatemodel(model);
	if (error)
	{
		printf("%s\n", GRBgeterrormsg(env));
		throw( DCR::DCRException( "DCR_SFSP_SOCP_GRB::addDelayRows(): failed to update GUROBI SOCP rows" ) );
	}
	
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
		rhs[r] = bigM; //M is set equal to the flow deadline
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
	
	error = GRBaddconstrs(model, rcnt, nzcnt, rmatbeg, rmatind, rmatval, sense, rhs, NULL);
	if (error) 
	{
		printf("%s\n", GRBgeterrormsg(env));	
		throw( DCR::DCRException( "DCR_SFSP_SOCP_GRB::addDelayRows(): failed to create GUROBI rows" ) );
	}
	
	error = GRBupdatemodel(model);
	if (error)
	{
		printf("%s\n", GRBgeterrormsg(env));
		throw( DCR::DCRException( "DCR_SFSP_SOCP_GRB::addDelayRows(): failed to update GUROBI SOCP rows" ) );
	}
	
	delete [] rmatbeg;
	delete [] rhs;
	delete [] sense;
	delete [] rmatind;
	delete [] rmatval;
		
}

/*--------------------------------------------------------------------------*/

/* This private method adds mutual capacity constraints in the SOCP model*/
/* \sum_{k \in FLows} r_ij^k <= c_ij \forall (i,j) \in Links*/
void DCR_MFSP_SOCP_GRB::addCapacityRows( void )
{	
	
	int r, c, rcnt, nzcnt, status, k, error;
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
		//FIXME: this is not exactly what the formula say, I should use c_ij(^k) ???
		rhs[r] = Links[r].capacity;
		sense[r] = GRB_LESS_EQUAL; 
		
		for(k = 0; k < numFlows; k++)
		{
			rmatind[c] = FOFFSET*k + r_offset + r;
			rmatval[c] = 1;
			c++;
		}
	}
	
	error = GRBaddconstrs(model, rcnt, nzcnt, rmatbeg, rmatind, rmatval, sense, rhs, NULL);
	if (error)
	{
		printf("%s\n", GRBgeterrormsg(env));
		throw( DCR::DCRException( "DCR_MFSP_SOCP_GRB::addCapacityRows(): failed to create GUROBI rows" ) );
	}
	
	error = GRBupdatemodel(model);
	if (error)
	{
		printf("%s\n", GRBgeterrormsg(env));
		throw( DCR::DCRException( "DCR_MFSP_SOCP_GRB::addCapacityRows(): failed to update GUROBI rows" ) );
	}
		
	delete [] rmatbeg;
	delete [] rhs;
	delete [] sense;
	delete [] rmatind;
	delete [] rmatval;	
	
}

/*--------------------------------------------------------------------------*/

/* This private method returns the number of entering arcs into node*/
int DCR_MFSP_SOCP_GRB::getNumInLinks(int node)
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
int DCR_MFSP_SOCP_GRB::getNumOutLinks(int node)
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
bool DCR_MFSP_SOCP_GRB::isLinkInNode(int node, int arc)
{
	if(Links[arc].endnode == node) return true;
	
	return false;
}

/*--------------------------------------------------------------------------*/

/* This private method checks if arc is exiting from node*/
bool DCR_MFSP_SOCP_GRB::isLinkOutNode(int node, int arc)
{
	if(Links[arc].startnode == node)return true;
	
	return false;
}

/*--------------------------------------------------------------------------*/

/* This private method checks if node is the source of flow*/
bool DCR_MFSP_SOCP_GRB::isFlowSource(int node, int flow)
{
	if(Flows[flow].sourcenode == node) return true;
	
	return false;
}

/*--------------------------------------------------------------------------*/

/* This private method checks if node is the sink of flow*/
bool DCR_MFSP_SOCP_GRB::isFlowSink(int node, int flow)
{
	if(Flows[flow].sinknode == node ) return true;
	
	return false;
}
/*--------------------------------------------------------------------------*/

/* This private method is used to copy array data (flows, links and nodes) from the load() function 
 * into the class as class members.*/
void DCR_MFSP_SOCP_GRB::copyDataArrays(DCRFlow *flows, DCRLink *links, DCRNode *nodes)
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
void DCR_MFSP_SOCP_GRB::getSol()
{
	int error;
	int i, k, c, nc;
	double *sol;

	/////////////////////////////////////////
	//GET OBJVAL
	//////////////////////////////////////////	
	error = GRBgetdblattr(model, "ObjVal", &objval);
	if (error)
	{		printf("%s\n", GRBgeterrormsg(env));
			throw( DCR::DCRException( "DCR_MFSP_SOCP_GRB::getSol: failed to get GUROBI objval" ) );
	}
	
	///////////////////////////////////////
	//GET SOLUTION
	////////////////////////////////////////
	
	GRBgetintattr(model, "NumVars", &(nc));	
	
	try
	{
		sol = new double[nc];		
	}catch(exception& e)
	{
		cout << "Standard exception: " << e.what() << endl;
	}
	
	//get the whole solution vector from cplex		
	error = GRBgetdblattrarray(model, "X", 0, nc, sol);
	if( error )
		throw( DCR::DCRException( "DCR_MFSP_SOCP_GRB::getSol(): failed to get GUROBI solution" ) );
	
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
		for(i = 0; i < numLinks; i++)
		{
			Rsol[c] = sol[k*FOFFSET + r_offset + i];
			c++;
		}
	}
		
	delete [] sol;
	
}

/*--------------------------------------------------------------------------*/
/* This private method deallocates dinamically allocated memory
 * */
void DCR_MFSP_SOCP_GRB::clean_up()
{
	int error;

	error = GRBfreemodel(model);
	if (error) 
	{
		printf("%s\n", GRBgeterrormsg(env));
		throw( DCR::DCRException( "DCR_MFSP_SOCP_GRB::destructor(): failed to free GRB model" ) );
	} 		

  	GRBfreeenv(env);
	
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
/*-------------------------- End File DCR_MFSP_SOCP_GRB.cpp-----------------*/
/*--------------------------------------------------------------------------*/
