#include "DCR.h"
#include "DCRGraph.h"
#include "DCR_MFSP_SOCP_CPX.h"
#include "DCR_MFSP_SOCP_GRB.h"
#include "DCR_SPT.h"
#include "DCR_INDI.h"
#include "DCR_3PRONG.h"
#include "EventSimulator.h"
#include "FlowSim.h"
#include "math.h"


using namespace std;

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
}


int main(int argc, char **argv)
{
	
	//set up command line parameters
	/////////////////////////
	char *filein; 
	char *fileout;  
	ofstream outfile;
	char format = 'm'; //instance format (e.g.,'m'="mnetgen", 's'="mulgen") 
	double percap = 0.0;
	int seed = 1;
	char alg = 's';//can be: ERA = e, INDI = i, 3-PRONGED = p, SOCP = s
	char heurtype = '1';// for ERA can be: 1 (ERA-I), 2 (ERA-H) 3 (ERA-H-ML); for INDI can be 1(INDI-WS), 2(INDI-SW)
	char solver = 'c';// for SOCP, can be: CPLEX = c, GUROBI = g;
	int wflows[1];
	
	switch( argc ) {
	
	case( 9 ): solver = argv[8][0];
	
	case( 8 ): heurtype = argv[7][0];
	
	case( 7 ): alg = argv[6][0];
	
	case( 6 ): seed = atoi(argv[5]);
	
	case( 5 ): percap = atof(argv[4]); 
	
	case( 4 ): format = argv[3][0];
	
	case( 3 ): fileout = argv[2];
			   filein = argv[1]; 
			   break;

	default:   cerr << "Usage: " << argv[ 0 ]
			   << " file_in file_out [format, percap, seed, alg, heurtype, solver]" << endl
			   << " format (default " << format << ")" << endl
			   << " percap (default " << percap << ")" << endl
			   << " seed (default " << seed << ")" << endl
			   << " alg (default " << alg << ")" << endl
			   << " heurtype (default " << heurtype << ")" << endl
			   << " solver (default " << solver << ")" << endl;
			   return( 1 );
	}
/*
	cout << "\n INPUT: " 
	<< "\n instance " << filein 
	<< "\n file-out " << fileout 
	<< "\n format " << format 
	<< "\n percap " << percap
	<< "\n seed " << seed 
	<< "\n alg " << alg  
	<< "\n heurtype " << heurtype 
	<< "\n solver " << solver << endl; 
*/	
	
	try
	{	
		
		//set up network data
		/////////////////////////
		int nnodes;
		int nlinks;
		int nflows;
		double MTU;
		DCR::DCRFlow *flows = 0;
		DCR::DCRLink *links = 0;
		DCR::DCRNode *nodes = 0;
		ReadData(filein, format, &(nnodes), &(nlinks), &(nflows), &(MTU), &(flows), &(links), &(nodes), percap, 0, seed);					
		FlowSim::SNetwork *net = new FlowSim::SNetwork();
		net->nnodes = nnodes;
		net->nlinks = nlinks;
		net->nflows =  nflows;
		net->MTU  = MTU;
		net->flows = flows;
		net->links = links;
		net->nodes = nodes;	
		
		//stats
		//////////////////////////
		double duration;
		int status;
		int numinf = 0;
			
		int stat;
		int ninf = 0;
		double obj;
		double sumtime = 0;
		double meantime;
		double maxtime = 0;
	
		int stat_alg;
		int ninf_alg = 0;
		double obj_alg;
		double sumdiff_alg = 0;
		double meandiff_alg;
		double maxdiff_alg = 0;
		double sumtime_alg = 0;
		double meantime_alg;
		double maxtime_alg = 0;
		int count = 0;

		//set up solver 
		////////////////////////
		
		DCR_MFSP_SOCP_CPX socp;
	    DCR_INDI erai;
		DCR * algo;

		//create solvers objects
		//if(alg == 'e')  {
		//	algo = new DCR_SPT;
		//	dynamic_cast<DCR_SPT*>(algo)->DCRsetHeur(heurtype);
		//}
		if(alg == 'i')  {
			algo = new DCR_INDI;
			dynamic_cast<DCR_INDI*>(algo)->DCRsetHeur(heurtype);
		}
		//else if(alg == 'p')  {
			//algo = new DCR_3PRONG;
			//dynamic_cast<DCR_3PRONG*>(algo)->DCRsetHeur(heurtype);
			//dynamic_cast<DCR_3PRONG*>(algo)->DCRsetSolver(solver);
		//}
		else if(alg == 's' && solver == 'c') {algo = new DCR_MFSP_SOCP_CPX;}		
		else if(alg == 's' && solver == 'g') {algo = new DCR_MFSP_SOCP_GRB;}
		else { cout << "\n Wrong parameters, exit... " << endl; return 1; }

		algo = new DCR_MFSP_SOCP_CPX;

		for(int fi = 0; fi < net->nflows; fi++)
		{
			
			//cout << "\n flow: " << fi << endl;
			
			//load algorithms
			//wflows[0]= fi;
			socp.DCRloadProblem(net->nnodes, net->nlinks, 1, net->flows + fi, net->links, net->nodes, net->MTU, DCR::SRP);	
			algo->DCRloadProblem(net->nnodes, net->nlinks, 1, net->flows + fi, net->links, net->nodes, net->MTU, DCR::SRP);	
			erai.DCRloadProblem(net->nnodes, net->nlinks, 1, net->flows + fi, net->links, net->nodes, net->MTU, DCR::SRP);	

			//ERA-I to check feasibility
			//////////////////////////////////
			
			erai.DCRsetHeur('1');
			status = erai.DCRsolve();
			if(status) { 
				numinf++; 
			}

			//SOCP
			//////////////////////////////////
			socp.DCRsetFsbEps(1e-06);
			socp.DCRsetOptEps(1e-04);
			//solve
			socp.DCRsetTime(true);
			socp.DCRstartTime();
			stat = socp.DCRsolve();
			socp.DCRstopTime();
			duration = socp.DCRgetTime();
			obj = socp.DCRgetObj();	

			//time
			sumtime += duration;
			if(duration > maxtime) maxtime = duration; 
			//solution
			if(stat) {
				ninf++; 
				obj = -1; 
				if(!status)
				{
					outfile.open("mylogSOCP.txt",ios::app);
					outfile << "\n SOCP failed but ERA-I did not for " << filein << " flow " << fi << " load " << percap << endl;
					outfile.close();
				}
			}
			else {
				if(status) 
				{
					outfile.open("mylogSOCP.txt",ios::app);
					outfile << "\n ERA-I failed but SOCP did not for " << filein << " flow " << fi << " load " << percap << endl;
					outfile.close();
				}
				else
				obj = socp.DCRgetObj();	
			}

			//ALGO
			////////////////////////////////////////
			//solve
			algo->DCRsetTime(true);			
			algo->DCRstartTime();			
			stat_alg = algo->DCRsolve();
			algo->DCRstopTime();
			duration = algo->DCRgetTime();
			obj = algo->DCRgetObj();	
			cout << obj << "\n";
			//time
			sumtime_alg += duration;
			if(duration > maxtime_alg) maxtime_alg = duration; 
			//solution
			if(stat_alg) {
				ninf_alg++;
				if(!status && alg != 'e')
				{
					outfile.open("mylogALGO.txt",ios::app);
					outfile << "\n ALGO failed but ERA-I did not  for " << filein << " algo " << alg << " flow " << fi << " load " << percap << endl;
					outfile.close();
				}
			} 
			else
			{
				if(status)
				{
					outfile.open("mylogALGO.txt",ios::app);
					outfile << "\n ERA-I failed but ALGO did not  for " << filein << " algo " << alg << " flow " << fi << " load " << percap << endl;
					outfile.close();

				}
				if(obj == -1)
				{
					outfile.open("mylogALGO.txt",ios::app);
					outfile << "\n SOCP failed but ALGO did not for " << filein << " algo " << alg << " flow " << fi << " load " << percap << endl;
					outfile.close();
				}
				else
				{
					obj_alg = algo->DCRgetObj(); 
					sumdiff_alg += (obj_alg - obj)/obj_alg; 
					if((obj_alg - obj)/obj_alg > maxdiff_alg) maxdiff_alg = (obj_alg - obj)/obj_alg;
					count ++;
				}
			}
		}
		
		meantime = sumtime/net->nflows;	
		meandiff_alg = sumdiff_alg/count;
		meantime_alg = sumtime_alg/net->nflows;
		
		outfile.open(fileout, ios::app);				
		outfile << "\n" << net->nflows << "\t" << numinf << "\t" << meantime << "\t" << maxtime << "\t" << ninf << "\t" << meantime_alg << "\t" << maxtime_alg << "\t" << meandiff_alg << "\t"  << maxdiff_alg << "\t" << ninf_alg;
		outfile.close();

		
	} catch(exception& e) { cout << "\nAn exception: \n" << e.what() << endl; }

	return 0;
	
	}
