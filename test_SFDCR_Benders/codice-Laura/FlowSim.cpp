#include "FlowSim.h"

#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <cassert>
#include <cmath>

//namespace flowsim {

///////////////////////////
//FlowArrivalEvent
///////////////////////////	

FlowSim::FlowArrivalEvent::FlowArrivalEvent(double time, FlowSim* fsim, int fidx) :
	EventSimulator::Event(time),
	fsim(fsim),
	flowidx(fidx){
}

/*---------------------------------------------------------*/

void FlowSim::FlowArrivalEvent::execute() {
	fsim->handleFlowArrivalEvent(flowidx);
}

///////////////////////////
//FlowDepartureEvent
///////////////////////////	

FlowSim::FlowDepartureEvent::FlowDepartureEvent(double time, FlowSim* fsim, int fidx, int *path, double *rates, int nhops) :
	EventSimulator::Event(time),
	fsim(fsim),
	flowidx(fidx),
	path(path),
	rates(rates),
	nhops(nhops) {
}

/*---------------------------------------------------------*/

void FlowSim::FlowDepartureEvent::execute() {
	fsim->handleFlowDepartureEvent(flowidx, path, rates, nhops);
}

///////////////////////////
//FlowSim
///////////////////////////	

FlowSim::FlowSim(EventSimulator& esim, DCR *solver, SNetwork* net, char *logfile) :
	simulationId(""),
	esim(esim),
	solver(solver),
	net(net),
	ARRrandgen(),
	DEPrandgen(),
	arrseed(1),
	depseed(1),
	arrrate(1.0),
	deprate(1.0),
	failed(0),
	maxtime(0),
	scheduled(),
	flowcount(0),
	sumtime(0),
	sumvartime(0),
	nsched(0),
	logfile(logfile){
}

FlowSim::~FlowSim() {
	
	delete solver;

	if (net->flows) {
		for(int i = 0; i < net->nflows; i++) {
			delete [] net->flows[i].caps;
			delete [] net->flows[i].costs;
		} 
	}	
	delete net->flows;
	delete net->links;
	delete net->nodes;
	delete net;
	delete [] scheduled;
	
}

///////////////////////////
//Public methods
///////////////////////////	

void FlowSim::setSimulationId(const std::string& simulation_id) {
	this->simulationId = simulation_id;
}

/*---------------------------------------------------------*/

void FlowSim::setArrSeed(int arrseed) {
	this->arrseed = arrseed;
}

/*---------------------------------------------------------*/

void FlowSim::setDepSeed(int depseed) {
	this->depseed = depseed;
}

/*---------------------------------------------------------*/

int FlowSim::getArrSeed(void) {
	return arrseed;
}

/*---------------------------------------------------------*/

int FlowSim::getDepSeed(void) {
	return depseed;
}

/*---------------------------------------------------------*/

void FlowSim::setArrRate(double arrrate) {
	this->arrrate = arrrate;
}

/*---------------------------------------------------------*/

double FlowSim::getArrRate(void) {
	return arrrate;
}

/*---------------------------------------------------------*/

void FlowSim::setDepRate(double deprate) {
	this->deprate = deprate;
}

/*---------------------------------------------------------*/

double FlowSim::getDepRate(void) {
	return deprate;
}

/*---------------------------------------------------------*/

int FlowSim::getNumFailed(void) {
	return failed;
}

double FlowSim::getMeanTime(void) {
	return sumtime/flowcount;
}

double FlowSim::getVarTime(void) {
	double sqrmean;

	sqrmean = pow(sumtime/flowcount, 2);

	return sumvartime/flowcount - sqrmean; 
}

double FlowSim::getMaxTime(void) {
	return maxtime;
}

int FlowSim::getFlowCount(void) {
	return flowcount;
}
/*---------------------------------------------------------*/

//FIXME: what happens if I call run more than once ???
void FlowSim::run() {

	if (net->nflows == 0)
		throw std::invalid_argument("no flows to simulate!");
	
	//initialize internal structures
	ARRrandgen.srand(arrseed);
	DEPrandgen.srand(depseed);	
	flowrandgen.srand(1);
	failed = 0;
	nsched = 0;
	scheduled = new int [net->nflows];
	for(int i = 0; i < net->nflows; i++)
		scheduled[i] = 0;
		
	//starts simulation with an arrival event
	generateFlowArrivalEvent();
	
	esim.run();		
}

///////////////////////////
//Private members
///////////////////////////	

/* FlowArrivalEvent callback */
void FlowSim::handleFlowArrivalEvent(int fidx) {	

	int nhops;
	int *X = new int[net->nnodes-1];
	double *R = new double[net->nnodes-1];	
	int status;	
	double duration;

	// solve path computation and resource allocation
	solver->DCRloadProblem(net->nnodes, net->nlinks, 1, net->flows + fidx, net->links, net->nodes, net->MTU, DCR::SRP);	
	
	//get time for stats
	solver->DCRsetTime(true);
	solver->DCRstartTime();				
	
	status = solver->DCRsolve();
	
	//get time for stats
	solver->DCRstopTime();			
	duration = solver->DCRgetTime();
	
	//update stats
	flowcount ++;
	sumtime += duration;
	sumvartime += pow(duration, 2);
	if(duration > maxtime) maxtime = duration;

	if(status) 
	{
		//update stats
		failed++;	
		
		delete [] X;
		delete [] R;
		
		//update scheduled events
		scheduled[fidx] = 0;
		nsched --;

		//update log 
		ofstream outfile;
		outfile.open(logfile, ios::app);
		if(status==4)
		{
		outfile << "\n***" << fidx << "***\t" << net->flows[fidx].sourcenode << "\t" << net->flows[fidx].sinknode << "\t" << duration << "\t" << nsched << "\t0\t" << esim.getTime(); 
		}
		else 
		{
		outfile << "\n" << fidx << "\t" << net->flows[fidx].sourcenode << "\t" << net->flows[fidx].sinknode << "\t" << duration << "\t" << nsched << "\t0\t" << esim.getTime(); 
		}
		outfile.close();
	}			
	else{
		//get solution
		nhops = solver->DCRgetUBSolPath(0, 0, X, R);	
								
		// occupy resources (resources)
		for(int i = 0; i < nhops; i++)
		{
			for(int k = 0; k < net->nflows; k++)
			{
				net->flows[k].caps[X[i]] -= R[i];
			}
		}

		//update log 
		ofstream outfile;
		outfile.open(logfile, ios::app);
		outfile << "\n" << fidx << "\t" << net->flows[fidx].sourcenode << "\t" << net->flows[fidx].sinknode << "\t" << duration << "\t" << nsched << "\t1\t" << esim.getTime();
		
		#ifdef PRINT_SOL
		outfile << "\t" << nhops;
		for(int i = 0; i < nhops; i++)
			outfile << "\t" << R[i];
		#endif 
		outfile.close();
		
		// generate new departure event (resources)
		generateFlowDepartureEvent(fidx, X, R, nhops); //FIXME: is this way to deal with X and R arrays okay???
	}
	
	// generate new arrival event
	generateFlowArrivalEvent();	
}

/*---------------------------------------------------------*/

/* FlowDepartureEvent callback */
void FlowSim::handleFlowDepartureEvent(int fidx, int *path, double *rates, int nhops) {
	
	//release network resources
	for(int i = 0; i < nhops; i++)
	{
		for(int k = 0; k < net->nflows; k++)
			{
				net->flows[k].caps[path[i]] += rates[i];
			}
	}
	
	delete [] path;
	delete [] rates;
	
	//update scheduled flows
	scheduled[fidx] = 0;
	nsched --;
}

/*---------------------------------------------------------*/

void FlowSim::generateFlowDepartureEvent(int fidx, int *path, double *rates, int nhops) {
	
	//FIXME: is this okay ???
	//schedule event according to exponential distribution
	double u = DEPrandgen.rand();
	double t = - log(1-u)/(deprate) + esim.getTime();
	//double t = DEPrandgen.rand() + esim.getTime();					
	esim.add(new FlowDepartureEvent(t, this, fidx, path, rates, nhops));
	
	//debug
	//cout << "\n generate flow departure event for flow " << fidx << " at time " << t << endl;

}

/*---------------------------------------------------------*/

void FlowSim::generateFlowArrivalEvent() {
	
	//pick a random flow among the unscheduled ones
	int rfidx = floor((net->nflows - nsched)*flowrandgen.rand());	
	
	int findx;
	
	int c = 0;	
	for(int i = 0; i < net->nflows; i++)
	{
		if(!scheduled[i])
		{
			if(c == rfidx)
			{
				findx = i; break;
			}//if 
			c++;
		}//if
	}
	
	//update scheduled flows
	scheduled[findx] = 1;
	nsched ++;
	
	//FIXME: is this okay ???
	//schedule event according to exponential distribution
	double u = ARRrandgen.rand();
	double t = - log(1-u)/(arrrate) + esim.getTime();
	//double t = ARRrandgen.rand() + esim.getTime();						
	esim.add(new FlowArrivalEvent(t, this, findx));
	
	//debug
	//cout << "\n generate flow arrival event for flow " << findx << " at time " << t << endl;
}
//} /* namespace flowsim */
