#ifndef FLOWSIM_H_
#define FLOWSIM_H_

#include <vector>
#include "EventSimulator.h"
#include "DCR.h"
#include "OPTUtils.h"
#include "DCR.h"
#include "DCRGraph.h"
#include "DCR_MFSP_SOCP_CPX.h"
#include "DCR_MFSP_SOCP_GRB.h"
#include "DCR_SPT.h"
#include "DCR_INDI.h"

//namespace flowsim {

class FlowSim {
		
public:

	///////////////////////////
	//Public types
	///////////////////////////	
	
//FIXME: is a network stucture really useful ???
//FIXME: a struct or a class would be better ???
	struct SNetwork {
		int nnodes;
		int nlinks;
		int nflows;
		double MTU;
		DCR::DCRFlow *flows;
		DCR::DCRLink *links;
		DCR::DCRNode *nodes;
	};
	
	/*---------------------------------------------------------*/

	class FlowArrivalEvent: public EventSimulator::Event {
		//FIXME: are these members okay?
		FlowSim* fsim;
		int flowidx;
	public:
		FlowArrivalEvent(double time, FlowSim* fsim, int fidx);
		virtual void execute();
	};
	
	/*---------------------------------------------------------*/

	class FlowDepartureEvent: public EventSimulator::Event {
		//FIXME: are these Event members okay?
		FlowSim* fsim;
		int flowidx;
		int *path;
		double *rates;
		int nhops;
	public:
		FlowDepartureEvent(double time, FlowSim* fsim, int fidx, int *path, double *rates, int nhops);
		virtual void execute();
	};

	///////////////////////////
	//Constructor and destructor
	///////////////////////////	

	FlowSim(EventSimulator& esim, DCR *solver, SNetwork* net, char *logfile);
	~FlowSim();
	
	///////////////////////////
	//Public methods
	///////////////////////////	

	void run();
	void setSimulationId(const std::string& simulation_id);
	void setArrSeed(int s);
	void setDepSeed(int s);
	int  getArrSeed(void);
	int  getDepSeed(void);
	void setArrRate(double ar);
	double getArrRate(void);
	void setDepRate(double dr);
	double getDepRate(void);
	int getNumFailed(void);
	int getFlowCount(void);
	double getMeanTime(void);
	double getVarTime(void);
	double getMaxTime(void);

	///////////////////////////
	//Private members
	///////////////////////////	

private:

	std::string simulationId;
	EventSimulator& esim;
	DCR* solver;
	SNetwork *net;
	OPTrand ARRrandgen; //FIXME: is it right to have 2 different ???
	OPTrand DEPrandgen; //FIXME: is it right to have 2 different ???
	OPTrand flowrandgen;
	int arrseed;
	int depseed; 
	double arrrate;
	double deprate;
	int *scheduled;
	int nsched;
	int failed;	
	double maxtime;
	int flowcount;
	double sumtime;
	double sumvartime;
	char *logfile;
	
	/*---------------------------------------------------------*/

	void generateFlowArrivalEvent();	
	void generateFlowDepartureEvent(int fidx, int *path, double *rates, int nhops);
	void handleFlowArrivalEvent(int fidx);
	void handleFlowDepartureEvent(int fidx, int *path, double *rates, int nhops);

};

//} /* namespace flowsim */
#endif /* FLOWSIM_H_ */
