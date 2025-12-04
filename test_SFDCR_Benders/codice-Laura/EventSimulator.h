#ifndef EVENTSIMULATOR_H_
#define EVENTSIMULATOR_H_

#include <cassert>
#include <iostream>
#include <queue>
#include <stdexcept>
#include <fstream>


using namespace std;

class EventSimulator {

public:

	///////////////////////////
	//Public types
	///////////////////////////	
	
	class Event {
	public:
		const double time;

		Event(double time) : time(time) { }
		virtual ~Event() { }
		virtual void execute() = 0;
	};

	/*---------------------------------------------------------*/

	class EventComparison {
	public:
		EventComparison() { }
		bool operator() (const Event* lhs, const Event* rhs) const {
			return (lhs->time > rhs->time);
		}
	};

	///////////////////////////
	//Constructor
	///////////////////////////	

	EventSimulator() :
		current(0),
		max_time(0),
		stopRequested(false),
		comparator(),
		events(),
		step(0)	{ };

	///////////////////////////
	//Destructor
	///////////////////////////	

	~EventSimulator() {
		/*foreach (Event* e, events) {
			delete e;
		}*/
		std::vector<Event *>::iterator it;
		for ( it = events.begin(); it != events.end(); ) {	
			delete (*it);  //FIXME: is this really destroying the object ???
			it = events.erase(it);
		}
	};

	///////////////////////////////
	//Getter and setter methods
	///////////////////////////////	

	unsigned int getStep() { return step; }
	void setStep(unsigned int s) { step = s; }  
	double getTime() { return current; }
	double getMaxTime() { return max_time; }
	void setMaxTime(double t) { max_time = t; }

	///////////////////////////
	//Other public methods
	///////////////////////////	

	// the Event object will be deleted when removed from the queue
	void add(Event* e) {	
		assert(e->time >= current);
		push(e);
	}

	/*---------------------------------------------------------*/

	void run() {
		while (!stopRequested && current < max_time && !events.empty()) {
			Event* e = top();
			pop();
			current = e->time;
			
			//debug
			//cout << "\n simulator: current time is " << current << endl;
			//getchar();
			
			e->execute();
			delete e;
			step++;
		}
	}

	/*---------------------------------------------------------*/

	// terminates the execution
	void stop() {
		stopRequested = true;
	}

	/*---------------------------------------------------------*/

	void increaseTime(double dt) {
		assert(dt >= 0);
		double newt = current + dt;
		if (newt == current)
			throw std::underflow_error("underflow in time increase");
		current = newt;
		assert(current <= top()->time);
	}

	/*---------------------------------------------------------*/

	void push(Event* e) {
		events.push_back(e);
		std::push_heap(events.begin(), events.end(), comparator);		
	}

	/*---------------------------------------------------------*/

	Event* top() {
		return events.front();
	}

	/*---------------------------------------------------------*/

	void pop() {
		std::pop_heap(events.begin(), events.end(), comparator);
		events.pop_back();
	}
	
///////////////////////////
//Protected members
///////////////////////////	

protected:

//FIXME: use a heap structure
//	typedef std::priority_queue<Event*, std::vector<Event*>, EventComparison> PriorityQueue;
	double current;
	double max_time;
	bool stopRequested;
	EventComparison comparator;
	std::vector<Event*> events;
	unsigned int step;

};

#endif /* EVENTSIMULATOR_H_ */
