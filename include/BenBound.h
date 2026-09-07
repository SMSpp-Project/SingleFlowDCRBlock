/*--------------------------------------------------------------------------*/
/*-------------------------- File BenBound.h --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class BenBound, which computes primal and dual bounds
 * for Single-Flow, Single-Path Delay-Constrained Routing (DCR) problems by
 * means of a Benders-like cutting-plane scheme applied to the (scalar)
 * Lagrangian dual function of the minimum-rate parameter r_min.
 *
 * The DCR problem, once the "burst size" / "peak rate" delay term is
 * dualized, reduces to the minimization over r_min > 0 of a one-dimensional,
 * in general nonconvex, piecewise-smooth function d( r_min ), each value of
 * which is obtained by solving a Lagrangian relaxation (via
 * DCRLagrangianSolver, which in turn relies on a Shortest Path Tree
 * computation) for the corresponding fixed r_min. BenBound repeatedly
 * evaluates d( . ) at a sequence of trial points and, exploiting the
 * (sub)gradient information obtained at each evaluation, builds a
 * piecewise-linear model of d( . ) out of the tangent/secant lines ("cuts")
 * collected so far; this is exactly the classical Benders/Kelley
 * cutting-plane master problem, here solved in closed form for the
 * univariate case by LineSearch(). The optimal value of the piecewise-linear
 * master problem provides a valid *lower* bound on the true optimum, while
 * the (feasible) function evaluations directly provide valid *upper* bounds;
 * the two bounds are refined and used to fathom the search until they meet
 * (within a tolerance) or a limit on the number of iterations is reached.
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
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef BENBOUND_H
#define BENBOUND_H

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "DCRLagrangianSolver.h" //serve?
#include "DCR.h"

#include <vector>


using namespace std;
//using namespace OPTtypes_di_unipi_it;
/*--------------------------------------------------------------------------*/
/*--------------------------- CLASS BenBound --------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// Benders-like bounding procedure for Single-Flow, Single-Path DCR problems
/** BenBound computes primal (upper) and dual (lower) bounds on the optimal
 * value of a Single-Flow, Single-Path Delay-Constrained Routing (DCR)
 * instance. The instance data (topology, link capacities/costs/delays, node
 * delays and the flow's source/sink/burst/rate/deadline) is loaded once via
 * LoadProblem(), after which Solve() runs the bounding procedure and the
 * get*() methods report the results (optimal/approximate objective value,
 * best primal solution found, iteration counters, elapsed time, ...).
 *
 * Internally, the DCR problem is attacked by dualizing (à la Lagrange) the
 * nonlinear part of the delay constraint that depends on the flow's burst
 * size and on a single scalar parameter r_min (a lower bound on the rate
 * reserved on every link of the path, which also determines the
 * corresponding Lagrangian multiplier lambda( r_min ) in closed form). For
 * any fixed r_min the resulting Lagrangian subproblem is solved by a
 * DCRLagrangianSolver, which reduces it to (a small number of) Shortest Path
 * Tree computations; the outcome is both a value d( r_min ) of the
 * (concave/convex-mixed) dual function and one or two supporting lines
 * ("cuts") of it. BenBound organizes the domain of r_min (the interval
 * where the Lagrangian subproblem is feasible, computed by Limitrmin()) into
 * a partition of subintervals Q[], each carrying the set of cuts collected
 * so far that are relevant to it; LineSearch() finds, in closed form, the
 * point where the current piecewise-linear upper envelope of the cuts is
 * minimized (this is the classical Benders/Kelley cutting-plane master
 * problem, solved exactly because it is one-dimensional and piecewise
 * linear), which both gives the current global lower bound and selects the
 * next r_min at which the Lagrangian subproblem is (re)solved. The whole
 * process, kicked off by Inizial(), iterates until the primal and dual
 * bounds (BestUB/BestLB, respectively tracking the best value of d( . ) ever
 * found and the value of the master problem) coincide up to a relative
 * tolerance, or a maximum number of iterations is exceeded. */

class BenBound
{

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

	public:

/*--------------------------------------------------------------------------*/
/*--------------------------- Inf() and Eps() ------------------------------*/
/*--------------------------------------------------------------------------*/
/** Very small class to simplify extracting the "+ infinity" value for a
    basic type; just use Inf<type>(). */

 template <typename T>
  class Inf {
   public:
  Inf() {}
  operator T() { return( std::numeric_limits<T>::max() ); }
  };

/*--------------------------------------------------------------------------*/
/** Very small class to simplify extracting the "machine epsilon" for a
    basic type; just use Eps<type>(). */

 template <typename T>
  class Eps {
   public:
  Eps() {}
  operator T() { return( std::numeric_limits<T>::epsilon() ); }
  };
	
/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public types
 *  @{ */

   /// status of the bounding procedure
   enum BndrStat
   {
	  OK, 	///<  solver found a feasible solution
	  Infeasible,   ///<  problem infeasible
	  Unbounded,    ///<  problem unbounded
	  Error         ///<  an error occurred
   };

/*--------------------------------------------------------------------------*/
   /// a subinterval of the domain of r_min, with the cuts relevant to it
   /** SubInterval represents one piece of the partition of the (feasible)
    * domain of the scalar parameter r_min into subintervals, each one
    * carrying the set of "cuts" (supporting lines of the Lagrangian dual
    * function d( r_min ), see the general notes of BenBound) that have been
    * collected so far and that are candidates for defining the piecewise-
    * linear upper envelope over that specific subinterval. Q[ 0 ] is a
    * "dummy" sentinel interval only used to store the left endpoint of the
    * overall domain; Q[ 1 ], ..., Q[ Q.size() - 1 ] are the actual
    * subintervals, listed in increasing order of their right endpoint
    * rmin. */
   struct SubInterval //is it necessary to save the two optimal cuts each time?
   {
     double rmin;    ///< right endpoint of the subinterval

     double inter;   ///< abscissa of the current candidate optimal point in
                      /// this subinterval (intersection of the two cuts
                      /// pCut/mCut, or of a cut with the subinterval border)

     double interVal; ///< ordinate of inter, i.e., the local lower bound on
                       /// d( . ) given by the cuts' intersection

     double Val;      ///< value of d( inter ), i.e., the local upper bound
                       /// on the optimum obtained by (re)solving the
                       /// Lagrangian subproblem at r_min == inter

     int bestCutpos;  ///< position (in Cuts) of the best cut, used when
                       /// infeasflag == 1 (all cuts in this subinterval have
                       /// positive slope, so no interior optimum exists yet)

     int branchedflag; ///< if nonzero, this subinterval has been fathomed
                        /// (its lower bound is already worse than the best
                        /// known upper bound BestUB, so it cannot improve it)

     int solflag;     ///< if nonzero, the Lagrangian subproblem has already
                       /// been (re)solved at r_min == inter

     int infeasflag;  ///< if nonzero, all the cuts collected so far in this
                       /// subinterval have positive slope (the subinterval
                       /// still needs to be restricted towards feasibility)

     DCRLagrangianSolver::LinearCut pCut; ///< the cut with nonnegative slope
                                           /// that, together with mCut,
                                           /// currently defines the optimum
     DCRLagrangianSolver::LinearCut mCut; ///< the cut with negative slope
                                           /// that, together with pCut,
                                           /// currently defines the optimum

     vector<DCRLagrangianSolver::LinearCut> Cuts; ///< all the cuts currently
                                                    /// available in this
                                                    /// subinterval
   };

/** @} ---------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor
 *  @{ */

	/// constructor of BenBound: initializes all data members to "empty"
	/** Void constructor: gives default values to all the data members of
	 * the class (in particular, no problem is loaded and the internal
	 * status is set to Error), so that LoadProblem() must be called before
	 * Solve() can be invoked. */

	BenBound(void);

/** @} ---------------------------------------------------------------------*/
/*---------------------------------- METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Loading the data of the problem
 *  @{ */

	/// loads a whole new DCR instance (topology, data and flow)
	/** Loads a complete Single-Flow DCR instance: the network topology
	 * (nnodes nodes, nlinks links), the per-link data (capacity, speed,
	 * delay, cost) in links, the per-node delay in nodes, the transmission
	 * unit size mtu, and the flow to be routed (source, sink, burst, rate
	 * and deadline) in flow. Any previously loaded instance is discarded. */

	void LoadProblem(int nnodes, int nlinks, DCR::DCRFlow flow, vector<DCR::DCRLink> links, vector<DCR::DCRNode> nodes, double mtu);

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	/// changes only the flow to be routed, keeping the topology
	/** Reloads only the flow data (source, sink, burst, rate, deadline) and
	 * the corresponding per-link costs, without re-reading the network
	 * topology; this is meant to be used after a first call to the other
	 * LoadProblem() overload, to solve a sequence of instances that only
	 * differ in the flow to be routed. */
    //FIXME: add a public method to modify only the flow
    // and one to modify myparam.
    void LoadProblem(DCR::DCRFlow flow);

/** @} ---------------------------------------------------------------------*/
/*----------------------------------SOLVE-----------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Solving the problem
 *  @{ */

	/// runs the Benders-like bounding procedure
	/** Runs the cutting-plane bounding procedure described in the general
	 * notes of the class: it repeatedly evaluates the Lagrangian dual
	 * function d( r_min ) (via Inizial() first, and then a loop driven by
	 * LineSearch()), refining the primal bound BestUB and the dual bound
	 * BestLB until they coincide (up to a relative tolerance) or the
	 * maximum number of iterations is reached. After the call, getUB(),
	 * getLB(), getObjVal(), getSolution() and the other get*() methods
	 * report the outcome. */

	void Solve();

/** @} ---------------------------------------------------------------------*/
/*----------------------------------GET RESULTS------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the results
 *  @{ */

    /// returns the value of d( . ) at the last evaluated r_min
    double getObjVal();
    /*<returns the optimal value found*/

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    /// returns the current cutting-plane approximation of the optimum
    double getApproxVal();
    /*<returns the approximation for this optimal value*/

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    /// returns the value of r_min at which the optimum was (best) found
    double getr_min();
    /*<optimal r_min*/

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    /// returns the value of the best primal-feasible solution found so far
    double getHeurVal();
    /*<value of the best primal-feasible solution found so far*/

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    /// returns the number of Benders (outer) iterations performed
    int getNumIterationBender();
    /*<number of points visited while solving the Benders procedure*/

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    /// returns the total number of Lagrangian (inner) iterations performed
    int getNumIterationLagr();
    /*<total sum of the number of Lagrangian iterations*/

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    /// returns true if at least one subinterval turned out to be convex
    /** Returns true if, during Solve(), at least one subinterval was found
     * for which the newly generated left and right cuts coincide, which
     * signals a (locally) convex piece of the dual function d( . ) and is
     * handled differently (see Solve()) from the general nonconvex case. */
    bool IsConvexInteration();
    /*<returns true if there was at least one convex iteration*/

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    /// returns d( . ) evaluated at the abscissae computed by getxPlotData()
    /** Debug/plotting helper: (re)solves the Lagrangian subproblem at each
     * of the numPoints abscissae returned by getxPlotData() (which must
     * therefore be called first) and returns the corresponding values of
     * d( r_min ), to be used together with getxPlotData() to plot (an
     * approximation of) the dual function over the first subinterval. */
    vector<double> getyPlotData();

    /// returns a discretization of the domain of the first subinterval
    /** Debug/plotting helper: calls Inizial() and returns numPoints equally
     * spaced abscissae spanning the first subinterval Q[ 1 ]. */
    vector<double> getxPlotData();

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    /// returns the elapsed running time, in seconds
    double getTime() const;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    /// returns the status of the bounding procedure
    BndrStat getStat();

/** @} ---------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS --------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Handling the timer
 *  @{ */

    /// sets a time limit for the solver, in seconds
    void DCRsetTimeLimit(long secs);

    /// if timeON is true creates/resets the timer, if false deletes it
    void DCRsetTime(bool timeON);

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    /// starts (or restarts) the timer, if it has been set
    void DCRstartTime();

    /// stops the timer, if it has been set
    void DCRstopTime();

/** @} ---------------------------------------------------------------------*/
/*----------------------------------GET RESULTS------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the primal solution
 *  @{ */

    /// returns the rate reserved on link i in the best solution found
    double getSolution(int i);
    /*<best primal-feasible solution found so far*/

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    /// returns the best available upper bound on the optimal value
    /** Returns min( BestUB , HeurVal ), i.e., the best between the best
     * value of d( . ) found while exploring the cutting-plane tree and the
     * value of the best primal-feasible (heuristic) solution found. */
    double getUB();

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    /// returns the best available lower bound on the optimal value
    /** Returns max( 0 , min( getUB() , BestLB ) ), i.e., the value of the
     * cutting-plane master problem, clipped from above by the current best
     * upper bound and from below by 0 (the DCR objective, being a sum of
     * nonnegative arc costs times nonnegative rates, cannot be negative). */
    double getLB();

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    /// returns the number of links of the loaded instance
    int get_Links(){ return(numLinks); };

    /// returns a pointer to the (sorted) array of link capacities
    double* get_caps(){ return(caps);}

/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @} ---------------------------------------------------------------------*/
/** @name Destructor
 *  @{ */

	/// destructor of BenBound: releases all dynamically allocated memory

	~BenBound();

/** @} ---------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/	


	private:

    // instance data - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    DCR::DCRFlow Flow;    ///< the flow to be routed
    DCR::DCRNode* Nodes;  ///< array of numNodes node data (node delays)
    DCR::DCRLink* Links;  ///< array of numLinks link data (capacity, speed,
                           /// delay, cost)
    int numNodes;         ///< number of nodes of the network
    int numLinks;         ///< number of links of the network
    double MTU;           ///< maximum transmission unit (packet size)
    double eps;           ///< relative tolerance used throughout the class

    // bounds and solution status - - - - - - - - - - - - - - - - - - - - - -

    double ObjVal;        ///< value of d( . ) at the last evaluated r_min
    int solvedflag;       ///< nonzero once the bounding procedure has met
                           /// its stopping criterion
    double rmin; //rmin corresponding to the optimal solution.
    double HeurVal; //value of the heuristic
    double ApproxVal; //value of the approximation.
    double mystep;        ///< step used while restricting towards a
                           /// feasible point of the Lagrangian subproblem
    double BestUB;        ///< best (smallest) value of d( . ) found so far
    double BestLB;        ///< value of the piecewise-linear master problem
                           /// at the current candidate optimum (valid lower
                           /// bound on the true optimal value)
    double crit_capc;     ///< "critical" link capacity value (see
                           /// Limitrmin()) used to restrict the initial
                           /// search interval for r_min
    int counter_ite_Ben; //counter for the number of Benders iterations.
    int counter_ite_Lag; //counter for the total number of Lagrangian iterations.
    bool is_convex_iteration; //flag to tell whether at least one convex cut was found

    double SOL_VALUE;     ///< (unused) cached copy of the optimal value
    vector<double> SOLUTION; ///< best primal-feasible rate solution found,
                              /// one entry per link (see getSolution())

    OPTtimers *timer;  ///< timer
    long tlimit; ///< time limit

    //data for plotting the function
    vector<double> xPlot; ///< abscissae computed by getxPlotData()
    int numPoints;        ///< number of points used by get[xy]PlotData()

    double* caps; //vector of the arc capacities.

    BndrStat BenStat;     ///< overall status of the bounding procedure

    vector<double> SPLabels; //dual values of SP
    DCRLagrangianSolver lagSol; ///< the Lagrangian subproblem solver used to
                                 /// evaluate d( r_min ) and obtain the cuts
    vector<SubInterval> Q; //subintervals of V considered.

    double myparam; //parameter in (0,1) that selects the point in cases of infeasibility

    // private methods - - - - - - - - - - - - - - - - - - - - - - - - - - -

    /// initializes the search: computes the feasible range of r_min and
    /// the first cut(s), possibly detecting an immediate optimum
    void Inizial();

    /// computes [ lower , upper ] bounds on the values of r_min for which
    /// the Lagrangian subproblem is feasible, via a binary search over the
    /// sorted array of link capacities (see capSort())
    double* Limitrmin();

    /// recursive randomized quicksort of caps[ sx .. dx ]
    void capSort(int sx, int dx);

    /// quicksort partition step around caps[ pivot ], on caps[ sx .. dx ]
    int Distrib(int sx, int pivot, int dx);

    /// swaps caps[ a ] and caps[ b ]
    void Swap(int a, int b);


    /// performs the exact line search over all subintervals of Q, finding
    /// the position in Q of the subinterval whose piecewise-linear cuts'
    /// envelope attains the smallest minimum, i.e., the next candidate
    /// optimum; also updates BestLB and fathoms (branchedflag) subintervals
    /// that cannot improve on BestUB
    int LineSearch();//returns the position in Q of the optimal intersection between the cuts

    /// among the cuts of subinterval Q[ i ], returns the position of the
    /// one attaining the highest value at abscissa Q[ i - 1 ].rmin
    int whchbest(int i); // support method for the LS.

    /// updates the nonnegative-slope cut pCut or the negative-slope cut
    /// mCut of subinterval Q[ i ] with the line of slope beta and
    /// intercept alpha, if it differs from both of them
    void UpdCut(double alpha, double beta, int i); // support method for the LS.


    /// deep-copies flow/links/nodes data into the Flow/Links/Nodes members
    void copyDataArray(DCR::DCRFlow flow, vector<DCR::DCRLink> links, vector<DCR::DCRNode> nodes);

    /// releases all dynamically allocated memory (Links, Nodes, caps, Q, ...)
    void clean_up();




};

#endif

/*--------------------------------------------------------------------------*/
/*-------------------------- End File BenBound.h ----------------------------*/
/*--------------------------------------------------------------------------*/
