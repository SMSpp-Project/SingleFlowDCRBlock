/*--------------------------------------------------------------------------*/
/*--------------------------- File BenBound.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the BenBound class.
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
 * \copyright &copy; by Antonio Frangioni, Laura Galli, Luca Mencarelli,
 *                      Enrico Sorbera
 */

/*--------------------------------------------------------------------------*/
/*--------------------------- IMPLEMENTATION -------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "BenBound.h"
#include "DCRLagrangianSolver.h"
#include "DCR.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <ctime>
#include <cstdlib>

#include <unistd.h>


/*--------------------------------------------------------------------------*/
/*--------------------- IMPLEMENTATION OF BenBound--------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------------- PUBLIC METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*----------------------------- CONSTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

// gives default ("empty") values to all data members; in particular no
// instance is loaded yet (Links/Nodes/caps are null) and BenStat == Error,
// so LoadProblem() must be called before Solve() can do anything useful

BenBound::BenBound()
{
 eps = 1e-6;
 BenStat = Error;

 timer = 0;

 Links = 0;
 Nodes = 0;
 caps = 0;
 solvedflag = 0;

 //FIXME: create a public method for the upd (before it was set to 0.995)
 myparam = 0.995; //for garr and sndlib 0.995
 //myparam = 1-1e-6;
 //myparam = 0.9995;

 SPLabels.resize( 1 );
 Q.resize( 1 );
 }

/*--------------------------------------------------------------------------*/
/*----------------------------------LOAD DATA ------------------------------*/
/*--------------------------------------------------------------------------*/

// loads a brand new instance: discards any previous one, deep-copies the
// topology/link/node/flow data and (re)initializes all bounds and
// counters to their "nothing solved yet" values; also (re)loads the same
// data into the Lagrangian subproblem solver lagSol

// void BenBound::LoadProblem (int nnodes, int nlinks, DCR::DCRFlow flow,
// DCR::DCRLink *links, DCR::DCRNode *nodes, double mtu)
void BenBound::LoadProblem( int nnodes , int nlinks , DCR::DCRFlow flow ,
                            std::vector< DCR::DCRLink > links ,
                            std::vector< DCR::DCRNode > nodes , double mtu )
{
 clean_up();

 numNodes = nnodes;
 numLinks = nlinks;
 MTU = mtu;
 ObjVal = Inf< double >();
 HeurVal = Inf< double >();
 BestUB = Inf< double >();
 BestUBrmin = -1;
 BestLB = -Inf< double >();
 ApproxVal = -Inf< double >();
 solvedflag = 0;
 counter_ite_Ben = 0;
 counter_ite_Lag = 0;
 is_convex_iteration = false;

 SPLabels.resize( numNodes );

 copyDataArray( flow , links , nodes );

 SOLUTION.resize( numLinks );
 HeurSOLUTION.assign( numLinks , 0.0 );

 lagSol.LoadProblem( numNodes , numLinks , Flow , Links , Nodes , MTU , 1 );
 }

/*--------------------------------------------------------------------------*/

// reloads only the flow (and the arc costs it implies) keeping the same
// topology, clearing the cutting-plane data structures (Q) and resetting
// the bounds/counters, then reloads lagSol with the new flow

void BenBound::LoadProblem( DCR::DCRFlow flow )
{
 int i;

 for( int i = 0 ; i < Q.size() ; i++ )
  Q[ i ].Cuts.clear();
 Q.clear();

 ObjVal = Inf< double >();
 HeurVal = Inf< double >();
 HeurSOLUTION.assign( numLinks , 0.0 );
 BestUB = Inf< double >();
 BestUBrmin = -1;
 BestLB = -Inf< double >();
 ApproxVal = -Inf< double >();
 solvedflag = 0;
 counter_ite_Ben = 0;
 counter_ite_Lag = 0;
 is_convex_iteration = false;

 Flow = flow;

 for( i = 0 ; i < numLinks ; i++ )
  Links[ i ].cost = Flow.costs[ i ];

 lagSol.LoadProblem( numNodes , numLinks , Flow , Links , Nodes , MTU , 1 );
 }

//method for loading, with the same topology, only the new flow to be served.

/*--------------------------------------------------------------------------*/
/*--------------------------- GET METHODS-----------------------------------*/
/*--------------------------------------------------------------------------*/

/// returns the value of d( . ) at the last evaluated r_min
double BenBound::getObjVal()
{
 //return(SOL_VALUE);
 return( ObjVal );
 }


/*--------------------------------------------------------------------------*/

/// returns r_min at which the (best) optimum was found
double BenBound::getr_min()
{
 return( rmin );
 }


/*--------------------------------------------------------------------------*/

/// best available upper bound: the smaller of BestUB and HeurVal

double BenBound::getUB()
{
 //std::cout << "UB=" << BestUB << std::endl;
 //return(ObjVal);
 return( std::min( BestUB , HeurVal ) );
 }

/*--------------------------------------------------------------------------*/

// best available lower bound: the master-problem value BestLB, clipped
// from above by getUB() (it can never exceed a known upper bound) and
// from below by 0 (the DCR objective cannot be negative)

double BenBound::getLB()
{
 //std::cout << "LB=" << BestLB << std::endl;
 //std::cout << HeurVal << "," << ObjVal << "," << ApproxVal << std::endl;
 //std::cout<<"BenStat="<<BenStat<<std::endl;
 //if(BestUB < BestLB || BestLB < -1e-6)
 //return(min(BestUB,HeurVal));
 return( std::max( 0.0 , std::min( std::min( BestUB , HeurVal ) , BestLB ) ) );
 //return(HeurVal);
 }


/*--------------------------------------------------------------------------*/


/// returns the overall status of the bounding procedure

BenBound::BndrStat BenBound::getStat()
{
 return( BenStat );
 }

/*--------------------------------------------------------------------------*/

/// returns the value of the best primal-feasible solution found so far

double BenBound::getHeurVal()
{
 return( HeurVal );
 }

/*--------------------------------------------------------------------------*/

/// returns the reserved rate for link i in the best solution found

double BenBound::getSolution( int i )
{
 return( SOLUTION[ i ] );
 }

/*--------------------------------------------------------------------------*/

/// returns the reserved rate for link i in the best verified
/// delay-feasible solution found (see getHeurVal())

double BenBound::getHeurSolution( int i )
{
 return( HeurSOLUTION.empty() ? 0.0 : HeurSOLUTION[ i ] );
 }

/*--------------------------------------------------------------------------*/

/// returns the current cutting-plane (master problem) approximate optimum

double BenBound::getApproxVal()
{
 return( ApproxVal );
 }

/*--------------------------------------------------------------------------*/

/// returns the number of Benders (outer) iterations performed by Solve()

int BenBound::getNumIterationBender()
{
 return( counter_ite_Ben );
 }

/*--------------------------------------------------------------------------*/

/// returns the total number of Lagrangian (inner) iterations performed

int BenBound::getNumIterationLagr()
{
 return( counter_ite_Lag );
 }

/*--------------------------------------------------------------------------*/

/// returns true if at least one (locally convex) subinterval was found

bool BenBound::IsConvexInteration()
{
 return( is_convex_iteration );
 }

/*--------------------------------------------------------------------------*/

// debug/plotting helper: (re)initializes the search via Inizial() and
// returns numPoints == 1000 abscissae equally spaced over the first
// subinterval Q[ 1 ], from its optimal point Q[ 1 ].inter to its right
// endpoint Q[ 1 ].rmin; the result is cached in xPlot for getyPlotData()

std::vector< double > BenBound::getxPlotData()
{
 int i;
 numPoints = 1000;

 xPlot.resize( numPoints );

 Inizial();

 for( i = 1 ; i < numPoints + 1 ; i++ ) {
  xPlot[ i - 1 ] =
   ( i * ( Q[ 1 ].rmin - Q[ 1 ].inter ) / numPoints ) +
   Q[ 1 ].inter; //store the subdivision of the interval into 100 points
  }

 return( xPlot );
 }


/*--------------------------------------------------------------------------*/

// debug/plotting helper: (re)solves the Lagrangian subproblem at each of
// the abscissae in xPlot (see getxPlotData(), which must be called
// first) and returns the corresponding values of the dual function
// d( r_min ) = lambda * burst / r_min - lambda * deadline + SPLabels[ t ]

std::vector< double > BenBound::getyPlotData()
{
 int i;
 double lambda;
 std::vector< double > yPlot;
 std::vector< int >
  RedGraph; //reduced graph that we will have the Lagrangian solver pass to us
 int cardRedGraph;

 yPlot.resize( numPoints );

 for( i = 0 ; i < numPoints ; i++ ) {
  lagSol.updrmin( xPlot[ i ] );

  lagSol.Solve();

  lambda = lagSol.getLambda();
  SPLabels = lagSol.getSPLabels();
  cardRedGraph = lagSol.getCardRedGraph();
  RedGraph.resize( cardRedGraph );
  RedGraph = lagSol.getRedGraphPos();

  yPlot[ i ] =
   lambda * Flow.burst / xPlot[ i ] - lambda * Flow.deadline +
   SPLabels[ Flow.sinknode ]; //value of the objective function at xplot[i]
  }

 return yPlot;
 }


/*--------------------------------------------------------------------------*/

/// returns the elapsed time, in seconds, or 0 if no timer was set

double BenBound::getTime( void ) const
{
 return( timer ? timer->Read() : 0 );
 }


/*--------------------------------------------------------------------------*/
/*------------------------- SET TIME METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

void BenBound::DCRsetTime( bool timeON )
{
 if( timeON )
  if( timer )
   timer->ReSet();
  else
   timer = new DCRtimer();
 else
  delete timer;
 }

/**< If timeON is true sets or resets the timer, if false deletes the timer.
 * \param timeON bool value */

/*--------------------------------------------------------------------------*/

void BenBound::DCRstartTime( void )
{
 if( ! timer )
  throw( std::logic_error( "BenBound::DCRstartTime: no timer was set" ) );

 timer->Start();
 }

/**< If timer was set or re-set, starts or re-starts timer ticking. */


/*--------------------------------------------------------------------------*/

void BenBound::DCRstopTime( void )
{
 if( ! timer )
  throw( std::logic_error( "BenBound::DCRstopTime: no timer was set" ) );

 timer->Stop();
 }

/**< If timer was set, stops timer ticking. */

/*--------------------------------------------------------------------------*/

void BenBound::DCRsetTimeLimit( long secs )
{
 tlimit = secs;
 }

/**< Sets a time limit for the solver.
 * \param secs long expressing timelimit in seconds */

/*--------------------------------------------------------------------------*/
/*----------------------------------SOLVE-----------------------------------*/
/*--------------------------------------------------------------------------*/

/// runs the Benders/Kelley cutting-plane bounding procedure
/** This is the core of the class: it repeatedly
 *
 * -# picks the current candidate optimum Q[ p ].inter (initialized by
 *    Inizial(), then chosen at each iteration by LineSearch());
 *
 * -# (re)solves the Lagrangian subproblem there (lagSol.Solve()), unless
 *    it had already been solved for that very point (Q[ p ].solflag);
 *    if the subproblem happens to be infeasible at that point (only
 *    possible while the left endpoint of the very first subinterval has
 *    not yet been made feasible, Q[ p ].infeasflag), the point is first
 *    moved towards feasibility by a geometric search controlled by
 *    mystep;
 *
 * -# out of the optimal Lagrangian multiplier lambda and dual solution
 *    (SPLabels, the shortest-path potentials), builds two new supporting
 *    lines ("cuts") of the dual function d( r_min ) at Q[ p ].inter: the
 *    left cut lc (obtained by summing, arc by arc, the contribution of
 *    each arc whose capacity exceeds the subinterval's left endpoint,
 *    following the closed-form Lagrangian analysis of the DCR problem)
 *    and the right cut rc (the actual tangent/derivative of the g( . )
 *    burst term at Q[ p ].inter, always available in closed form). If
 *    lc and rc coincide, the dual function is locally *convex* around
 *    Q[ p ].inter and no new breakpoint needs to be introduced (only the
 *    cuts of all existing subintervals are updated); otherwise it is
 *    (possibly) *nonconvex* and a new subinterval is spliced into Q at
 *    position p, splitting the old one, with lc added to the cuts of
 *    all subintervals to its left and rc to all subintervals to its
 *    right;
 *
 * -# updates the running primal bound BestUB with the newly evaluated
 *    d( Q[ p ].inter ) (== Q[ p ].Val);
 *
 * -# calls LineSearch() to find, over the (possibly updated) partition
 *    Q, the position p of the subinterval whose piecewise-linear cuts'
 *    envelope has the smallest minimum: this both gives the new value of
 *    BestLB (the master-problem/dual bound) and selects the next trial
 *    point Q[ p ].inter for the following iteration.
 *
 * The loop stops (solvedflag = 1) as soon as BestUB and BestLB coincide
 * up to the relative tolerance eps, or after at most 5000 iterations; in
 * the latter case ObjVal is set to the current cutting-plane
 * approximation "approx" so that a valid (if not provably optimal)
 * bound is always returned. Finally, the primal solution corresponding
 * to the last solved Lagrangian subproblem is retrieved from lagSol via
 * getRSOLS(). */

void BenBound::Solve()
{
 int i , j;
 int Qsize;
 double approx;
 int p = 1;         //insertion position in the new subdivision of Q
 int nonConvexflag; //flag to distinguish the case of convex cuts
 int triedflag =
  0; //flag to tell whether we have already tried the jump to the optimum
 double lambda;
 int noLSneeded = 0;

 DCRLagrangianSolver::LAGStatus lgstat;
 DCRLagrangianSolver::LinearCut lc; //left cut
 DCRLagrangianSolver::LinearCut rc; //right cut

 double HeurApprox;

 double sqr;
 double barl;
 double SPcost;
 double releps;
 double zetam;
 double zetap;
 double zetaunico;
 double d;
 double tildez;
 double discrim; //if <=0 that arc gives a zero contribution
 double ms;      //helper for the computation, partial value of the slope

 Inizial();

 //std::cout << "Inizial()" << std::endl;

 // cout<<"inizializzato"<<endl;
 Qsize = Q.size(); // = 2, limit values.

 //cout<<"limiti: "<<Q[0].rmin<<"-"<<Q[1].rmin<<endl;

 //before next line was commented
 //solvedflag = 0;

 if( solvedflag == 0 ) {
  if( BenStat == 0 ) {
   do {
    // for(int count_feas = Q[0].rmin; count_feas < Q[1].rmin; count_feas +=
    // 1000){

    //           lagSol.updrmin(count_feas);
    //           cout<<" : "<<lagSol.isFeasible()<<endl;
    //           }

    //         cout<<"num of iteration"<< counter_ite_Ben<<endl;
    // cout<<"lagstat = "<<lgstat<<" lambda = " <<lambda<<" Q[p]inter = "<<
    // Q[p].inter<< " value = "<<Q[p].Val<<endl;
    // cout<<"Objval = "<<ObjVal<<" Heurval = "<<HeurVal<<endl;
    // cout<<"LB = "<<BestLB<<" UB = "<< BestUB<<endl;

    double BLB = BestLB;
    double HUB = HeurVal;

    bool is_feas;
    double Qnew;

    lagSol.updrmin(
     Q[ p ].inter ); //let the new candidate show us its potential :)
    lagSol.Solve();
    counter_ite_Lag += lagSol.getNumIte();
    lgstat = lagSol.getStatus();

    if( Q[ p ].solflag == 1 && counter_ite_Ben > 2 ) {
     // the point has been solved already: move it a little, towards the
     // left endpoint of the window of the r_min that can be feasible at
     // all, exactly as the preliminary phase moves it towards the right
     // one [see Inizial()]. Note that the endpoint is where it stops: a
     // multiplicative shrink of the point alone walks past it, and an
     // r_min below Q[ 0 ].rmin reserves less than the sustained rate of
     // the flow, i.e. it is not a solution of the problem
     Q[ p ].inter = myparam * Q[ p ].inter +
                    ( 1 - myparam ) * Q[ 0 ].rmin;

     lagSol.updrmin(
      Q[ p ].inter ); //let the new candidate show us its potential :)
     is_feas = lagSol.isFeasible();

     if( lgstat == 0 )
      Q[ p ].solflag = 0;
     }

    if( lgstat == 1 ) {

     for( int i = 0 ; i < 1000 ; i++ ) {

      Qnew = ( i * ( Q[ 1 ].rmin - Q[ 0 ].rmin ) / 1000 ) + Q[ 0 ].rmin;

      lagSol.updrmin( Qnew );
      lagSol.Solve();
      counter_ite_Lag += lagSol.getNumIte();
      lgstat = lagSol.getStatus();

      if( lgstat == 0 ) {
       break;
       }
      }
     Q[ p ].inter = Qnew;
     }


    if( Q[ p ].solflag == 0 ) {
     if(
      Q[ p ].infeasflag ==
      1 ) //fix things up for the Lagrangian by finding a new feasibility point
     {
      //cout<<"entering infeasflag ==1"<<endl;
      double
       temppoint;  // shrink the interval, using the last infeasible value as
                   // the lower bound
      double oldinterx;
      int feas = 1;

      oldinterx = Q[ 1 ].inter;

      //look for the other feasibility point
      int count_iter_feas = 0;

      while( feas != 0 && count_iter_feas <= 20 ) {
       temppoint =
        mystep * Q[ 0 ].rmin +
        ( 1 - mystep ) *
         Q[ 1 ].inter; //let's see which one changes based on feasibility

       lagSol.updrmin( temppoint );
       feas = lagSol.isFeasible();

       if( feas != 0 )
        Q[ 0 ].rmin = temppoint;
       else
        Q[ 1 ].inter = temppoint;
       count_iter_feas++;
       }

      Q[ 1 ].interVal =
       Q[ 1 ].Cuts[ Q[ 1 ].bestCutpos ].q +
       Q[ 1 ].inter *
        Q[ 1 ]
         .Cuts[ Q[ 1 ].bestCutpos ]
         .m;  // FIXME: actually it is not certain that bestCutpos is the
              // optimal position
      Q[ 1 ].solflag = 0;
      Q[ 1 ].infeasflag = 0;

      //if I no longer move, I settle for this.

      if( Q[ 1 ].inter < 1 )
       releps = 1;
      else
       releps = Q[ 1 ].inter;

      //before next two lines were commented
      //if(abs(Q[1].inter - oldinterx) < eps * releps/100)
      //   noLSneeded = 1;

      } //if infeasflag = 1

     lagSol.updrmin(
      Q[ p ].inter ); //let the candidate show us its potential :)
     lagSol.Solve();
     counter_ite_Lag += lagSol.getNumIte();

     lgstat = lagSol.getStatus();
     lambda = lagSol.getLambda();

     HeurApprox = lagSol.getHeurVal();

     if( HeurApprox <
         HeurVal ) { //if we found a better feasible solution, we save it
      HeurVal = HeurApprox;
      // lagSol's own HeurVal (and hence the routing behind it) is reset
      // on every updrmin() call, so getHeurRSOLS() here is guaranteed to
      // be the routing that attains *this* HeurApprox, from the solve
      // just performed at Q[ p ].inter above -- not a stale one from a
      // different candidate rmin
      HeurSOLUTION = lagSol.getHeurRSOLS();
      }

     SPLabels = lagSol.getSPLabels();

     if(
      lgstat ==
      0 )  // this should be guaranteed by limitrmin and by the check just
           // above the re-optimization
     {

      lc.m = 0;
      lc.q = 0;

      ///////COMPUTATION OF NEW CUTS
      // compute the two new cuts (supporting lines) of d( . ) at
      // Q[p].inter, out of the just-computed optimal multiplier
      // lambda and shortest-path potentials SPLabels: the closed-
      // form Lagrangian relaxation of the DCR problem decomposes
      // additively over the arcs, so both the slope (lc.m) and
      // the intercept (lc.q) of the left cut are accumulated arc
      // by arc below (only arcs whose capacity exceeds the left
      // endpoint Q[0].rmin can contribute; each arc's closed-form
      // contribution depends on where Q[p].inter falls relative
      // to the arc's own capacity and to the breakpoints
      // zetam/zetap/zetaunico of its (piecewise) cost function)

      ///Left cut

      // computation of the value taking into account the contribution of each
      // arc
      // the intercepts are also summed arc by arc, since some of them cause
      // discontinuities.

      sqr = sqrt( lambda * MTU / Links[ 1 ].cost );

      for( i = 0 ; i < numLinks ; i++ ) {
       if( Links[ i ].capacity > Q[ 0 ].rmin ) {
        sqr =
         sqrt( lambda * MTU / Links[ i ].cost ); //sqrt(lambda * L / f_ij)

        //set to 0 the SP duals that are at +\infty
        if( SPLabels[ Links[ i ].endnode ] > 1e200 )
         SPLabels[ Links[ i ].endnode ] = 0;
        if( SPLabels[ Links[ i ].startnode ] > 1e200 )
         SPLabels[ Links[ i ].startnode ] = 0;

        d = -SPLabels[ Links[ i ].endnode ] +
            SPLabels[ Links[ i ].startnode ]; //d_i-d_j

        if( d < -1e-10 ) //otherwise the contribution is certainly zero
        {
         if( lambda == 0 ) //first, simplified case
         {
          zetaunico = ( -d ) / Links[ i ].cost; //(d_j - d_i)/f_ij
          if( zetaunico - Q[ 0 ].rmin > zetaunico * eps &&
              zetaunico - Links[ i ].capacity < eps * Links[ i ].capacity &&
              zetaunico - Q[ p ].inter < eps * zetaunico ) {
           ms = ( Links[ i ].cost * Q[ 0 ].rmin + d ) /
                ( Q[ 0 ].rmin - Q[ p ].inter );
           lc.m += ms;
           lc.q += -ms * Q[ p ].inter;
           }

          if( Links[ i ].capacity - Q[ 0 ].rmin > Links[ i ].capacity * eps &&
              zetaunico - Links[ i ].capacity > eps * zetaunico &&
              Q[ p ].inter - Links[ i ].capacity > Q[ p ].inter * eps ) {
           if( zetaunico - Q[ p ].inter > zetaunico * eps ) {
            ms = ( Links[ i ].cost * Links[ i ].capacity + d ) /
                 ( Links[ i ].capacity - Q[ p ].inter );
            lc.m += ms;
            lc.q += -ms * Q[ p ].inter;
            }
           if( zetaunico - Q[ p ].inter < Q[ p ].inter * eps ) {
            ms = ( Links[ i ].cost * Q[ 0 ].rmin + d ) /
                 ( Q[ 0 ].rmin - Q[ p ].inter );
            lc.m += ms;
            lc.q += -ms * Q[ p ].inter;
            }
           }
          } //lambda = 0
         else //general case
         {
          barl = MTU / Links[ i ].speed + Links[ i ].delay +
                 Nodes[ Links[ i ].startnode ].delay;
          discrim =
           pow( lambda * barl + d , 2 ) - 4 * Links[ i ].cost * lambda * MTU;

          if( discrim > 0 ) //then zeta_\pm are well defined
          {
           zetap = ( -lambda * barl - d + sqrt( discrim ) ) /
                   ( 2 * Links[ i ].cost ); //larger root
           zetam = ( -lambda * barl - d - sqrt( discrim ) ) /
                   ( 2 * Links[ i ].cost ); //smaller root

           if(
            zetam - Links[ i ].capacity < eps * Links[ i ].capacity &&
            Links[ i ].capacity - sqr < eps * Links[ i ].capacity &&
            Q[ p ].inter - Links[ i ].capacity >
             eps *
              Links[ i ]
               .capacity ) {  // && Q[p].inter - Links[i].capacity >
                              // eps*Links[i].capacity
            ms = ( Links[ i ].cost * Links[ i ].capacity +
                   lambda * MTU / Links[ i ].capacity + lambda * barl + d ) /
                 ( Links[ i ].capacity - Q[ p ].inter );
            lc.m += ms;
            // c.q += (Links[i].cost * Links[i].capacity + lambda * MTU /
            // Links[i].capacity + lambda * barl + d) - ms * Links[i].capacity;
            lc.q += -ms * Q[ p ].inter;
            }

           if( sqr - Links[ i ].capacity < eps * sqr &&
               Links[ i ].capacity - zetap < eps * zetap &&
               Q[ p ].inter >=
                Links[ i ].capacity ) { //&& Q[p].inter >= Links[i].capacity
            if( Q[ p ].inter - zetap < Q[ p ].inter * eps ) {
             if( std::abs( sqr - Q[ p ].inter ) > eps * sqr ) {
              ms = ( Links[ i ].cost * sqr + lambda * MTU / sqr +
                     lambda * barl + d ) /
                   ( sqr - Q[ p ].inter );
              if( Links[ i ].cost - ms > eps * Links[ i ].cost &&
                  ms >= 0.1 ) {
               tildez = sqrt( ( lambda * MTU ) / ( Links[ i ].cost - ms ) );
               lc.q += Links[ i ].cost * tildez + lambda * MTU / tildez +
                       lambda * barl + d - tildez * ms;
               lc.m += ms;
               //std::cout << ms << std::endl;
               }
              else {
               lc.q +=
                ( Links[ i ].cost * Links[ i ].capacity +
                  lambda * MTU / Links[ i ].capacity + lambda * barl + d ) -
                ms * Links[ i ].capacity;
               lc.m += ms;
               //std::cout << ms << std::endl;
               }
              }
             }
            else {
             if( std::abs( sqr - Q[ p ].inter ) > eps * sqr ) {
              ms = ( Links[ i ].cost * sqr + lambda * MTU / sqr +
                     lambda * barl + d ) /
                   ( sqr - Q[ p ].inter );
              if( Links[ i ].cost - ms > eps * Links[ i ].cost &&
                  ms >= 0.1 ) {
               ms = ( Links[ i ].cost * sqr + lambda * MTU / sqr +
                      lambda * barl + d ) /
                    ( sqr - Q[ p ].inter );
               lc.m += ms;
               tildez = sqrt( ( lambda * MTU ) / ( Links[ i ].cost - ms ) );
               lc.q += Links[ i ].cost * tildez + lambda * MTU / tildez +
                       lambda * barl + d - tildez * ms;
               //std::cout << ms << std::endl;
               }
              }
             }
            }

           if( Links[ i ].capacity - zetap > eps * Links[ i ].capacity ) {
            if( Links[ i ].capacity - Q[ p ].inter >
                eps * Links[ i ].capacity ) {
             if( std::abs( sqr - Q[ p ].inter ) > eps * sqr ) {
              ms = ( Links[ i ].cost * sqr + lambda * MTU / sqr +
                     lambda * barl + d ) /
                   ( sqr - Q[ p ].inter );
              if( Links[ i ].cost - ms > eps * Links[ i ].cost &&
                  ms >= 0.1 ) { //-10
               lc.m += ms;
               tildez = sqrt( ( lambda * MTU ) / ( Links[ i ].cost - ms ) );
               lc.q += Links[ i ].cost * tildez + lambda * MTU / tildez +
                       lambda * barl + d - tildez * ms;
               //std::cout << ms << std::endl;
               }
              }
             }
            else {
             if( std::abs( sqr - Q[ p ].inter ) > eps * sqr ) {
              ms = ( Links[ i ].cost * sqr + lambda * MTU / sqr +
                     lambda * barl + d ) /
                   ( sqr - Q[ p ].inter );
              if( Links[ i ].cost - ms > eps * Links[ i ].cost &&
                  ms >= 0.1 ) {
               lc.m += ms;
               tildez = sqrt( ( lambda * MTU ) / ( Links[ i ].cost - ms ) );
               lc.q += Links[ i ].cost * tildez + lambda * MTU / tildez +
                       lambda * barl + d - tildez * ms;
               //std::cout << ms << std::endl;
               }
              }
             }
            }
           }  // if the points where the arc's contribution vanishes are well
             // defined
          //otherwise the contribution will be zero.
          } //case lambda > 0
         } //if d < 0
        } //only the arcs with capacity greater than the minimum feasible one!
       } //for all the arcs.

      //finally we have to sum the contributions of the function g(z).
      // g( r_min ) = lambda * burst / r_min - lambda * deadline +
      // SPLabels[ sink ] is the (always concave) burst term of
      // d( . ); its slope ms is added to the left cut lc, while
      // the right cut rc below is exactly its tangent at
      // Q[p].inter (a valid supporting line for the whole d( . )
      // on the right of Q[p].inter)

      ms = -lambda * Flow.burst / pow( Q[ p ].inter , 2 );
      lc.m += ms;
      lc.q += lambda * Flow.burst / Q[ p ].inter - lambda * Flow.deadline +
              SPLabels[ Flow.sinknode ] -
              Q[ p ].inter * ms; //g(r_min) - r_min * ms

      //right cut

      rc.m = -lambda * Flow.burst / pow( Q[ p ].inter , 2 );
      rc.q = lambda * Flow.burst / Q[ p ].inter - lambda * Flow.deadline +
             SPLabels[ Flow.sinknode ] - Q[ p ].inter * rc.m;

      //cout<< " the cuts produced: lc.m = "<<lc.m<<" lc.q = "<<lc.q<<endl;
      //cout<<" rc.m = "<<rc.m<<" rc.q = "<<rc.q<<endl;

      Q[ p ].Val = rc.q + Q[ p ].inter * rc.m; //take the value at the point.
      Q[ p ].solflag = 1;

      ObjVal = Q[ p ].Val;

      if( Q[ p ].Val < BestUB ) {
       BestUB = Q[ p ].Val;
       BestUBrmin = Q[ p ].inter;
       }

      /*<convexity holds when lc.m <= rc.m; now rc.m is negative and lc.m
        is obtained by summing positive terms, so in fact the only case in
        which convexity can hold is when they are equal.*/
      // NB: despite its name, nonConvexflag == 1 here signals the
      // (locally) *convex* case (lc and rc coincide): see the
      // "CONVEX CASE" branch below, which is taken precisely
      // when nonConvexflag != 0

      if( lc.m == rc.m && lc.q == rc.q )
       nonConvexflag = 1; //FIXME: equality up to rounding errors
      else
       nonConvexflag = 0;

      /////////NON-CONVEX CASE
      // a genuinely new breakpoint is needed: split Q[p] into two
      // subintervals at Q[p].inter (the new one, I, inherits the
      // "left" part), and propagate lc/rc as new candidate cuts
      // to, respectively, every subinterval to the left/right of
      // the split point

      if( nonConvexflag ==
          0 ) //add the new point to Q and the cuts it produces
      {

       BenBound::SubInterval I;

       I.Cuts.resize( Q[ p ].Cuts.size() );

       for( i = 0 ; i < Q[ p ].Cuts.size() ; i++ ) {
        I.Cuts[ i ].m =
         Q[ p ]
          .Cuts[ i ]
          .m; //initialize with the cuts of the subinterval it was contained in
        I.Cuts[ i ].q = Q[ p ].Cuts[ i ].q;
        }

       I.rmin = Q[ p ].inter;
       I.Val = Q[ p ].Val;
       I.inter = Q[ p ].inter; //for now, the solution is exactly this!
       I.interVal = Q[ p ].interVal;
       I.pCut.m = Q[ p ].pCut.m; //we also copy the cuts that define it
       I.pCut.q = Q[ p ].pCut.q; //from which the LS will start.
       I.mCut.m = Q[ p ].mCut.m;
       I.mCut.q = Q[ p ].mCut.q;
       I.solflag = 1;
       I.branchedflag = 0;
       I.infeasflag = 0;

       if( Q[ p ].Cuts.size() == 1 ) {
        Q[ p ].mCut.q = rc.q;
        Q[ p ].mCut.m = rc.m;

        I.mCut.q = lc.q;
        I.mCut.m = lc.m;
        }  // if we were in the case defined by a single cut, now we also have
          // the other one

       Q.insert( Q.begin() + p , I );
       Qsize++;

       //add lc to all points on the left:

       for(
        i = 1; i <= p;
        i++ ) //not at pos 0, since it does not really represent a subinterval
       {
        Q[ i ].Cuts.push_back( lc );
         }

       //add rc to all points on the right:

       for( i = p + 1 ; i < Qsize ; i++ ) {
        Q[ i ].Cuts.push_back( rc );
        }

       } //if non convex.

      ///  ///CONVEX CASE: no new point is added to Q, only the cuts
      ///  ///produced at
      /// each pre-existing position
      else {
       /*cout<<"convex case"<<endl;*/
       is_convex_iteration = true;

       if( Q[ p ].Cuts.size() == 1 ) {
        Q[ p ].mCut.q = rc.q;
        Q[ p ].mCut.m = rc.m;
        }  // if we were in the case defined by a single cut, now we also have
          // the other one

       for( i = 1 ; i < Qsize ; i++ ) {
        Q[ i ].Cuts.push_back(
         lc ); //also at position p, since it had not been added before
        }
       } // if convex
      } //if the Lagrangian finds a solution
     else //otherwise it means we are getting something wrong
     {
      throw( std::logic_error( "BenBound::Solve: the Lagrangian "
                               "relaxation found no solution" ) );
      }
     } //if we had not already solved for that value

    //next line originally not in the code
    //if(Q[p].solflag == 1){noLSneeded = 1;}

    if( noLSneeded == 0 ) //if(noLSneeded == 0)
    {
     p = LineSearch();
     //if(Q.size()==2||Q.size()==3) cout<<"P = "<<p<<endl;
     }

    if( Q[ p ].interVal > BestLB )
     BestLB = Q[ p ].interVal;

    ObjVal = Q[ p ].Val;

    approx =
     Q[ p ].interVal; //improve our estimates with the values given by the LS

    //std::cout << "OV = " << ObjVal << std::endl;

    if( approx > 1 )
     releps = approx; //relative precision
    else
     releps = 1;

    // stopping criteria: either this iteration made no progress at all on
    // both bounds (HUB/BLB are the bounds at loop entry), or the
    // relative gap between the best upper and lower bound has closed to
    // within eps.
    //
    // The upper bound used here is HeurVal, *not* BestUB: BestUB is
    // Q[ p ].Val, a Lagrangian-relaxation value alpha + lambda * beta for
    // whatever candidate attained it, which equals the candidate's actual
    // (achievable) cost alpha only when beta == 0 there. Whenever the line
    // search lands on a "convex case" point (a single, exact supporting
    // line, no further subdivision needed), Q[ p ].interVal trivially
    // equals Q[ p ].Val in that very same iteration, so BestUB and BestLB
    // coincide immediately regardless of whether lambda is anywhere near
    // the true dual-maximizing multiplier for that r_min -- e.g. a
    // candidate with beta strictly negative (feasible, but with slack, so
    // its lambda has NOT yet been pushed down to the point of exact
    // complementary slackness) can trigger this "gap closed" check while
    // its real cost is nowhere near BestLB. HeurVal cannot suffer from
    // this: it is only ever updated from a candidate DCRLagrangianSolver
    // itself verified delay-feasible (beta < 0), so it is always alpha,
    // the genuine achievable cost, never a dual value in disguise.
    //
    // Both criteria below also double as safety valves bounding how far
    // LineSearch() can push r_min: without them it can keep proposing
    // candidates past the point where the Lagrangian subproblem is still
    // well-posed, reaching an r_min the reduced graph cannot support at
    // all and running straight into the "no solution" logic_error thrown
    // further down. So, imprecise as they are as *convergence* signals,
    // they are load-bearing as *bracket* ones and must stay
    if( HeurVal == HUB && BestLB == BLB ) {
     solvedflag = 1;
     }

    if( ( HeurVal - BestLB ) / HeurVal < eps ) {
     solvedflag = 1;
     }

    counter_ite_Ben++; //number of points visited

    // if(Q.size()==2||Q.size()==3) cout<<"number of points = "<<Q.size()<<"
    // with rmin = "<<rmin<<endl;

    // }while(solvedflag == 0 && abs(ObjVal - approx) > eps * releps); (before
    // on the code)
    } while( solvedflag == 0 &&
            ( std::abs( ObjVal - approx ) > eps * releps / 100 ) &&
            counter_ite_Ben < 5000 );
   //std::cout << "HV = " << HeurVal << std::endl;
   //}while(counter_ite_Ben<100 && abs(BestUB-BestLB)>eps*BestUB);

   if( solvedflag == 0 )
    ObjVal = approx; //in any case we must return an LB!

   rmin = Q[ p ].inter;

   // Q[ p ] here is simply wherever the last LineSearch() call landed,
   // which is *not* necessarily the point that attained BestUB: the
   // stopping criteria above can (and on some instances do) fire while
   // the line search is sitting on a different, worse-or-infeasible-
   // looking subinterval than the one that produced the best primal
   // value seen so far. Reporting rmin (and hence extracting SOLUTION
   // below) from the wrong point makes get_ub()/has_var_solution() see
   // an inconsistent, possibly DCR-infeasible routing even though a
   // genuinely feasible one at value BestUB was found earlier in this
   // very call: use that point instead, whenever one was recorded
   if( BestUBrmin > 0 )
    rmin = BestUBrmin;
   //cout << "Q[p].Val = " << Q[p].Val << endl;
   //cout << "rmin = " << rmin << endl;

   } //if there are no errors
  } //if the problem was not already solved

 //SOL_VALUE = ObjVal;//lagSol.getObjVal();
 //std::cout << "BUB=" << BestUB << " BLB=" << BestLB << std::endl;
 //std::cout << "GAP=" << std::abs(getUB() - getLB())/getUB() << std::endl;

 // the line search above can leave lagSol's internal state (path[],
 // Linksp[].rstar, ...) corresponding to an rmin different from the one
 // just reported: the "already solved, nudge the point and just check
 // isFeasible()" shortcut a few lines above updates lagSol's rmin and
 // probes feasibility without a full Solve(), so if that happens to be
 // the last thing done before this point is accepted, lagSol's solution
 // reflects Inizial()'s quick initial trial rather than a converged
 // shortest path at the *reported* rmin. Force one final, full solve at
 // the exact value being reported so the primal solution extracted below
 // is always internally consistent with it (updrmin() is a cheap no-op,
 // and Solve() likewise, when lagSol is already solved for this rmin, so
 // this only does real work in the mismatched case)
 lagSol.updrmin( rmin );
 lagSol.Solve();

 SOLUTION = lagSol.getRSOLS();
 }

/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

// releases all dynamically allocated memory

BenBound::~BenBound()
{
 clean_up();
 }

/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

// deep-copies the link/node/flow data passed to LoadProblem() into the
// Links/Nodes/Flow data members (numLinks/numNodes must already be set)

// void BenBound::copyDataArray(DCR::DCRFlow flow, DCR::DCRLink *links,
// DCR::DCRNode *nodes)
void BenBound::copyDataArray( DCR::DCRFlow flow ,
                              std::vector< DCR::DCRLink > links ,
                              std::vector< DCR::DCRNode > nodes )
{

 int i;
 Links = new DCR::DCRLink[ numLinks ];
 Nodes = new DCR::DCRNode[ numNodes ];

 Flow = flow;

 for( i = 0 ; i < numLinks ; i++ ) {
  Links[ i ].startnode = links[ i ].startnode;
  Links[ i ].endnode = links[ i ].endnode;
  Links[ i ].speed = links[ i ].speed;
  Links[ i ].capacity = links[ i ].capacity;
  Links[ i ].delay = links[ i ].delay;
  Links[ i ].cost = links[ i ].cost; //Flow.costs[i];
  }

 for( i = 0 ; i < numNodes ; i++ ) {
  Nodes[ i ].delay = nodes[ i ].delay;
  }
 }


/*--------------------------------------------------------------------------*/

/// exact line search over all subintervals of Q
/** For each subinterval Q[ i ] (i >= 1) that has not been fathomed
 * (branchedflag == 0), finds the point interx in [ Q[ i - 1 ].rmin ,
 * Q[ i ].rmin ] where the upper envelope of the two "running" cuts
 * pCut/mCut (updated via UpdCut() as the search progresses) meets the
 * highest of all the cuts currently stored in Q[ i ].Cuts; this is the
 * exact minimizer, restricted to the subinterval, of the piecewise-
 * linear model of the (locally convex, once restricted to a single
 * subinterval) dual function. If the intersection falls outside the
 * subinterval, the corresponding endpoint is used instead. Special-cased
 * is the leftmost subinterval (i == 1) when all its cuts still have
 * positive slope: in that case no interior minimizer exists yet
 * (infeasflag is set) and the subinterval's candidate value is simply
 * its left endpoint.
 *
 * Once every subinterval has been processed, the one with the smallest
 * interVal is the best current candidate for the global optimum: its
 * value updates BestLB (a valid lower bound, since interVal is the value
 * of a piecewise-linear *under*-estimator... actually an upper envelope
 * of supporting lines, hence itself a lower bound on d( . )) and
 * ApproxVal, and its position in Q is returned so that Solve() knows at
 * which r_min to (re)solve the Lagrangian subproblem next. While
 * scanning, subintervals whose interVal already exceeds the best known
 * upper bound BestUB are fathomed (branchedflag = 1) since they cannot
 * possibly contain a better solution. */

int BenBound::LineSearch()
{
 int i , j;
 double interx;
 double interxApprox;
 double interxVal;
 double someinter;
 int maxpos;
 double
  max; //in reality this will then be the minimum of the maximum function
 double min;
 int minpos;
 double releps;
 int counter = 0;

 int corrflag = 0;
 int minconflag = 0;
 int isize;
 int branchcounter = 0;

 for(
  i = 1; i < Q.size();
  i++ ) //for every subinterval except the one at position 0, we perform an LS
 {
  //poscounter = 0;
  isize = Q[ i ].Cuts.size();
  max = -Inf< double >(); //reinitialize it to be safe.
  // maxpos must be reset here too: unlike max, it used to carry over from
  // the previous subinterval, whose Cuts vector can have a different size,
  // so a degenerate search below (e.g. a NaN intersection making every
  // comparison false) left Q[ i ].Cuts[ maxpos ] reading out of bounds
  maxpos = 0;
  counter = 0;

  if( Q[ i ].branchedflag == 0 ) {
   // let's run a few checks to see whether we are in the case with all lines
   // positive:
   if( Q[ i ].pCut.m > 0 && Q[ i ].mCut.m > 0 &&
       Q[ i ].Cuts[ isize - 1 ].m > 0 ) //then they are all positive
   {
    maxpos = whchbest( i );
    interx = Q[ i - 1 ].rmin;
    interxVal = Q[ i ].Cuts[ maxpos ].q + interx * Q[ i ].Cuts[ maxpos ].m;
    if( std::abs( Q[ i - 1 ].inter - Q[ i - 1 ].rmin ) >
        0 ) //if it was not its optimal value, it gets reinitialized
    {
     Q[ i ].solflag = 0;
     Q[ i ].Val = Inf< double >();
     }

    if( i == 1 ) //otherwise we need to keep its value
    {

     Q[ i ].bestCutpos = maxpos;
     Q[ i ].infeasflag = 1;
     Q[ i ].interVal = Q[ i ].pCut.q + Q[ i - 1 ].rmin * Q[ i ].pCut.m;
     }
    else {
     Q[ i ].inter = Q[ i - 1 ].rmin;
     Q[ i ].interVal = Q[ i ].pCut.q + Q[ i - 1 ].rmin * Q[ i ].pCut.m;
     }
    }
   else {
    if( Q[ i ].pCut.m > 0 && Q[ i ].mCut.m > 0 &&
        Q[ i ].Cuts[ isize - 1 ].m <= 0 )
     UpdCut( Q[ i ].Cuts[ isize - 1 ].q , Q[ i ].Cuts[ isize - 1 ].m ,
             i ); //if the last one is negative we update

    do //here comes the LS at last!  //then in any case we start the LS
    {
     //poscounter = 0;
     interx =
      ( Q[ i ].pCut.q - Q[ i ].mCut.q ) /
      ( Q[ i ].mCut.m - Q[ i ].pCut.m ); //compute the new intersection

     interxApprox =
      Q[ i ].mCut.q + interx * Q[ i ].mCut.m; //and the height on the old cuts
     //double test = Q[i].pCut.q + interx * Q[i].pCut.m;

     for( j = 0 ; j < isize ;
          j++ ) //look for the cut that gives the maximum height at that point
     {
      someinter = Q[ i ].Cuts[ j ].q + interx * Q[ i ].Cuts[ j ].m;

      if( max <= someinter ) {
       max = someinter;
       maxpos = j;
       }
           }

     interxVal =
      Q[ i ].Cuts[ maxpos ].q +
      interx *
       Q[ i ].Cuts[ maxpos ].m; //finding the true value of the approximation

     UpdCut( Q[ i ].Cuts[ maxpos ].q , Q[ i ].Cuts[ maxpos ].m ,
             i ); //so we update one of the cuts that defines the solution

     counter++;

     if( std::abs( interxVal ) > 1 )
      releps = std::abs( interxVal );
     else
      releps = 1;

     //}while(abs(interxVal - interxApprox) > 1e-3);
     } while(
     std::abs( interxVal - interxApprox ) > eps * releps / 100 &&
     counter <=
      100 );  // since these are exact LS anyway, a releps could optionally be
              // added

    //at this point we must check whether we stayed inside the subinterval,
    //otherwise we will take the closest endpoint as the value.

    if( interx <= Q[ i ].rmin &&
        interx >
         Q[ i - 1 ].rmin ) //if we are inside, the optimum is the one found
    {
     if( std::abs( Q[ i ].inter - interx ) >
         eps * interx / 100 ) //new value, reinitialize everything.
     {
      Q[ i ].solflag = 0;
      Q[ i ].Val = Inf< double >();
      }

     Q[ i ].inter = interx;
     Q[ i ].interVal = interxVal;

     }
    else {
     if( interx <= Q[ i - 1 ].rmin ) //if I am in the previous interval
     {
      if( std::abs( Q[ i - 1 ].inter - Q[ i - 1 ].rmin ) >
          eps * Q[ i - 1 ].rmin /
           100 ) //if it was not its optimal value, it gets reinitialized
      {
       Q[ i ].solflag = 0;
       Q[ i ].Val = Inf< double >();
       }

      if( i == 1 ) //otherwise we need to keep its value
      {
       Q[ i ].bestCutpos = maxpos;
       Q[ i ].infeasflag = 1;
       Q[ i ].interVal = Q[ i ].pCut.q + Q[ i - 1 ].rmin * Q[ i ].pCut.m;
       }
      else {
       Q[ i ].inter = Q[ i - 1 ].rmin;
       Q[ i ].interVal = Q[ i ].pCut.q + Q[ i - 1 ].rmin * Q[ i ].pCut.m;
       }
      }

     else //if I am in the following interval
     {
      if( std::abs( Q[ i ].inter - Q[ i ].rmin ) >
          eps * Q[ i ].rmin /
           100 ) //if it was not its optimal value, it gets reinitialized
      {
       Q[ i ].solflag = 0;
       Q[ i ].Val = Inf< double >();
       }

      Q[ i ].inter = Q[ i ].rmin;
      Q[ i ].interVal = Q[ i ].mCut.q + Q[ i ].inter * Q[ i ].mCut.m;
      }
     }
    } //if they were not all positive
   }  // if it had not already been excluded that the optimum was in this
      //  interval
   } //for all the subintervals

 // once all the LS have been carried out, we must determine which is the
 // minimum among the minima and return its position
 // it will be the most promising value, for which we will compute the DL
 // solution.
 // This could be done inside the loop; for now it is placed here to avoid
 // errors.
 // While scanning we can also do some branching: if for a given subinterval
 // the value found by the LS is >= ObjVal, i.e. the one currently considered
 // optimal,
 // we can eliminate that subinterval

 min = Inf< double >();
 // defensive default: if Q has no real subinterval (size <= 1), or every
 // Q[ i ].interVal fails the comparison below (e.g. because it is NaN),
 // the loop never assigns minpos; falling back to the last valid index
 // keeps Q[ minpos ] below in bounds instead of reading garbage
 minpos = ( Q.size() > 1 ) ? 1 : 0;

 for( i = 1 ; i < Q.size() ; i++ ) {
  if( min >= Q[ i ].interVal ) {
   min = Q[ i ].interVal;
   minpos = i;
   }

  if( Q[ i ].interVal >= BestUB ) {
   Q[ i ].branchedflag = 1;
   branchcounter++;
   }
  }

 //next line not in the original code
 //if(Q[minpos].interVal > BestLB)
 BestLB = Q[ minpos ].interVal;
 ApproxVal = Q[ minpos ].interVal;

 //std::cout << "LB=" << BestLB << std::endl;
 //std::cout << "minpos=" << minpos << std::endl;

 if( corrflag == 0 )
  return minpos;
 else
  return minconflag;
 }


/*--------------------------------------------------------------------------*/

/// updates the pair of cuts (pCut, mCut) defining the optimum of Q[ i ]
/** Replaces either Q[ i ].pCut (if the new line has nonnegative slope
 * beta) or Q[ i ].mCut (if beta is negative) with the line of slope beta
 * and intercept alpha, unless that line already coincides with one of
 * the two currently stored cuts (in which case nothing is done). This is
 * the elementary step of the exact line search performed by
 * LineSearch(): pCut and mCut are, respectively, the best "ascending"
 * and "descending" supporting lines found so far, whose intersection is
 * refined at each call until it stabilizes. */

void BenBound::UpdCut( double alpha , double beta , int i )
{
 if( Q[ i ].pCut.m != beta || Q[ i ].pCut.q != alpha ) {
  if(
   Q[ i ].mCut.m != beta ||
   Q[ i ].mCut.q !=
    alpha ) //if it already matched one of the two, updating it is pointless
  {
   if( beta >= 0 ) {
    Q[ i ].pCut.q = alpha;
    Q[ i ].pCut.m = beta;
    }
   else {
    Q[ i ].mCut.q = alpha;
    Q[ i ].mCut.m = beta;
    }
   }
  }
 }
/*<modifies one of the two optimal cuts that currently define the
  solution, depending on the slope*/

/*--------------------------------------------------------------------------*/

/// returns the position, in Q[ i ].Cuts, of the highest cut at rmin
/** Evaluates every cut currently stored in Q[ i ].Cuts at the abscissa
 * Q[ i - 1 ].rmin (the left endpoint of subinterval i) and returns the
 * position of the one attaining the largest value there, i.e., the cut
 * that is "active" (part of the upper envelope) immediately to the
 * right of the previous subinterval; used by LineSearch() in the
 * special case where all of Q[ i ]'s cuts have positive slope. */

int BenBound::whchbest( int i )
{

 int j;
 double max = -Inf< double >();
 int maxpos;
 double interc;

 for( j = 0 ; j < Q[ i ].Cuts.size() ; j++ ) {
  interc = Q[ i ].Cuts[ j ].q + Q[ i - 1 ].rmin * Q[ i ].Cuts[ j ].m;

  if( interc >= max ) {
   max = interc;
   maxpos = j;
   }
  }

 //cout<<"this is the return value of LineSearch "<<maxpos<<endl;

 return maxpos;
 }

/*--------------------------------------------------------------------------*/

/// initializes the search: sets up Q[] and the first cut(s)
/** Prepares the bounding procedure for a fresh Solve() call:
 *
 * - computes, via Limitrmin(), the interval [ limits[ 0 ] , limits[ 1 ] ]
 *   of values of r_min for which the Lagrangian subproblem is feasible;
 *   if the two bounds coincide, the whole problem collapses to a single
 *   point and is solved immediately (solvedflag = 1);
 *
 * - otherwise sets up the (initially trivial) partition Q = { Q[0], Q[1] }
 *   with Q[0].rmin/Q[1].rmin the left/right endpoints of the feasible
 *   range, and searches (by a geometric shrink controlled by myparam) for
 *   a point close to the left endpoint that is still feasible, used to
 *   initialize Q[1].inter;
 *
 * - tries to shrink the right endpoint Q[1].rmin towards the "critical"
 *   capacity crit_capc (the smallest link capacity beyond which lambda,
 *   the optimal Lagrangian multiplier, becomes 0, computed by
 *   Limitrmin()) via a binary search on the "is lambda == 0 optimal?"
 *   predicate (lagSol.is0opt()), since to its right the dual function is
 *   guaranteed nonincreasing and hence cannot improve on the value found
 *   at the shrunk endpoint;
 *
 * - (re)solves the Lagrangian subproblem at the (possibly shrunk) right
 *   endpoint Q[1].rmin and builds the corresponding left cut exactly as
 *   done in the main loop of Solve() (same arc-by-arc closed-form
 *   computation), storing it as both the first entry of Q[1].Cuts and as
 *   the initial pCut;
 *
 * - if that cut already has nonpositive slope, the right endpoint is
 *   itself the optimum and the procedure terminates immediately;
 *   otherwise Q[1].mCut is reset and the first call to LineSearch() (from
 *   Solve()) will locate the next candidate point. */

void BenBound::Inizial()
{
 int i;
 double * limits;
 double lambda;
 double sqr;
 double barl;   //\bar{l}_ij
 double SPcost; //lambda barl + f_ij rmin + lambda L /rmin
 int feas = 1;
 int in0opt;
 double zetam;
 double zetap;
 double zetaunico;
 double d;
 double tildez;
 double discrim; //if <= 0 that arc gives a zero contribution
 double ms;      //helper for the computation, partial value of the slope
 double counter = 0;
 double pointeps = 1e-9;

 DCRLagrangianSolver::LAGStatus lgstat;
 DCRLagrangianSolver::LinearCut c;
 c.m = 0;
 c.q = 0;


 limits = new double[ 2 ];

 limits =
  Limitrmin(); //limit values for which feasibility of r_min is guaranteed
 // cout<<"limit values for rmin: "<<limits[0]<<"-"<<limits[1]<<endl;

 mystep = 0.9; // before it was 0.9
 /*cout<<"I am here"<<endl;*/

 if( BenStat == 0 ) {
  if( limits[ 0 ] == limits[ 1 ] ) {
   lagSol.updrmin( limits[ 1 ] );
   //  lambda = lagSol.getLambda();
   lagSol.Solve();
   counter_ite_Lag += lagSol.getNumIte();
   lgstat = lagSol.getStatus();
   lambda = lagSol.getLambda();

   SPLabels = lagSol.getSPLabels();
   ObjVal = lambda * Flow.burst / limits[ 0 ] - lambda * Flow.deadline +
            SPLabels[ Flow.sinknode ];

   rmin = limits[ 1 ];

   solvedflag = 1;

   //cout<<"solved because limitrmin values are equal"<<endl;
   }
  else //general case where a non-trivial interval of feasibility exists
  {

   Q.resize( 2 );

   Q[ 0 ].solflag = 0;
   Q[ 1 ].solflag = 0;
   Q[ 1 ].branchedflag = 0;
   Q[ 1 ].infeasflag = 0;


   Q[ 0 ].rmin = limits[ 0 ];
   Q[ 1 ].rmin = limits[ 1 ];

   // first we look for a left endpoint close to the one found, but guaranteed
   // to be feasible,
   //we save it for now in Q[1].inter.

   double
    temppoint;  // shrink the interval, using the last infeasible value as the
                // lower bound

   Q[ 1 ].inter = Q[ 0 ].rmin;

   //look for the other feasibility point

   while( feas != 0 ) {
    temppoint = Q[ 1 ].inter;
    Q[ 1 ].inter =
     myparam * Q[ 1 ].inter +
     ( 1 - myparam ) *
      Q[ 1 ]
       .rmin;  // the left endpoint varies, the right one stays fixed at
               // Q[1].rmin

    lagSol.updrmin( Q[ 1 ].inter );
    feas = lagSol.isFeasible();
    counter++;
    //cout<<"searching multiple times?"<<endl;
    }

   Q[ 0 ].rmin = temppoint;

   //cout<<"temppoint = "<<Q[0].rmin<<endl;

   //////PRELIMINARY PHASE OF INTERVAL RESTRICTION.

   if(
    crit_capc <
    limits
     [ 1 ] )  // case where the critical point falls inside the interval, let's
              // start from there!
   {

    //cout<<"I am here in crit_capc"<<endl;
    lagSol.updrmin( crit_capc );

    in0opt = lagSol.is0opt();

    if( in0opt ==
        0 ) //in this case we can already cut away everything to its right
    {
     /*cout<<"I am inside in0opt crit_cap"<<endl;*/
     double tempd = crit_capc;
     /*cout<<"starting from"<<crit_capc<<endl;*/
     double temps = Q[ 1 ].inter;
     double cent = ( tempd + temps ) / 2;

     while( tempd - temps >
            pointeps *
             temps ) //use an LS to find the first point for which lambda = 0
     {
      lagSol.updrmin( cent );

      in0opt = lagSol.is0opt();

      if( in0opt == 0 )
       tempd = cent;
      else
       temps = cent;

      cent = ( tempd + temps ) / 2;
      }
     /*cout<<"Q[1].rmin = "<<tempd<<endl;*/
     Q[ 1 ].rmin = tempd;

     //cout<<"Q0rmin-Q1rmin = "<<Q[0].rmin<<Q[1].rmin<<endl;

     }
    else  // even if lambda != 0, we want to start the iteration just before
            // the critical point.
    {
     mystep = 0.75; //mystep = 0.75;
     Q[ 1 ].inter =
      crit_capc - 2; //stop slightly before it to avoid computation errors

     //cout<<"Q[1].inter  = "<<Q[1].inter<<endl;
     }
    } //critcap < Q[1].rmin
   else //otherwise let's still try to shrink the interval!
   {
    //cout<<"I am over here"<<endl;
    lagSol.updrmin( Q[ 1 ].rmin );

    in0opt = lagSol.is0opt();

    if( in0opt ==
        0 ) //in this case we can already cut away everything to its right
    {
     /*cout<<"I am inside in0opt"<<endl;*/

     //cout<<"entering the cut?"<<endl;
     /*cout<<"starting from"<<Q[1].rmin<<endl;*/
     double tempd = Q[ 1 ].rmin;
     double temps = Q[ 1 ].inter;
     double cent = ( tempd + temps ) / 2;

     while( tempd - temps >
            pointeps *
             temps ) //use an LS to find the first point for which lambda = 0
     {
      lagSol.updrmin( cent );

      in0opt = lagSol.is0opt();

      if( in0opt == 0 )
       tempd = cent;
      else
       temps = cent;

      cent = ( tempd + temps ) / 2;
      }
     /*cout<<"Q[1].rmin = "<<tempd<<endl;*/

     Q[ 1 ].rmin = tempd;

     } //if it fails to do so, nothing happens

    } //critcap = Q[1].rmin

   /////COMPUTATION OF THE LEFT CUT AT THE RIGHT ENDPOINT OF THE INTERVAL

   /*<in any case, after having shrunk the interval, we start populating
     it with an extremal cut*/

   lagSol.updrmin( Q[ 1 ].rmin ); //left cut for Q[1]
   lagSol.Solve();

   counter_ite_Lag = lagSol.getNumIte();
   lgstat = lagSol.getStatus();
   lambda = lagSol.getLambda();

   //cout<<"the Lagrangian solution was at lambda = "<<lambda<<endl;
   if( lgstat == 0 ) //this should actually be guaranteed by Limitrmin.
   {

    SPLabels = lagSol.getSPLabels();

    ObjVal = lambda * Flow.burst / Q[ 1 ].rmin - lambda * Flow.deadline +
             SPLabels[ Flow.sinknode ];
    HeurVal = lagSol.getHeurVal();
    HeurSOLUTION = lagSol.getHeurRSOLS();

    //cout<<ObjVal<<HeurVal<<endl;

    //before the next line was *not* commented
    //if(abs(ObjVal - HeurVal) < eps*abs(ObjVal)) {solvedflag = 1;}

    //         cout<<"objval = "<<<<" Heurval ="<<lagSol.getHeurVal()<<endl;

    //left cut
    /// computation of the value taking into account the contribution of each
    /// arc
    // the intercepts are also summed arc by arc, since some of them cause
    // discontinuities.

    for( i = 0 ; i < numLinks ; i++ ) {

     if( Links[ i ].capacity > Q[ 0 ].rmin ) {
      if( SPLabels[ Links[ i ].endnode ] > 1e200 )
       SPLabels[ Links[ i ].endnode ] = 0;
      if( SPLabels[ Links[ i ].startnode ] > 1e200 )
       SPLabels[ Links[ i ].startnode ] = 0;

      sqr = sqrt( lambda * MTU / Links[ i ].cost ); //sqrt(lambda * L / f_ij)
      d = SPLabels[ Links[ i ].startnode ] -
          SPLabels[ Links[ i ].endnode ]; //d_i-d_j
      //cout<<"d = "<<d<<endl;
      if( d < 0 ) {
       if( lambda == 0 ) //first, simplified case
       {
        zetaunico = ( -d ) / Links[ i ].cost; //(d_j - d_i)/f_ij
        // std::cout << "zetaunico - Q[0].rmin = " << zetaunico - Q[0].rmin <<
        // std::endl;
        // std::cout << "zetaunico - Links[i].capacity = " << zetaunico -
        // Links[i].capacity << std::endl;
        // std::cout << "Links[i].capacity - Q[0].rmin = " << Links[i].capacity
        // - Q[0].rmin << std::endl;
        // std::cout << "zetaunico - Links[i].capacity = " << zetaunico -
        // Links[i].capacity << std::endl;
        // std::cout << "zetaunico - Q[1].rmin = " << zetaunico - Q[1].rmin <<
        // std::endl;
        // std::cout << "Q[1].rmin - Links[i].capacity = " << Q[1].rmin -
        // Links[i].capacity << std::endl;

        if( zetaunico - Q[ 0 ].rmin > zetaunico * eps &&
            zetaunico - Links[ i ].capacity < eps * Links[ i ].capacity &&
            zetaunico - Q[ 1 ].rmin < eps * zetaunico ) {
         ms = ( Links[ i ].cost * Q[ 0 ].rmin + d ) /
              ( Q[ 0 ].rmin - Q[ 1 ].rmin );
         c.m += ms;
         c.q += -ms * Q[ 1 ].rmin;
         }

        if( Links[ i ].capacity - Q[ 0 ].rmin > Links[ i ].capacity * eps &&
            zetaunico - Links[ i ].capacity > eps * zetaunico &&
            Q[ 1 ].rmin - Links[ i ].capacity > Q[ 1 ].rmin * eps ) {
         if( zetaunico - Q[ 1 ].rmin > zetaunico * eps ) {
          ms = ( Links[ i ].cost * Links[ i ].capacity + d ) /
               ( Links[ i ].capacity - Q[ 1 ].rmin );
          c.m += ms;
          c.q += -ms * Q[ 1 ].rmin;
          }
         if( zetaunico - Q[ 1 ].rmin < Q[ 1 ].rmin * eps ) {
          ms = ( Links[ i ].cost * Q[ 0 ].rmin + d ) /
               ( Q[ 0 ].rmin - Q[ 1 ].rmin );
          c.m += ms;
          c.q += -ms * Q[ 1 ].rmin;
          }
         }
        } //lambda = 0
       else //general case
       {
        barl = MTU / Links[ i ].speed + Links[ i ].delay +
               Nodes[ Links[ i ].startnode ].delay;
        discrim =
         pow( lambda * barl + d , 2 ) - 4 * Links[ i ].cost * lambda * MTU;

        if( discrim > 0 ) //then zeta_\pm are well defined
        {
         zetap = ( -lambda * barl - d + sqrt( discrim ) ) /
                 ( 2 * Links[ i ].cost ); //larger root
         zetam = ( -lambda * barl - d - sqrt( discrim ) ) /
                 ( 2 * Links[ i ].cost ); //smaller root

         if(
          zetam - Links[ i ].capacity < eps * Links[ i ].capacity &&
          Links[ i ].capacity - sqr < eps * Links[ i ].capacity &&
          Q[ 1 ].rmin - Links[ i ].capacity >
           eps *
            Links[ i ]
             .capacity ) {  // && Q[p].inter - Links[i].capacity >
                            // eps*Links[i].capacity
          ms = ( Links[ i ].cost * Links[ i ].capacity +
                 lambda * MTU / Links[ i ].capacity + lambda * barl + d ) /
               ( Links[ i ].capacity - Q[ 1 ].rmin );
          c.m += ms;
          // c.q += (Links[i].cost * Links[i].capacity + lambda * MTU /
          // Links[i].capacity + lambda * barl + d) - ms * Links[i].capacity;
          c.q += -ms * Q[ 1 ].rmin;
          }

         if( sqr - Links[ i ].capacity < eps * sqr &&
             Links[ i ].capacity - zetap < eps * zetap &&
             Q[ 1 ].rmin >= Links[ i ].capacity ) {
          if( Q[ 1 ].rmin - zetap < Q[ 1 ].rmin * eps ) {
           if( std::abs( sqr - Q[ 1 ].rmin ) > eps * sqr ) {
            ms = ( Links[ i ].cost * sqr + lambda * MTU / sqr +
                   lambda * barl + d ) /
                 ( sqr - Q[ 1 ].rmin );
            if( Links[ i ].cost - ms > eps * Links[ i ].cost && ms >= 0.1 ) {
             tildez = sqrt( ( lambda * MTU ) / ( Links[ i ].cost - ms ) );
             c.q += Links[ i ].cost * tildez + lambda * MTU / tildez +
                    lambda * barl + d - tildez * ms;
             c.m += ms;
             }
            else {
             c.q +=
              ( Links[ i ].cost * Links[ i ].capacity +
                lambda * MTU / Links[ i ].capacity + lambda * barl + d ) -
              ms * Links[ i ].capacity;
             c.m += ms;
             }
            }
           }
          else {
           if( std::abs( sqr - Q[ 1 ].rmin ) > eps * sqr ) {
            ms = ( Links[ i ].cost * sqr + lambda * MTU / sqr +
                   lambda * barl + d ) /
                 ( sqr - Q[ 1 ].rmin );
            if( Links[ i ].cost - ms > eps * Links[ i ].cost && ms >= 0.1 ) {
             ms = ( Links[ i ].cost * sqr + lambda * MTU / sqr +
                    lambda * barl + d ) /
                  ( sqr - Q[ 1 ].rmin );
             c.m += ms;
             tildez = sqrt( ( lambda * MTU ) / ( Links[ i ].cost - ms ) );
             c.q += Links[ i ].cost * tildez + lambda * MTU / tildez +
                    lambda * barl + d - tildez * ms;
             }
            }
           }
          }

         if( Links[ i ].capacity - zetap > eps * Links[ i ].capacity ) {
          if( Links[ i ].capacity - Q[ 1 ].rmin >
              eps * Links[ i ].capacity ) {
           if( std::abs( sqr - Q[ 1 ].rmin ) > eps * sqr ) {
            ms = ( Links[ i ].cost * sqr + lambda * MTU / sqr +
                   lambda * barl + d ) /
                 ( sqr - Q[ 1 ].rmin );
            if( Links[ i ].cost - ms > eps * Links[ i ].cost && ms >= 0.1 ) {
             c.m += ms;
             tildez = sqrt( ( lambda * MTU ) / ( Links[ i ].cost - ms ) );
             c.q += Links[ i ].cost * tildez + lambda * MTU / tildez +
                    lambda * barl + d - tildez * ms;
             }
            }
           }
          else {
           if( std::abs( sqr - Q[ 1 ].rmin ) > eps * sqr ) {
            ms = ( Links[ i ].cost * sqr + lambda * MTU / sqr +
                   lambda * barl + d ) /
                 ( sqr - Q[ 1 ].rmin );
            if( Links[ i ].cost - ms > eps * Links[ i ].cost && ms >= 0.1 ) {
             c.m += ms;
             tildez = sqrt( ( lambda * MTU ) / ( Links[ i ].cost - ms ) );
             c.q += Links[ i ].cost * tildez + lambda * MTU / tildez +
                    lambda * barl + d - tildez * ms;
             }
            }
           }
          }
         }  // if the points where the arc's contribution vanishes are well
           // defined
        //otherwise the contribution will be zero.
        } //case lambda > 0
       } //if d < 0
      } //only the arcs with capacity greater than the minimum feasible one!
     } //for all the arcs.

    //finally we have to sum the contributions of the function g(z).

    ms = -lambda * Flow.burst / pow( Q[ 1 ].rmin , 2 );
    c.m += ms;
    c.q += lambda * Flow.burst / Q[ 1 ].rmin - lambda * Flow.deadline +
           SPLabels[ Flow.sinknode ] -
           Q[ 1 ].rmin * ms; //g(r_min) - r_min * ms
    // cout<<"Q[1].rmin = "<<Q[1].rmin<<" SpLabels[sinknode] =
    // "<<SPLabels[Flow.sinknode]<<endl;
    // cout<<"ms ="<<ms<<" c.m ="<<c.m<<" c.q ="<<c.q<<endl;

    //end of left cut computation, let's store it.

    Q[ 1 ].Cuts.push_back( c );

    Q[ 1 ].pCut.q =
     c.q; ///one of the cuts that defines the current point from which
    Q[ 1 ].pCut.m = c.m; // the LS starts
    }
   else
    throw( std::logic_error( "BenBound::Inizial: no cut defines the "
                             "starting point" ) );

   if(
    Q[ 1 ].pCut.m <=
    0 )  // if the slope at the endpoint is already negative we are already at
         // the optimum!
   {
    //cout<<"am I in the case pcut.m<0?"<<endl;
    Q[ 1 ].solflag = 1;

    solvedflag = 0;

    //cout<<"solved because the slope at the endpoint is negative"<<endl;

    ObjVal = lambda * Flow.burst / Q[ 1 ].rmin - lambda * Flow.deadline +
             SPLabels[ Flow.sinknode ];

    ApproxVal = ObjVal;

    rmin = Q[ 1 ].rmin;

    HeurVal = lagSol.getHeurVal();
    HeurSOLUTION = lagSol.getHeurRSOLS();
    }
   // otherwise we initialize Q[1].interval with its value and pass everything
   // to solve

   Q[ 1 ].mCut.q = 0;
   Q[ 1 ].mCut.m = 0; //give it a recognizable initialization

   } //if a non-trivial interval of feasibility existed
  } //if Benstat = 0
 //cout<<"solved flag = "<<solvedflag<<endl;
 //exit(1);
 }


/*--------------------------------------------------------------------------*/


/// computes the feasible range of r_min, and the "critical" capacity
/** Finds [ bounds[ 0 ] , bounds[ 1 ] ], a valid interval of values of
 * r_min for which the Lagrangian subproblem (solved by lagSol) is
 * feasible. Note that feasibility is *not* monotone in r_min: too small
 * a r_min makes the transmission-delay term MTU / r_min too large to
 * meet the deadline, while too large a r_min excludes links whose
 * capacity is below r_min from the (reduced) graph used by the shortest
 * path subproblem, up to possibly disconnecting it; hence feasibility
 * only holds on a "window" of r_min values, bracketed by the (sorted)
 * link capacities.
 *
 * The array of link capacities caps[] is first sorted (capSort()); a
 * linear scan from the smallest capacity upwards then finds crit_capc,
 * the smallest capacity value at which the subproblem is feasible
 * (skipping over ties): bounds[0] is set to
 * max( Flow.rate , Flow.burst / Flow.deadline ), a valid lower bound on
 * any feasible r_min. Starting from that point, a binary search over the
 * remaining (sorted) capacities then finds bounds[1], the largest
 * capacity value up to which feasibility persists (if feasibility
 * already holds at the very largest capacity, bounds[1] is simply set
 * to it without a binary search).
 *
 * Sets BenStat to OK if a feasible range was found, or to Infeasible
 * (returning nullptr) if the Lagrangian subproblem is infeasible even
 * at the largest capacity. */

double * BenBound::Limitrmin()
{
 int i , j = 0;
 int dx , sx , cent;
 int feas; //feasibility check
 //double rp, rm;

 sx = 0;
 dx = numLinks - 1;
 cent = ( dx + sx ) / 2;
 caps = new double[ numLinks ];

 for( i = 0 ; i < numLinks ; i++ )
  caps[ i ] = Links[ i ].capacity;

 capSort( sx , dx );

 ////DEBUG_DATA
 // double* red_caps = new double[numLinks];

 // red_caps[0] = caps[0];
 // j = 0;

 // for(i = 0; i < numLinks; i++){

 //   if(caps[i] != red_caps[j]){
 //     red_caps[j+1] = caps[i];
 //     j++;
 //   }
 // }
 // cout<<"types of caps:";
 // for(i = 0; i < j+1; i++) cout<<", "<<red_caps[i];

 // cout<<";"<<endl;
 // j = 0;
 //END_DEBUG_DATA

 ///It is not true that if there is no solution for the minimum capacity,
 //then the problem is necessarily empty (Waxman100_4).
 // so we take as the left endpoint the value max{Flow.rate,
 // Flow.burst/Flow.deadline}

 lagSol.LoadProblem( numNodes , numLinks , Flow , Links , Nodes , MTU ,
                     caps[ 0 ] );
 feas = lagSol.isFeasible();
 /*
    if(feas == 1){
      lagSol.updrmin(std::max(Flow.rate, Flow.burst/Flow.deadline));
      feas = lagSol.isFeasible();
    }
  */

 // caps[] holds numLinks entries (indices 0 .. numLinks - 1, i.e. 0 .. dx):
 // the loop condition must therefore guard on "j < dx", not "j < numLinks",
 // since the body increments j *before* indexing caps[ j ]. With the
 // off-by-one that used to be here, an instance whose Lagrangian
 // subproblem is infeasible at every single capacity (feas stays nonzero
 // through the last valid index, numLinks - 1) drove j one step further,
 // to numLinks, and read caps[ numLinks ] / used it via updrmin(): one
 // double past the end of the "new double[ numLinks ]" block allocated
 // above. That is a heap out-of-bounds read, so whatever happened to sit
 // in the adjacent heap memory (which differs from run to run depending
 // on allocator/ASLR state, not on the instance data) was fed to
 // isFeasible() as if it were a real capacity; on the runs where that
 // garbage looked "feasible", crit_capc below was set from it instead of
 // correctly falling through to the BenStat = Infeasible branch, silently
 // corrupting the whole Inizial()/line-search window that crit_capc seeds
 // for this flow -- while a different run, with different garbage, could
 // still stumble onto a correct answer, exactly the non-reproducibility
 // this hunts
 while( feas != 0 && j < dx ) {
  j++;

  if( caps[ j ] != caps[ j - 1 ] ) //same rmin implies same result
  {
   lagSol.updrmin( caps[ j ] );
   feas = lagSol.isFeasible();
   }
  }

 crit_capc = caps[ j ];

 if( feas == 0 ) {
  BenStat = OK;
  double * bounds = new double[ 2 ];

  if( Flow.rate < Flow.burst / Flow.deadline )
   bounds[ 0 ] =
    Flow.burst /
    Flow.deadline; //in any case this is the lower bound for r_min;
  else
   bounds[ 0 ] = Flow.rate;

  lagSol.updrmin( caps[ dx ] );
  feas = lagSol.isFeasible();

  if( feas == 0 ) //if it is fine for the last capacity we have the range;
  {
   bounds[ 1 ] = caps[ dx ];

   return bounds;
   }
  else  // otherwise binary search over the capacity vector to find the
          // minimum one
  {
   sx = j; //starting from the lower bound;

   //FIXME: this could be improved, observing that capacities often repeat
   //and by doing a check on sx+1, since often only one feasible value exists.
   while( dx - sx > 1 ) {
    cent = ( sx + dx ) / 2;
    lagSol.updrmin( caps[ cent ] );
    feas = lagSol.isFeasible();

    if( feas == 0 )
     sx = cent;
    else
     dx = cent;
    }

   lagSol.updrmin( caps[ dx ] );
   feas = lagSol.isFeasible();
   if( feas == 0 ) {
    bounds[ 1 ] = caps[ dx ];
    }
   /*lagSol.updrmin(caps[dx]+1);
        feas = lagSol.isFeasible();
        if(feas==0)
        {
          cout<<"limitrmin"<<endl;
          exit(1);
        } }*/


   else {
    bounds[ 1 ] = caps[ sx ];
    }

   /* lagSol.updrmin(caps[sx]+1);
          feas = lagSol.isFeasible();
          if(feas==0)
          {
            cout<<"sx limitrmin"<<endl;
            exit(1);
          } }*///the doubt remains that in these cases there could be
                //  some margin left

   /**FIXME: if the SPT becomes disconnected when r_min equals a certain
          capacity, it is clear that we cannot go any higher, but if the
          Lagrangian simply becomes infeasible, is it equally clear that,
          in the large gap between the last feasible capacity
          and the first infeasible one, there are no salvageable values?*/

   return bounds;
   }

  } //if the problem is feasible

 else {
  BenStat = Infeasible;
  return NULL;
  }
 }

/*--------------------------------------------------------------------------*/

/// recursively sorts caps[ sx .. dx ] in nondecreasing order
/** Randomized recursive quicksort (a random pivot in [ sx , dx ] is
 * chosen at each call to avoid worst-case behaviour on already-sorted
 * input) of the array of link capacities caps[], used by Limitrmin() to
 * enable a binary search over the feasibility of the Lagrangian
 * subproblem as a function of r_min. */

void BenBound::capSort( int sx , int dx ) //Randomized recursive quicksort.
{
 int pivot , rango;
 srand( (unsigned)time( NULL ) );

 if( sx < dx ) {
  pivot = rand() % ( dx - sx ) + sx;

  rango = Distrib( sx , pivot , dx );
  capSort( sx , rango - 1 );
  capSort( rango + 1 , dx );
  }
 }

/*--------------------------------------------------------------------------*/

/// Lomuto/Hoare-style partition of caps[ sx .. dx ] around caps[ pv ]
/** Partition step of capSort()'s quicksort: moves the pivot element
 * caps[ pv ] to the end of the range, then partitions caps[ sx .. dx ]
 * so that all elements <= the pivot come before it and all elements
 * >= the pivot come after it, finally placing the pivot in its sorted
 * position and returning that position. */

int BenBound::Distrib( int sx , int pv , int dx )  // private method for
                                                   // quicksort
{
 int i , j;

 if( pv != dx )
  Swap( pv , dx );

 i = sx;
 j = dx - 1;

 while( i <= j ) {
  while( i <= j && caps[ i ] <= caps[ dx ] )
   i++;

  while( i <= j && caps[ j ] >= caps[ dx ] )
   j--;

  if( i < j )
   Swap( i , j );
  }

 if( i != dx )
  Swap( i , dx );

 return i;
 }

/*--------------------------------------------------------------------------*/

/// swaps caps[ a ] and caps[ b ]

void BenBound::Swap( int a , int b ) //private method for quicksort
{
 double temp;

 temp = caps[ a ];
 caps[ a ] = caps[ b ];
 caps[ b ] = temp;
 }

/*--------------------------------------------------------------------------*/

/// releases all dynamically allocated memory and clears container members

void BenBound::clean_up()
{
 if( Links != 0 ) {
  delete[] Links;
  Links = 0;
  }
 if( Nodes != 0 ) {
  delete[] Nodes;
  Nodes = 0;
  }
 if( caps != 0 ) {
  delete[] caps;
  caps = 0;
  }

 for( int i = 0 ; i < Q.size() ; i++ )
  Q[ i ].Cuts.clear();
 Q.clear();

 xPlot.clear();

 SOLUTION.clear();
 HeurSOLUTION.clear();
 }


/*--------------------------------------------------------------------------*/
/*-------------------------- End File BenBound.cpp -------------------------*/
/*--------------------------------------------------------------------------*/
