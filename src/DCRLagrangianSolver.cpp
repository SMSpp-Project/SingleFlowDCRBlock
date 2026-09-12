/*--------------------------------------------------------------------------*/
/*-------------------------- File DCRLagrangianSolver.cpp ------------------*/
/*--------------------------------------------------------------------------*/

/** @file
 * Implementation of the DCRLagrangianSolver class [see
 * DCRLagrangianSolver.h], which computes the Lagrangian dual of the
 * Single-Flow Single-Path SRP Delay Constrained Routing problem, for a
 * given minimum per-hop rate r_min, by a two-cut line search over the
 * Lagrangian multiplier driving repeated Shortest Path Tree solves.
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

#include "DCR.h"
#include "SPT.h"
#include "DCRLagrangianSolver.h"

#include <algorithm>
#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>


/*--------------------------------------------------------------------------*/
/*--------------------- IMPLEMENTATION OF DCRLagrangianSolver---------------*/
/*--------------------------------------------------------------------------*/
// Solves, by a two-cut line search, the Lagrangian dual of DCR w.r.t. the
// delay constraint. The class keeps a copy of the network structure on
// which arcs can be closed/reopened via the ModCaps vector; out of this it
// builds a "reduced graph" that discards all arcs whose capacity is below
// r_min (given as input at load time), and out of the reduced graph it
// builds the simplified arc list that is actually fed to the SPT solver.

/*--------------------------------------------------------------------------*/
/*--------------------------- PUBLIC METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*----------------------------- CONSTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/
// constructs an "empty" solver: all fields are given their default/neutral
// value, and all the vectors used to store the reduced graph and the
// solution are given a placeholder size of 1; the object becomes usable
// only after LoadProblem() has been called.

DCRLagrangianSolver::DCRLagrangianSolver( void )
{
 eps = 1e-10; //Eps<double>();
 //cout<<"precision: "<<eps<<endl;
 ObjVal = -Inf< double >();
 InterVal = 0;
 lambda = 0;
 HeurVal = 0;

 maxCutSize = 50;

 LimitVal = 0;
 cardRedGraph = 0;

 num_ite = 0;

 lagstat = Error;
 solvedflag = 0;
 inizialflag = 0;
 reoptflag = 0;

 timer = 0;
 RedGraLinks.resize( 1 );
 RedGraPos.resize( 1 );
 Linksp.resize( 1 );
 Links = 0;
 Nodes = 0;
 ModCaps.resize( 1 );
 //XSol = 0;
 numHops = 0;
 RSol.resize( 1 );
 RSolCosts.resize( 1 );
 SPLabels.resize( 1 );
 Cuts.resize( 1 );
 solneg.resize( 1 );
 solpos.resize( 1 );
 checkSolPos.resize( 1 );
 checkSolNeg.resize( 1 );
 OptPath.resize( 1 );
 }

/*--------------------------------------------------------------------------*/
/*----------------------------------LOAD DATA ------------------------------*/
/*--------------------------------------------------------------------------*/
// loads a full instance (network + single flow + r_min) from memory: see
// LoadProblem() in DCRLagrangianSolver.h

void DCRLagrangianSolver::LoadProblem( int nnodes , int nlinks ,
                                       DCR::DCRFlow flow ,
                                       DCR::DCRLink * links ,
                                       DCR::DCRNode * nodes , double mtu ,
                                       double r_min )
{

 clean_up();

 // once any data from previous calls has been cleared, save the new
 // data and reset all the solution/status flags - - - - - - - - - - -
 numNodes = nnodes;
 numLinks = nlinks;
 MTU = mtu;
 rmin = r_min;
 ObjVal = -Inf< double >();
 lambda = 0;
 HeurVal = 0;

 num_ite = 0;

 inizialflag = 0;
 solvedflag = 0;
 reoptflag = 0;
 nonposflag = 0;
 copyDataArray( flow , links , nodes );

 RSOLS.resize( nlinks );
 HeurRSOLS.assign( nlinks , 0.0 );
 }

/*--------------------------------------------------------------------------*/

// loads a new flow while keeping the previously loaded network data: see
// LoadProblem( DCR::DCRFlow ) in DCRLagrangianSolver.h

void DCRLagrangianSolver::LoadProblem( DCR::DCRFlow flow )
{

 ////FIXME: should this reopen all the arcs that have been closed?
 //currently it does; comment out the assignment to ModCaps in the
 //second for loop below if this is not desired.
 int i;

 ///since clean_up() cannot be called here (it would destroy the
 //network structure), only what is strictly necessary is cleared:
 for( i = 0 ; i < Cuts.size() ; i++ ) {
  Cuts[ i ].RSol.clear();
  }

 Cuts.clear();


 /// install the new data and reset the solution/status flags
 ObjVal = -Inf< double >();
 lambda = 0;
 HeurVal = Inf< double >();
 HeurRSOLS.assign( numLinks , 0.0 );

 inizialflag = 0;
 solvedflag = 0;
 reoptflag = 0;

 num_ite = 0;

 ModCaps.resize( numLinks );

 Flow = flow;

 for( i = 0 ; i < numLinks ; i++ ) {
  Links[ i ].cost = Flow.costs[ i ];
  ModCaps[ i ] = Links[ i ].capacity; //reopens all the arcs
  }
 }
/*--------------------------------------------------------------------------*/
/*----------------------------------MODIFIERS-------------------------------*/
/*--------------------------------------------------------------------------*/
// changes r_min without reloading the rest of the problem data: see
// updrmin() in DCRLagrangianSolver.h

void DCRLagrangianSolver::updrmin( double r_min )

{

 // if(rmin-r_min > 0)
 //cout<<"it has decreased!"<<endl;

 if( rmin != r_min ||
     solvedflag ==
      0 ) //otherwise the problem had already been solved for this rmin
 {
  lambda = 0;
  inizialflag = 0;

  solvedflag = 0;
  reoptflag = 1; //try Reopt() in Solve()

  ObjVal = -Inf< double >();
  HeurVal = Inf< double >();
  HeurRSOLS.assign( numLinks , 0.0 );

  /*  for(int i = 0; i < Cuts.size(); i++)
     cout<<"slope of cut: "<<i<<" = "<<Cuts[i].c.m;

   cout<<endl;*/

  //discard the arcs used for SP and the reduced-graph structure
  /* ModCaps.clear();
     SPLabels.clear();
     RedGraLinks.clear();
     RedGraPos.clear();
     RSol.clear();
     RSolCosts.clear();
     Linksp.clear();*/

  rmin = r_min;
  }
 //else reoptflag = 0;
 // cout<<"lambda after updrmin"<<lambda<<endl;
 //one could require else solvedflag = 1, but it should already have
 // remained so from the previous iteration.
 }

/*--------------------------------------------------------------------------*/

void DCRLagrangianSolver::closeArcs( std::vector< int > arcs , int na )
{
 int i;

 inizialflag = 0;

 for( i = 0 ; i < na ; i++ )
  ModCaps[ arcs[ i ] ] = 0;

 //we avoid modifying the capacities in Links, because in that case
 //we would lose the information about what that capacity used to be,
 //which would make it impossible to later reopen the closed arcs.
 //The ModCaps vector keeps track of every modification to the capacities.
 }

/*--------------------------------------------------------------------------*/
// reopens a set of previously closed arcs, restoring the original capacity
// stored in Links into the corresponding entry of ModCaps

void DCRLagrangianSolver::openArcs( int * arcs , int na )
{
 int i;

 inizialflag = 0;

 for( i = 0 ; i < na ; i++ )
  ModCaps[ arcs[ i ] ] = Links[ arcs[ i ] ].capacity;
 }

/*--------------------------------------------------------------------------*/
/*--------------------------- GET METHODS-----------------------------------*/
/*--------------------------------------------------------------------------*/
// The idea behind this method is that checking feasibility of the
// Lagrangian relaxation only requires the initialization phase of the
// problem; hence, if only this is of interest and the answer turns out to
// be positive, we avoid generating all the information that solving the
// full problem would produce.
bool DCRLagrangianSolver::isFeasible( void )
{
 Inizial();

 inizialflag = 1; //signals that the problem has already been initialized,
 //so that, should it later be solved again with the same data, this
 //phase can be skipped.

 if( lagstat == 0 )
  return 0;
 else
  return 1;
 }

/*--------------------------------------------------------------------------*/
// solves SP for lambda = 0 and, if needed, for a small lambda = eps; by
// comparing the two objective values it is decided whether lambda = 0 is
// already the optimal multiplier (see is0opt() in DCRLagrangianSolver.h)

bool DCRLagrangianSolver::is0opt( void )
{
 std::vector< int > path;
 double beta_0 = Flow.burst / rmin - Flow.deadline;
 double alpha = 0;     //initialized
 double beta = beta_0; //to the constant term.
 double nulobj = 0;
 double epsobj = 0;
 int i;

 // solve the Shortest Path Tree subproblem at lambda = 0 - - - - - - - -

 lambda = 0;

 setReducedGraph();

 setSPTcosts();

 spt.updCosts( Linksp );
 spt.Solve(); //solve SP with this new lambda.

 nhops = spt.getNHops();

 path.resize( nhops );
 path = spt.getPath(); //path holds the arc indices of the optimal SP tour.

 for( i = 0 ; i < nhops ; i++ ) //compute the values of beta and alpha.
 {
  beta += MTU / Linksp[ path[ i ] ].rstar +
          MTU / RedGraLinks[ path[ i ] ].speed +
          RedGraLinks[ path[ i ] ].delay +
          Nodes[ RedGraLinks[ path[ i ] ].startnode ].delay;
  alpha += RedGraLinks[ path[ i ] ].cost * Linksp[ path[ i ] ].rstar;
  }

 if( beta <= 0 ) {
  //the path is already delay-feasible: save it as a heuristic solution
  if( HeurVal >= alpha ) {
   HeurVal = alpha;

   HeurRSOLS.assign( numLinks , 0.0 );
   for( int hi = 0 ; hi < nhops ; hi++ )
    HeurRSOLS[ RedGraPos[ path[ hi ] ] ] = Linksp[ path[ hi ] ].rstar;
   }

  return 0; //then lambda = 0 is optimal
  }
 else //otherwise try with a small lambda = eps - - - - - - - - - - - - -
 {
  nulobj = alpha;
  lambda = 1e-6;

  alpha = 0;
  beta = beta_0;

  setSPTcosts();

  spt.updCosts( Linksp );
  spt.Solve(); //solve SP with this new lambda.

  nhops = spt.getNHops();

  path.resize( nhops );
  path = spt.getPath(); //path holds the arc indices of the optimal SP tour.

  for( i = 0 ; i < nhops ; i++ ) //compute the values of beta and alpha.
  {
   beta += MTU / Linksp[ path[ i ] ].rstar +
           MTU / RedGraLinks[ path[ i ] ].speed +
           RedGraLinks[ path[ i ] ].delay +
           Nodes[ RedGraLinks[ path[ i ] ].startnode ].delay;
   alpha += RedGraLinks[ path[ i ] ].cost * Linksp[ path[ i ] ].rstar;
   }

  epsobj = alpha + lambda * beta;

  //cout<< " value at eps "<<epsobj<< " value at 0 "<<nulobj<<endl;

  // if the Lagrangian value increases with lambda, lambda = 0 is NOT
  // optimal (the dual function is still rising); otherwise it is
  if( epsobj > nulobj )
   return 1;
  else
   return 0;
  }
 }

/*--------------------------------------------------------------------------*/

double DCRLagrangianSolver::getLambda( void )
{
 //cout<<"the lambda passed by the lagrangian is"<<lambda<<endl;
 return( lambda );
 }

/*--------------------------------------------------------------------------*/


// returns the convex combination of the positive- and negative-slope
// rate solutions that (approximately) zeroes out the delay slack: the
// combination weight molt is chosen so that betaneg and betapos, the
// slopes of the two solutions, average out to (approximately) zero
std::vector< double > DCRLagrangianSolver::getCheckSol( void )
{
 double molt;
 std::vector< double > checkSol;
 int i;


 if( nonposflag == 0 )
  molt = betaneg / ( betaneg - betapos );
 else
  molt = 0;

 //molt = 0;

 checkSol.resize( solneg.size() );

 // cout<<"solution of size"<< checkSol.size();

 for( i = 0 ; i < checkSol.size() ; i++ )
  checkSol[ i ] = ( 1 - molt ) * solpos[ i ] + (molt)*solneg[ i ];

 return( checkSol );
 }


/*--------------------------------------------------------------------------*/

std::vector< int > DCRLagrangianSolver::getCheckSolPos()
{
 return( checkSolPos );
 }


/*--------------------------------------------------------------------------*/

std::vector< int > DCRLagrangianSolver::getCheckSolNeg()
{
 return( checkSolNeg );
 }


/*--------------------------------------------------------------------------*/

double
DCRLagrangianSolver::getInterVal( void ) //shouldn't really be necessary
{
 return( InterVal );
 }

/*--------------------------------------------------------------------------*/

double DCRLagrangianSolver::getOptBeta()
{
 return( optBeta );
 }

/*--------------------------------------------------------------------------*/

double DCRLagrangianSolver::getObjVal( void )
{
 return( ObjVal );
 }

/*--------------------------------------------------------------------------*/

std::vector< double > DCRLagrangianSolver::getSPLabels( void )
{
 return( SPLabels );
 }

/*--------------------------------------------------------------------------*/

int DCRLagrangianSolver::getnumHops()
{
 return( numHops );
 }

/*--------------------------------------------------------------------------*/


std::vector< int > DCRLagrangianSolver::getPath()
{
 return( OptPath );
 }


/*--------------------------------------------------------------------------*/

std::vector< double > DCRLagrangianSolver::getRSol()
{
 return( RSol );
 }

/*--------------------------------------------------------------------------*/

std::vector< double > DCRLagrangianSolver::getRSolCosts()
{
 return( RSolCosts );
 }


/*--------------------------------------------------------------------------*/


std::vector< int > DCRLagrangianSolver::getRedGraphPos()
{
 return( RedGraPos );
 }

/*--------------------------------------------------------------------------*/

int DCRLagrangianSolver::getCardRedGraph()
{
 return( cardRedGraph );
 }


/*--------------------------------------------------------------------------*/


std::vector< int > DCRLagrangianSolver::getComplGraph()
{
 return( ComplGraph );
 }


/*--------------------------------------------------------------------------*/

double DCRLagrangianSolver::getHeurVal()
{
 return( HeurVal );
 }

/*--------------------------------------------------------------------------*/
/* vector<DCRLagrangianSolver::Cut_Val> DCRLagrangianSolver::getCuts(void)
   {
        return(Cuts);
   }*/


/*--------------------------------------------------------------------------*/


double DCRLagrangianSolver::DCRgetTime( void )
{
 return( timer ? timer->Read() : 0 );
 }


/*--------------------------------------------------------------------------*/

DCRLagrangianSolver::LAGStatus DCRLagrangianSolver::getStatus( void )
{
 return( lagstat );
 }

/*--------------------------------------------------------------------------*/

int DCRLagrangianSolver::getNumIte( void )
{
 return( num_ite );
 }

/*--------------------------------------------------------------------------*/
/*------------------------- SET TIME METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

void DCRLagrangianSolver::DCRsetTime( bool timeON )
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

void DCRLagrangianSolver::DCRstartTime( void )
{
 if( ! timer )
  throw( std::logic_error( "DCRLagrangianSolver::DCRstartTime: no timer "
                           "was set" ) );

 timer->Start();
 }

/**< If timer was set or re-set, starts or re-starts timer ticking. */


/*--------------------------------------------------------------------------*/

void DCRLagrangianSolver::DCRstopTime( void )
{
 if( ! timer )
  throw( std::logic_error( "DCRLagrangianSolver::DCRstopTime: no timer "
                           "was set" ) );

 timer->Stop();
 }

/**< If timer was set, stops timer ticking. */

/*--------------------------------------------------------------------------*/

void DCRLagrangianSolver::DCRsetTimeLimit( long secs )
{
 tlimit = secs;
 }

/**< Sets a time limit for the solver.
 * \param secs long expressing timelimit in seconds */

/*--------------------------------------------------------------------------*/
/*----------------------------------SOLVE-----------------------------------*/
/*--------------------------------------------------------------------------*/

void DCRLagrangianSolver::Solve( void )
{
 int i;
 int opt = 1; //signals whether Reopt() has replaced the initialization
 std::vector< int > path;


 double releps; //max between ObjVal and 1, for the relative precision.

 double alpha;
 double beta;
 double beta_o =
  Flow.burst / rmin - Flow.deadline; //to avoid recomputing this value

 reoptflag = 0; //DEBUG: disables reoptimization
 //inizialflag == 0;

 // if a warm-start is requested, try Reopt() first; it may or may not
 // succeed in producing a valid pCut/mCut pair without a full Inizial()
 if( reoptflag == 1 ) {
  opt = Reopt();
  }

 //cout<<"opt = "<<opt<<" inizialflag = "<<inizialflag<<endl;
 if( opt == 1 ) {
  // Reopt() did not warm-start (or was not attempted): a full
  // initialization is needed, unless it was already performed (e.g.
  // by a previous call to isFeasible() on the same data)
  if( inizialflag == 0 ) {
   Inizial();
   inizialflag = 1;
   }
  }
 else {
  // Reopt() succeeded: pCut/mCut/lambda are already set up, we only
  // need the reduced graph to be current for the line search below
  setReducedGraph();
  //setSPTcosts();
  // spt.LoadProblem(numNodes, cardRedGraph, Linksp, Flow.sourcenode,
  // Flow.sinknode);
  }
 //cout<<"initialized"<<endl;
 //cout<<"lagstat ="<<lagstat<<endl;

 // two-cut line search: iterate as long as the problem is known to be
 // feasible and not already solved - - - - - - - - - - - - - - - - - -
 // use_bisection, set at the end of an iteration, forces the *next* one
 // to evaluate at the plain bisection midpoint of [ pLambda , mLambda ]
 // instead of the naive line-intersection lambda: see the note on
 // pLambda / mLambda in the header for why this fallback exists.
 // bisect_retries counts how many times this has fired so far, capped at
 // max_bisect_retries: bracket_has_room already bounds this to at most
 // ~60 halvings on its own relative-precision grounds, so the cap is a
 // pure safety margin against a case that unexpectedly fails to shrink,
 // preventing an unbounded blow-up in solve time on larger instances
 bool use_bisection = false;
 int bisect_retries = 0;
 static constexpr int max_bisect_retries = 100;
 if( lagstat == OK && solvedflag != 1 ) {
  do {
   ////variables to save the new cut: generated at every iteration
   ///because extending their lifetime would cause a deterioration of
   //the data previously stored in the Cuts vector.

   DCRLagrangianSolver::Cut_Val
    cut; //cut and solution to be inserted into the saved cuts

   // vector<double> tempRsol;  //generated solutions to be saved together with
   // the cuts
   // cout<<"beta_0 = "<<beta<<" "<<betainizial<<" "<<Flow.burst / rmin -
   // Flow.deadline<<endl;
   //intercept and slope of the new cut
   alpha = 0;     //initialized
   beta = beta_o; //to the constant term.

   if( use_bisection ) {
    // the previous iteration's naive intersection converged to a lambda
    // where beta is not actually close to 0: fall back to plain
    // bisection of the true [ pLambda , mLambda ] bracket for this one
    // trial, in case the spurious convergence was only an artifact of a
    // loose bracket rather than a genuine gap (see the note on
    // pLambda / mLambda)
    lambda = ( pLambda + mLambda ) / 2;
    use_bisection = false;
    }
   else
    lambda =
     ( mCut.q - pCut.q ) /
     ( pCut.m -
       mCut.m ); //lambda is set to the intersection point of the two cuts.

   if( lambda <
       0 ) //also add to the line search the linear constraint lambda = 0;
   {
    lambda = 0; /*cout<<"I am forcing lambda = 0";*/
    }

   //cout<<"inside lagrange, how does lambda evolve?"<<lambda<<endl;
   InterVal =
    mCut.q +
    lambda *
     mCut
      .m;  // old line (either one, it makes no difference) evaluated at the
           // new point.
          //cout<<"InterVal = "<<InterVal<<endl;
          //vector<double> tempRsol;

   //cout<<"1beta_0 ="<<beta_o<<endl;
   setSPTcosts(); //new lambda, new costs.
   spt.updCosts( Linksp );
   spt.Solve(); //solve SP with this new lambda.

   nhops = spt.getNHops();

   path.resize( nhops );
   cut.RSol.resize( nhops );
   //tempRsol.resize(nhops);


   // cout<<"2beta_0 ="<<beta_o<<endl;
   path = spt.getPath(); //path holds the arc indices of the optimal SP tour.

   for( i = 0 ; i < nhops ; i++ ) //compute the new values of beta and alpha.
   {
    beta += MTU / Linksp[ path[ i ] ].rstar +
            MTU / RedGraLinks[ path[ i ] ].speed +
            RedGraLinks[ path[ i ] ].delay +
            Nodes[ RedGraLinks[ path[ i ] ].startnode ].delay;
    alpha += RedGraLinks[ path[ i ] ].cost * Linksp[ path[ i ] ].rstar;
    //tempRsol[i] = Linksp[path[i]].rstar;
    cut.RSol[ i ] = Linksp[ path[ i ] ].rstar;
    }

   optBeta = beta;

   if( beta < -1e-20 && alpha <= HeurVal ) {
    HeurVal = alpha;

    // save the routing that attains this HeurVal *now*, from the path/
    // nhops of *this* iteration: HeurRSOLS must not be filled in after
    // the do-while loop below exits, the way RSOLS is, because a later
    // iteration can overwrite path/nhops with a different (possibly
    // delay-infeasible) routing before the loop ends, leaving HeurVal
    // and the routing that justifies it permanently out of sync
    HeurRSOLS.assign( numLinks , 0.0 );
    for( int hi = 0 ; hi < nhops ; hi++ )
     HeurRSOLS[ RedGraPos[ path[ hi ] ] ] = Linksp[ path[ hi ] ].rstar;
    } //save the best delay-feasible primal solution found so far

   /*
         /////save the data for the CHECK against CPLEX!!!!
         if(beta < 0)
         {
           solneg.resize(nhops);

           checkSolNeg.resize(nhops);

           for(i = 0; i < nhops; i++)
           {
            solneg[i] = Linksp[path[i]].rstar;
            checkSolNeg[i] = RedGraPos[path[i]];
           }
           betaneg = beta;
         }
         else
         {
           solpos.resize(nhops);
           checkSolPos.resize(nhops);

           for(i = 0; i < nhops; i++)
            {
              solpos[i] = Linksp[path[i]].rstar;
             checkSolPos[i] = RedGraPos[path[i]];
           }
           betapos = beta;

         }
         ///end of the saving, comment out in the future to save memory
    */
   cut.c.m = beta;
   cut.c.q = alpha;

   cut.RSolsize = nhops;
   cut.rmin = rmin;

   if( Cuts.size() < maxCutSize )
    Cuts.push_back( cut ); //save the new cut

   UpdCut( alpha , beta ,
          lambda ); //update the value of one of the two optimal cuts.

   ObjVal =
    alpha +
    lambda *
     beta;  // new line evaluated at the point, i.e. the objective function
            // value.

   if( InterVal > 1 )
    releps = InterVal;
   else
    releps = 1;

   //cout<<"in solve absolute difference"<<(InterVal - ObjVal)<<endl;
   //cout<<"relative threshold"<<eps*releps<<endl<<endl;
   num_ite++;

   // the two-cut model (InterVal) and the fresh evaluation (ObjVal) can
   // agree -- satisfying the convergence test below -- even though beta
   // is nowhere near 0: e.g. once every arc on the shortest path has its
   // rate clamped at its own capacity, lambda no longer moves (alpha,
   // beta) at all, so re-evaluating there trivially "confirms" a bracket
   // that was never actually narrowed onto the true beta == 0 crossing
   // (see the note on pLambda / mLambda). Detect this and, as long as the
   // true bracket still has room, retry once at its plain bisection
   // midpoint instead of accepting the spurious convergence -- this
   // recovers cases where the bracket was simply too loose; a further
   // "converged with room left" event on the very next iteration (the
   // bracket having genuinely halved) instead signals that the gap is
   // real (typically a discontinuous switch to a cheaper-but-infeasible
   // shortest path) and there is nothing more bisection can do about it
   if( ( InterVal - ObjVal ) <= eps * releps ) {
    const bool beta_not_tight =
     std::abs( beta ) > 1e-6 * std::max( 1.0 , std::abs( Flow.deadline ) );
    const bool bracket_has_room =
     std::abs( pLambda - mLambda ) >
     1e-9 * std::max( 1.0 , std::abs( mLambda ) );
    use_bisection = beta_not_tight && bracket_has_room &&
                    ( bisect_retries < max_bisect_retries );
    if( use_bisection )
     bisect_retries++;
    }

   } while(
   ( ( InterVal - ObjVal ) >
     ( eps *
       releps ) ) ||  /// gap > epsilon * max{1,InterVal}: convergence test
                      /// of the line search
   use_bisection );

  //cout<<" eps = "<<eps<<endl;
  //cout<<"in solve relative difference"<<(InterVal - ObjVal)/InterVal<<endl;
  //cout<<endl<<endl<<"beta = "<<beta<<endl<<endl;
  //cout<<endl<<endl<<"alpha = "<<alpha<<endl<<endl;

  //cout<<"number of line search iterations"<<counter<<endl;
  //cout<<"past the iteration?"<<endl;

  // save the node potentials (dual solution) of the last SP solve - - - -
  SPLabels.resize( numNodes );

  std::vector< double > labels;

  labels = spt.getLabel(); //save the potentials.

  for( i = 0 ; i < numNodes ; i++ )
   SPLabels[ i ] = labels[ i ];

  //cout<<"print?"<<endl;

  RSol.resize( nhops );
  RSolCosts.resize( nhops );

  for( i = 0 ; i < numLinks ; i++ )
   RSOLS[ i ] = 0.0;

  for( i = 0 ; i < nhops ;
       i++ ) //save the (r_ij, f_ij) pairs of the optimal solution
  {
   RSol[ i ] = Linksp[ path[ i ] ].rstar;
   //std::cout << i << "," << RSol[i] << std::endl;
   RSolCosts[ i ] = Linksp[ path[ i ] ].cost;
        }

  // also fill in RSOLS, the rate solution indexed over *all* the arcs
  // of the (unreduced) network, leaving 0 wherever an arc is unused.
  // path[ j ] is an index into the REDUCED graph (Linksp / RedGraLinks,
  // sized cardRedGraph: only the arcs whose capacity supports the current
  // rate); it must be translated back to the original arc index via
  // RedGraPos[] before comparing against i, which ranges over the full,
  // unreduced graph -- comparing i directly to path[ j ] matches an
  // original arc index against a reduced-graph one, which are unrelated
  // numbering spaces, and reports essentially arbitrary arcs as "used"
  for( int i = 0 ; i < numLinks ; i++ )
   for( int j = 0 ; j < nhops ; j++ )
    if( i == RedGraPos[ path[ j ] ] )
     RSOLS[ i ] = RSol[ j ];


  numHops = nhops;
  solvedflag = 1;
  } //if it had not already been solved and is not already deemed infeasible
 //if(solvedflag == 1)
 //cout<<"number of saved cuts"<<Cuts.size()<<endl;
 //cout<<"optimal lambda: "<<lambda<<" solved with status: "<<lagstat<<endl;

 path.clear();
 //cout<<"I did not perform the iteration"<<endl;
 }

/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/
// destructor of DCRLagrangianSolver: just calls clean_up()

DCRLagrangianSolver::~DCRLagrangianSolver()
{
 clean_up();
 }

/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/
// deep-copies the network and flow data into the internal Links/Nodes
// arrays, and initializes ModCaps to the original arc capacities: see
// copyDataArray() in DCRLagrangianSolver.h

void DCRLagrangianSolver::copyDataArray( DCR::DCRFlow flow ,
                                         DCR::DCRLink * links ,
                                         DCR::DCRNode * nodes )
{
 int i;
 Links = new DCR::DCRLink[ numLinks ];
 ModCaps.resize( numLinks );
 Nodes = new DCR::DCRNode[ numNodes ];

 Flow = flow;

 for( i = 0 ; i < numLinks ; i++ ) {
  Links[ i ].startnode = links[ i ].startnode;
  Links[ i ].endnode = links[ i ].endnode;
  Links[ i ].speed = links[ i ].speed;
  Links[ i ].capacity = links[ i ].capacity;
  Links[ i ].delay = links[ i ].delay;
  Links[ i ].cost = links[ i ].cost; //Flow.costs[i];
  ModCaps[ i ] = Links[ i ].capacity;
  }

 //cout<<"a capacity="<<Links[0].capacity<<endl;
 //cout<<"a cost="<<Links[0].cost<<endl;
 for( i = 0 ; i < numNodes ; i++ ) {
  Nodes[ i ].delay = nodes[ i ].delay;
  }
 }


/*--------------------------------------------------------------------------*/

// computes the trivial worst-case upper bound LimitVal: see getLimitVal()
// in DCRLagrangianSolver.h
void DCRLagrangianSolver::getLimitVal( void )
{
 int i;

 LimitVal = 0;

 for( i = 0 ; i < cardRedGraph ; i++ ) {
  if( RedGraLinks[ i ].cost > 0 )
   LimitVal +=
    RedGraLinks[ i ].cost *
    RedGraLinks[ i ].capacity; //saturate all the positive-cost arcs
  } //and sum up their costs.
 }
/*upper bound on the value of any feasible solution*/

/*--------------------------------------------------------------------------*/
// returns the Lagrangian cost of reserving rate r on arc lindex of the
// reduced graph, at the current lambda: see getCost() in
// DCRLagrangianSolver.h

double DCRLagrangianSolver::getCost( double r , int lindex )
{
 double c;
 double barl = MTU / RedGraLinks[ lindex ].speed +
               RedGraLinks[ lindex ].delay +
               Nodes[ RedGraLinks[ lindex ].startnode ].delay;
 c = lambda * barl + RedGraLinks[ lindex ].cost * r + lambda * MTU / r;

 return c;
 }
/*<returns the costs for the arc lindex of SP as a function of r^*_ij,
  passed as input as r*/

/*--------------------------------------------------------------------------*/

// tries to warm-start the line search from the cached cuts after r_min has
// changed: see Reopt() in DCRLagrangianSolver.h
int DCRLagrangianSolver::Reopt( void )
{

 /*FIXME: find values of lambda for which the intersection between the
     cuts lies below the objective function value.
     Are we sure the cuts are saved in a safe way? (it would seem so)*/

 int isnegcut = 1; //flag telling whether a negative-slope cut is available
 int negpos;       //position of the negative-slope cut in Cuts
 int isposcut = 1; //same for positive-slope cuts
 int pospos;
 double epsrel; //factor for the relative precision
 int i , k;
 int delflag = 0; //flag for discarding solutions that are no longer feasible
 double someinter;
 double value;
 double approx;
 double min = Inf< double >();
 int minpos;
 int counterdeleted = 0;
 int maxiter = 0;

 //cout<<"cuts present"<<Cuts.size()<<endl;
 if( Cuts.size() <= 1 )
  return 1; //a single cut, or none, is too little to reoptimize on

 for( i = 0 ; i < Cuts.size() ; i++ )  // check which solutions are still
                                       // feasible
 {
  //cout<<"nhops of the solution of cut "<<i<<" = "<<Cuts[i].RSolsize<<endl;
  for( k = 0 ; k < Cuts[ i ].RSolsize ; k++ ) {
   if( Cuts[ i ].RSol[ k ] <
       rmin ) //it is enough that this happens for a single r_ij
   {
    //cout<<Cuts[i].RSol[k]<<"-"<<rmin<<endl;
    delflag = 1; //in which case the solution is no longer feasible
    }
   }
  if( delflag == 1 ) {
   counterdeleted++;
   Cuts.erase( Cuts.begin() + i ); //so the cut is discarded
   i--;
   delflag = 0;
   }
  else //if the cut is still valid, check its slope
  //FIXME: a sensitivity margin could be added here
  {
   //correct the slope of the cut for the new value of rmin
   Cuts[ i ].c.m =
    Cuts[ i ].c.m - Flow.burst / Cuts[ i ].rmin + Flow.burst / rmin;
   //and record that the correction for this cut has been applied
   Cuts[ i ].rmin = rmin;

   if( Cuts[ i ].c.m < 0 ) {
    isnegcut = 0;
    negpos = i;
    }
   else {
    if( Cuts[ i ].c.m >= 0 ) //should be redundant
    {
     isposcut = 0;
     pospos = i;
     }
    }
   }
  }
 //cout<<"cuts deleted"<<counterdeleted<<endl;
 if( Cuts.size() <= 1 )
  return 1;
 //cout<<"there remain "<<Cuts.size()<<endl;

 if( isnegcut == 1 || isposcut == 1 )
  return 1; //if no positive- or no negative-slope cut survived, give up

 ////line search for the optimal intersection among the still-feasible cuts

 pCut.m = Cuts[ pospos ].c.m;
 pCut.q = Cuts[ pospos ].c.q;
 mCut.m = Cuts[ negpos ].c.m;
 mCut.q = Cuts[ negpos ].c.q;

 do {
  lambda = ( mCut.q - pCut.q ) / ( pCut.m - mCut.m );

  ////FIXME: understand why this happens!!!!!!!
  if( lambda < 0 ||
      std::abs( pCut.m - mCut.m ) < 1e-20 ) { //cout<<"Bug in reopt"<<endl;
   lambda = 0;
   return 1;
   }
  approx = mCut.q + lambda * mCut.m;

  for( i = 0 ; i < Cuts.size() ;
       i++ ) //to compute the function value at that lambda,
  {          //find which cut gives the lowest value.
   someinter = Cuts[ i ].c.q + lambda * Cuts[ i ].c.m;

   if( min >= someinter ) {
    min = someinter;
    minpos = i;
    }
        }

  value = Cuts[ minpos ].c.q + lambda * Cuts[ minpos ].c.m;

  UpdCut( Cuts[ minpos ].c.q , Cuts[ minpos ].c.m , lambda );

  if( approx > 1 )
   epsrel = approx;
  else
   epsrel = 1;

  //cout<<"here is the loop"<<endl;
  maxiter++;
  //cout<<"approx - value"<<approx - value<<endl;
  } while( approx - value > eps * epsrel &&
          maxiter <= 30 ); // convergence test, capped at 30 iterations
 //FIXME: loops with garr199904, launched by BenBound, at flow 78
 // cout<<"approx || value : "<<approx<<"||"<<value<<endl;
 //cout<<"niter = "<<maxiter<<endl;
 if( maxiter >= 30 ) {
  std::cout << "watch out, reopt has gone into a loop" << std::endl;
  return 1; // did not converge: fall back to a full Inizial()
  }
 // cout<<"reopt ends with value = "<<value<<endl;
 //cout<<"public cuts:"<<pCut.m<<"||"<<mCut.m<<endl;
 return 0;
 }

/*--------------------------------------------------------------------------*/

// initializes the line search: builds the reduced graph, solves the
// Shortest Path Tree subproblem at lambda = 0 and, depending on the sign
// of the resulting delay slack, either declares that lambda = 0 is
// already optimal or searches (by doubling lambda) for a first
// negative-slope cut; sets lagstat accordingly. See Inizial() in
// DCRLagrangianSolver.h
void DCRLagrangianSolver::Inizial( void )
{
 int i , j;
 std::vector< int > path;
 // vector<double> tempRsol;
 double beta_0 = Flow.burst / rmin - Flow.deadline;
 double alpha = 0;     //initialized
 double beta = beta_0; //to the constant term.

 lambda = 0; // the first trial of the two-cut search is always at
             // lambda = 0, made explicit here (rather than relying on the
             // caller to have left it there) since it is recorded as
             // pLambda below

 //cout<<"beta = "<<beta<<endl;
 setReducedGraph(); //work on the reduced graph.

 setSPTcosts();

 //cout<<"cardRedGraph="<<cardRedGraph<<endl;
 //cout<<Flow.sourcenode<<Flow.sinknode<<endl;

 spt.LoadProblem( numNodes , cardRedGraph , Linksp , Flow.sourcenode ,
                  Flow.sinknode );
 spt.Solve();

 spstat = spt.getStatus();
 //cout<<"spstat is "<<spstat<<endl;
 if(
  spstat ==
  0 ) //the reduced graph has no negative-cost cycle and is not disconnected.
 {
  getLimitVal();

  nhops = spt.getNHops();
  path.resize( nhops );
  path = spt.getPath(); //path holds the arc indices of the optimal SP tour

  // cout<<"path in the lagrangian and respective rstar: ";

  for( i = 0 ; i < nhops ; i++ ) //compute the values of beta and alpha.
  {
   beta += MTU / Linksp[ path[ i ] ].rstar +
           MTU / RedGraLinks[ path[ i ] ].speed +
           RedGraLinks[ path[ i ] ].delay +
           Nodes[ RedGraLinks[ path[ i ] ].startnode ].delay;
   alpha += RedGraLinks[ path[ i ] ].cost * Linksp[ path[ i ] ].rstar;
   // cout<<"barl in lagr for arc "<<i<<" is "<<MTU /
   // RedGraLinks[path[i]].speed + RedGraLinks[path[i]].delay +
   // Nodes[RedGraLinks[path[i]].startnode].delay<<endl;
   // cout<<"and its rstar is "<<Linksp[path[i]].rstar<<endl;
   }
  // cout<<endl;
  //cout<<"alpha= "<<alpha<<endl;
  //cout<<"beta = "<<beta<<endl;

  UpdCut( alpha , beta , lambda ); //update the first optimal cut.

  /////check whether 0 is already the optimal lambda
  if( beta <= 0 ) {
   DCRLagrangianSolver::Cut_Val cut;
   //cout<<"negative slope"<<endl;
   lagstat = OK;
   InterVal = alpha;
   ObjVal = alpha;
   lambda = 0;
   optBeta = beta;

   HeurVal =
    alpha; // save the objective function value at this beta<=0 solution

   cut.RSol.resize( nhops );
   path.resize( nhops );
   OptPath.resize( nhops );

   path = spt.getPath();


   RSol.resize( nhops );
   RSolCosts.resize( nhops );

   /*
          ////SAVE THE DATA FOR THE CHECK AGAINST CPLEX
          solneg.resize(nhops);
          checkSolNeg.resize(nhops);

          for(i = 0; i < nhops; i++)
           {
            solneg[i] = Linksp[path[i]].rstar;
            checkSolNeg[i] = RedGraPos[path[i]];
           }

          betaneg = beta;
          nonposflag = 1;
          //end of saving
    */

   for( i = 0 ; i < nhops ;
        i++ ) //save the (r_ij, f_ij) pairs of the optimal solution
   {
    RSol[ i ] = Linksp[ path[ i ] ].rstar;
    cut.RSol[ i ] = Linksp[ path[ i ] ].rstar;
    RSolCosts[ i ] = Linksp[ path[ i ] ].cost;
    OptPath[ i ] = path[ i ];
         }

   for( i = 0 ; i < numLinks ; i++ )
    RSOLS[ i ] = 0.0;

   // path[ j ] is a reduced-graph index: translate it back to the
   // original arc index via RedGraPos[] before comparing against i,
   // exactly as in the analogous loop at the end of Solve()
   for( int i = 0 ; i < numLinks ; i++ )
    for( int j = 0 ; j < nhops ; j++ )
     if( i == RedGraPos[ path[ j ] ] )
      RSOLS[ i ] = RSol[ j ];

   // this beta <= 0 solution is delay-feasible by construction, and it is
   // exactly the one HeurVal above was just set from: HeurRSOLS must be
   // kept in lockstep with it, the same as RSOLS is with the (possibly
   // different, not necessarily feasible) value returned by getRSOLS()
   HeurRSOLS = RSOLS;

   cut.RSolsize = nhops;
   cut.c.m = beta;
   cut.c.q = alpha;

   cut.rmin = rmin;

   if( Cuts.size() < maxCutSize )
    Cuts.push_back( cut ); //save the negative-slope cut.

   SPLabels.resize( numNodes );
   std::vector< double > labels;

   labels = spt.getLabel(); //save the potentials.

   for( i = 0 ; i < numNodes ; i++ )
    SPLabels[ i ] = labels[ i ];

   solvedflag = 1;
   }
  else //////search for a negative-slope cut, since beta > 0 at lambda = 0
  {
   //cout<<"searching for the cut with negative slope"<<endl;
   lambda = 1;

   // repeatedly double lambda and re-solve SP, until either a
   // delay-feasible (beta <= 0) path is found, the running
   // Lagrangian value exceeds the trivial upper bound LimitVal
   // (declaring infeasibility), or SP itself fails
   while( beta > 0 && ObjVal <= LimitVal && spstat == 0 ) {
    alpha = 0; //initialized at every iteration
    beta = beta_0;

    setSPTcosts();
    spt.updCosts( Linksp );
    spt.Solve();
    spstat = spt.getStatus();
    nhops = spt.getNHops();

    path.resize( nhops );
    path = spt.getPath(); //path holds the arc indices of the optimal SP tour
    //cout<<path[0]<<path[1]<<endl;
    for( i = 0 ; i < nhops ; i++ ) //compute the values of beta and alpha.
    {
     beta += MTU / Linksp[ path[ i ] ].rstar +
             MTU / RedGraLinks[ path[ i ] ].speed +
             RedGraLinks[ path[ i ] ].delay +
             Nodes[ RedGraLinks[ path[ i ] ].startnode ].delay;
     alpha += RedGraLinks[ path[ i ] ].cost * Linksp[ path[ i ] ].rstar;
     }

    //cout<<"beta = "<<beta<<endl;
    //cout<<"ObjVal ="<<ObjVal<<endl;

    ObjVal = alpha + lambda * beta;
    lambda = 2 * lambda; //other growth sequences could be chosen as well
    } // while searching for a negative beta

   DCRLagrangianSolver::Cut_Val cut;
   //tempRsol.resize(nhops);
   cut.RSol.resize( nhops );

   for( i = 0 ; i < nhops ; i++ ) {
    //tempRsol[i] = Linksp[path[i]].rstar;
    cut.RSol[ i ] = Linksp[ path[ i ] ].rstar;
    }

   /*/
          ////SAVE THE DATA FOR THE CHECK AGAINST CPLEX
          solneg.resize(nhops);
          checkSolNeg.resize(nhops);

          for(i = 0; i < nhops; i++)
           {
            solneg[i] = Linksp[path[i]].rstar;
            checkSolNeg[i] = RedGraPos[path[i]];
           }

          betaneg = beta;
          nonposflag = 1;
          //end of the saving
    */

   cut.c.m = beta;
   cut.c.q = alpha;

   /*for(i = 0; i < nhops; i++)
           cut.RSol[i] = tempRsol[i];*/

   cut.RSolsize = nhops;
   cut.rmin = rmin;

   if( Cuts.size() < maxCutSize )
    Cuts.push_back( cut ); //save the negative-slope cut.

   // lambda has already been doubled (line above) past the value that
   // actually produced this (alpha, beta): halve it back to record the
   // multiplier that was really used
   UpdCut( alpha , beta , lambda / 2 ); //update the second optimal cut.
   //cout<<"LimitVal="<<LimitVal<<endl;
   //cout<<"ObjVal="<<ObjVal<<endl;

   if( ObjVal > LimitVal ) {
    lagstat = Infeasible;
    } // then the Lagrangian relaxation is empty. ?? but if we
     //started with a cut, there should be at least one point;
     //this may signal an increasing function
   else
    lagstat = OK;
   }
  //cout<<"beta = "<<beta<<endl;
  } // if SP succeeded
 else
  lagstat = Infeasible;
 path.clear();
 }


/*--------------------------------------------------------------------------*/

/// /FIXME: it would be better to do this with a single scan, but it is only
/// done once
// (re)builds the reduced graph, keeping only the arcs whose (working)
// capacity is >= rmin: see setReducedGraph() in DCRLagrangianSolver.h
void DCRLagrangianSolver::setReducedGraph()
{
 int i , j = 0 , h = 0;
 cardRedGraph = 0; //for the subsequent calls

 for( i = 0 ; i < numLinks ; i++ ) {
  if( rmin <=
      ModCaps[ i ] ) //these are the ones we will work on (r_min <= c_ij)
   cardRedGraph++;
  }

 //FIXME: do I really need RedGraLinks, now that I save its positions?
 Linksp.resize( cardRedGraph );
 RedGraLinks.resize( cardRedGraph );
 RedGraPos.resize( cardRedGraph );
 //ComplGraph.resize(numLinks - cardRedGraph);

 for( i = 0 ; i < numLinks ; i++ ) {
  if( rmin <=
      ModCaps[ i ] ) //otherwise SP must not consider them (i.e. x_ij = 0).
  {
   RedGraLinks[ j ].startnode = Links[ i ].startnode;
   RedGraLinks[ j ].endnode = Links[ i ].endnode;
   RedGraLinks[ j ].capacity = Links[ i ].capacity;
   RedGraLinks[ j ].speed = Links[ i ].speed;
   RedGraLinks[ j ].delay = Links[ i ].delay;
   RedGraLinks[ j ].cost = Links[ i ].cost;
   Linksp[ j ].startnode = Links[ i ].startnode;
   Linksp[ j ].endnode = Links[ i ].endnode;
   RedGraPos[ j ] = i;
   j++;
   }
  /* else
       {
         ComplGraph[h] = i;
         h++;
       }*/
  }
 }

/*--------------------------------------------------------------------------*/
// for the current lambda, computes the locally optimal rate r* and the
// corresponding Lagrangian cost of every arc of the reduced graph: see
// setSPTcosts() in DCRLagrangianSolver.h

void DCRLagrangianSolver::setSPTcosts()
{
 int i;
 double sqr;

 for( i = 0 ; i < cardRedGraph ; i++ ) {
  ///choice of the costs according to the various cases, see paragraph 2.
  //cout<<"set STPcost"<<endl;
  if( RedGraLinks[ i ].cost < 0 ) {
   // negative-cost arc: the Lagrangian cost is decreasing in r,
   // so the minimum over [rmin, capacity] is at r* = capacity
   Linksp[ i ].rstar = RedGraLinks[ i ].capacity;
   Linksp[ i ].cost = getCost( RedGraLinks[ i ].capacity , i );
   }
  else {
   // non-negative-cost arc: the unconstrained minimizer of
   // cost*r + lambda*MTU/r is sqrt(lambda*MTU/cost); clip it to
   // the feasible range [rmin, capacity]
   sqr = sqrt( lambda * MTU / RedGraLinks[ i ].cost );

   if( sqr < rmin ) {
    Linksp[ i ].rstar = rmin;
    Linksp[ i ].cost = getCost( rmin , i );
    }
   else {
    if( RedGraLinks[ i ].capacity < sqr ) {
     Linksp[ i ].rstar = RedGraLinks[ i ].capacity;
     Linksp[ i ].cost = getCost( RedGraLinks[ i ].capacity , i );
     }
    else { //cout<<"I choose"<<sqr<<endl;
     Linksp[ i ].rstar = sqr;
     Linksp[ i ].cost = getCost( sqr , i );
     } //case rmin < sqr < c_ij
    }
   } //end of the case split on the sign of the cost.
  } //for all reduced graph arcs
 }

/*--------------------------------------------------------------------------*/
// releases all dynamically allocated memory and resets the solution data:
// see clean_up() in DCRLagrangianSolver.h

void DCRLagrangianSolver::clean_up()
{
 int i;

 if( Links != 0 ) {
  delete[] Links;
  Links = 0;
  }
 if( Nodes != 0 ) {
  delete[] Nodes;
  Nodes = 0;
  }
 //if(timer != 0) {delete [] timer; timer = 0;}
 ModCaps.clear();
 SPLabels.clear();
 RedGraLinks.clear();
 RedGraPos.clear();
 RSol.clear();
 RSolCosts.clear();
 Linksp.clear();
 solpos.clear();
 solneg.clear();
 checkSolPos.clear();
 checkSolNeg.clear();
 OptPath.clear();
 RSOLS.clear();
 HeurRSOLS.clear();

 for( i = 0 ; i < Cuts.size() ; i++ ) {
  Cuts[ i ].RSol.clear();
  }

 Cuts.clear();
 }

/*--------------------------------------------------------------------------*/
/*----------------------------------MODIFIERS-------------------------------*/
/*--------------------------------------------------------------------------*/


// updates pCut or mCut, depending on the sign of the slope, with the new
// (alpha, beta) cut: see UpdCut() in DCRLagrangianSolver.h

void DCRLagrangianSolver::UpdCut( double alpha , double beta ,
                                  double at_lambda )
{
 if( beta >= 0 ) {
  pCut.m = beta;
  pCut.q = alpha;
  pLambda = at_lambda;
  }

 else {
  mCut.m = beta;
  mCut.q = alpha;
  mLambda = at_lambda;
  }
 }
/*<updates one of the two optimal cuts defining the current solution,
 * depending on its slope*/


/*--------------------------------------------------------------------------*/
/*-------------------- End File DCRLagrangianSolver.cpp --------------------*/
/*--------------------------------------------------------------------------*/
