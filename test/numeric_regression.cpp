// Numerical and certificate regressions; doubles in the fixture retain 17 digits.
#include "SPT.h"
#include "DCRLagrangianSolver.h"
#ifndef DCR_BASELINE
#include "DCRCertificate.h"
#include <boost/multiprecision/cpp_int.hpp>
#endif
#include <chrono>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <thread>

static int failures=0;
static void check(bool ok,const char * name) {
 std::cout << name << ',' << (ok?"PASS":"FAIL") << '\n';
 failures+=!ok;
}
using Edge=SPT::SPTLink;
static bool solve(std::vector<Edge> e,int n,int s,int t,SPT::Status expected) {
 SPT p; p.LoadProblem(n,int(e.size()),e,s,t);
 DCRDeadline deadline(.05); p.setDeadline(&deadline);
 try {p.Solve();} catch(const DCRDeadline::Expired &) {return false;}
 if(p.getStatus()!=expected) return false;
 if(expected!=SPT::OK) return p.getPath().empty() && p.getNHops()==0 && p.getLabel().empty();
 int at=s;
 for(auto a:p.getPath()) {if(a<0 || a>=int(e.size()) || e[a].startnode!=at) return false; at=e[a].endnode;}
 return at==t && int(p.getPath().size())==p.getNHops();
}
static void paths(const char * fixture) {
 std::ifstream input(fixture); std::string magic; int n,m,s,t;
 input>>magic>>n>>m>>s>>t; std::vector<Edge> e;
 for(int j=0;j<m;++j) {int index; Edge a{}; input>>index>>a.startnode>>a.endnode>>a.cost;e.push_back(a);}
 check(input && solve(e,n,s,t,SPT::Unbounded),"exact_43_27_graph");
 check(solve({{0,1,0,0},{1,2,-1,0},{2,1,1,0},{2,3,1,0}},4,0,3,SPT::OK),"exact_zero_cycle");
 const double tiny=std::numeric_limits<double>::denorm_min();
 check(solve({{0,1,0,0},{1,2,tiny,0},{2,1,0,0},{2,3,0,0}},4,0,3,SPT::OK),"tiny_positive_cycle");
 check(solve({{0,1,0,0},{1,2,-1,0},{2,1,0,0},{2,3,1,0}},4,0,3,SPT::Unbounded),"reachable_negative_cycle");
 check(solve({{0,3,1,0},{1,2,-1,0},{2,1,0,0}},4,0,3,SPT::OK),"unreachable_negative_cycle");
 check(solve({{0,1,0,0}},3,0,2,SPT::Infeasible),"unreachable_sink");
 check(solve({{0,0,-tiny,0}},1,0,0,SPT::Unbounded),"negative_self_loop");
 check(solve({{0,0,0,0}},1,0,0,SPT::OK),"zero_self_loop");
 #ifndef DCR_BASELINE
 check(solve({{0,2,0,0}},2,0,1,SPT::Error),"invalid_endpoint");
 check(solve({{0,1,NAN,0}},2,0,1,SPT::Error),"nonfinite_cost");
 #endif
 SPT reuse; reuse.LoadProblem(2,1,{{0,1,1,0}},0,1);reuse.Solve();
 DCRDeadline expired(0);reuse.setDeadline(&expired);bool interrupted=false;
 try {reuse.Solve();} catch(const DCRDeadline::Expired &) {interrupted=true;}
#ifndef DCR_BASELINE
 check(interrupted && reuse.getStatus()==SPT::Interrupted && reuse.getPath().empty() && reuse.getLabel().empty() && reuse.getNHops()==0,"deadline_clears_stale_solution");
#endif
}
static void oracle() {
 DCR::DCRNode nodes[4]{};
 DCR::DCRLink links[]={{0,1,10,10,0,0,0},{1,2,10,10,0,0,-1e-13},
                       {2,1,10,10,0,0,0},{2,3,10,10,0,0,0}};
 DCR::DCRFlow flow{0,3,0,1,10,nullptr,nullptr};
 DCRLagrangianSolver pr;pr.LoadProblem(4,4,flow,links,nodes,1,1);
 DCRDeadline deadline(.05);pr.setDeadline(&deadline);bool unknown=false, timed_out=false;
 try {pr.Solve();}
#ifndef DCR_BASELINE
 // Oracle failure is an ordinary status; it must not require an exception.
#endif
 catch(const DCRDeadline::Expired &) {timed_out=true;}
 unknown=!timed_out && pr.getStatus()==DCRLagrangianSolver::Error;
 check(unknown,"oracle_cycle_is_uncertified");
}
#ifndef DCR_BASELINE
static void rewards() {
 using boost::multiprecision::cpp_rational;
 const double a=781.5417213022222,c=1.0000000003236618,f=1000;
 const double old=f-(f+c*a);
 check(c*a+old==-1.1368683772161603e-13,"historical_cancellation_reproduced");
 const auto r=DCRCertificate::mandatory(f,c,a,40000);
 check(c*a+r.fixed>=0 && cpp_rational(r.fixed)+cpp_rational(c)*cpp_rational(a)>=0 &&
  cpp_rational(r.shift)+cpp_rational(r.fixed)<=cpp_rational(f),"reward_exact_cycle_bound");
 std::mt19937_64 gen(903);bool valid=true;int count=0;
 for(int j=0;j<10000;++j) {
  const double cost=std::ldexp(.5+double(gen()%10000)/20000,int(gen()%1801)-900)*(j%2?-1:1);
  const double lo=std::ldexp(.5+double(gen()%10000)/20000,int(gen()%1801)-900);
  const double cap=lo*2, fixed=std::ldexp(.5,int(gen()%1801)-900)*(j%3?-1:1);
  const auto r=DCRCertificate::mandatory(fixed,cost,lo,cap);
  if(!std::isfinite(r.fixed) || !std::isfinite(r.shift)) continue;
  const double endpoint=cost>=0?lo:cap;
  valid &= cpp_rational(r.fixed)+cpp_rational(cost)*cpp_rational(endpoint)>=0;
  valid &= cpp_rational(r.shift)+cpp_rational(r.fixed)<=cpp_rational(fixed);
  valid &= cost*endpoint+r.fixed>=0;
  ++count;
 }
 check(valid && count>8000,"reward_10000_exact_rational_properties");
 const double sub=std::numeric_limits<double>::denorm_min();
 const auto z=DCRCertificate::mandatory(1000,sub,sub,1);
 check(cpp_rational(z.fixed)+cpp_rational(sub)*cpp_rational(sub)>=0 &&
  cpp_rational(z.shift)+cpp_rational(z.fixed)<=1000,"reward_underflow");
}
static void partitions() {
 struct Interval {double a,b,lb;};
 auto cmp=[](const Interval & a,const Interval & b){return a.lb>b.lb;};
 std::vector<Interval> q={{1,9,4},{9,12,6}};std::make_heap(q.begin(),q.end(),cmp);
 int calls=0;bool interrupted=false;
 try {DCRCertificate::split(q,5,[&](double,double)->double {
   if(++calls==2) throw DCRDeadline::Expired{};
   return 100;
  },cmp);} catch(const DCRDeadline::Expired &) {interrupted=true;}
 check(interrupted && calls==2 && q.size()==2 && q.front().a==1 && q.front().b==9 && q.front().lb==4,"second_child_interruption_preserves_cover");
 calls=0;
 bool done=DCRCertificate::split(q,5,[&](double,double){return ++calls==1?100:-INFINITY;},cmp);
 check(!done && q.size()==2 && q.front().lb==4,"uncertified_child_preserves_cover");
 calls=0;
 done=DCRCertificate::split(q,5,[&](double,double){return ++calls==1?100:5;},cmp);
 check(done && q.size()==3 && q.front().lb==5,"both_children_commit_cover");
 done=DCRCertificate::split(q,7,[](double,double){return 0;},cmp);
 check(done && q.front().lb==5,"children_inherit_valid_parent_bound");
 DCRDeadline parent(0),child(INFINITY,&parent);bool expired=false;
 try {child.check();} catch(const DCRDeadline::Expired &) {expired=true;}
 check(expired,"nested_deadline_preserves_node_limit");
 DCRDeadline unlimited(INFINITY),shorter(0,&unlimited);expired=false;
 try {shorter.check();} catch(const DCRDeadline::Expired &) {expired=true;}
 check(expired,"nested_deadline_preserves_tranche_limit");
 // A slow nested operation must see the tranche deadline while it runs.
 DCRDeadline tranche(.001);expired=false;
 try {for(;;) tranche.check();} catch(const DCRDeadline::Expired &) {expired=true;}
 check(expired,"deadline_inside_oracle");
}
#endif
int main(int argc,char ** argv) {
 if(argc!=2) return 2;
 std::cout<<"test,result\n";
 paths(argv[1]);oracle();
#ifndef DCR_BASELINE
 rewards();partitions();
#endif
 return failures?1:0;
}
