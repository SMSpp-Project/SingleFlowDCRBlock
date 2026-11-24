#ifndef DCR_CERTIFICATE_H
#define DCR_CERTIFICATE_H

#include <algorithm>
#include <cmath>
#include <limits>

namespace DCRCertificate {
// Outward rounding by one adjacent value encloses the exact result of each
// individual IEEE-754 operation under round-to-nearest, including subnormals.
// This is an arithmetic error bound, not a tolerance or a fixed reward nudge.
inline double down(double x) {
 return std::nextafter(x,-std::numeric_limits<double>::infinity());
}
inline double add(double a,double b) { return down(a+b); }

struct Reward { double fixed, shift; };
inline Reward mandatory(double fixed,double cost,double a,double capacity) {
 const double endpoint = cost >= 0 ? a : capacity;
 const double product = cost * endpoint;
 // Use a lower enclosure of c*r_min directly, avoiding the cancellation
 // fixed - (fixed + c*r_min). Then fixed' + c*r >= 0 on the entire domain.
 const double lower = (cost == 0 || endpoint == 0) ? 0 : down(product);
 const double transformed = -lower;
 // Account for the ACTUAL transformed coefficient, not the nominal reward:
 // shift <= fixed - fixed'. For mandatory x=1 this preserves the bound.
 return {transformed,down(fixed-transformed)};
}

// Strong certificate guarantee: do not remove a covering parent until BOTH
// children are certified. Exceptions, NaNs and -infinity leave it untouched.
// Children inherit the parent's valid bound as well as their new certificate.
template<class Heap,class Oracle,class Compare>
bool split(Heap & heap,double middle,Oracle oracle,Compare compare) {
 const auto parent = heap.front();
 const double left = oracle(parent.a,middle);
 if(std::isnan(left) || left == -std::numeric_limits<double>::infinity()) return false;
 const double right = oracle(middle,parent.b);
 if(std::isnan(right) || right == -std::numeric_limits<double>::infinity()) return false;
 heap.reserve(heap.size()+1);
 std::pop_heap(heap.begin(),heap.end(),compare);
 heap.pop_back();
 heap.push_back({parent.a,middle,std::max(parent.lb,left)});
 std::push_heap(heap.begin(),heap.end(),compare);
 heap.push_back({middle,parent.b,std::max(parent.lb,right)});
 std::push_heap(heap.begin(),heap.end(),compare);
 return true;
}
}
#endif
