#ifndef DCR_DEADLINE_H
#define DCR_DEADLINE_H

#include <chrono>
#include <cmath>

/// Cooperative wall-clock budget shared by a node and its nested oracles.
class DCRDeadline {
public:
 struct Expired {};
 explicit DCRDeadline( double seconds , const DCRDeadline * parent = nullptr ) :
  expires(std::chrono::steady_clock::now()+std::chrono::duration<double>(seconds)) {
  // Materialize the tighter absolute limit once: checking a deeply nested
  // oracle still needs only one clock read, independent of nesting depth.
  if(parent && parent->expires < expires) expires = parent->expires;
 }
 void check() const {
  if(std::chrono::steady_clock::now() >= expires) throw Expired{};
 }
private:
 std::chrono::time_point<std::chrono::steady_clock,
                         std::chrono::duration<double>> expires;
};

/// Non-owning deadline; unset for standalone legacy solvers by default.
class DCRDeadlineCheck {
public:
 void setDeadline( const DCRDeadline * value ) { deadline = value; }
 const DCRDeadline * getDeadline() const { return deadline; }
protected:
 void checkDeadline() const { if( deadline ) deadline->check(); }
private:
 const DCRDeadline * deadline = nullptr;
};

#endif
