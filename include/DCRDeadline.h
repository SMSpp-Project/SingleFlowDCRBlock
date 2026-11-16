#ifndef DCR_DEADLINE_H
#define DCR_DEADLINE_H

#include <chrono>
#include <cmath>

/// Cooperative wall-clock budget shared by a node and its nested oracles.
class DCRDeadline {
public:
 struct Expired {};
 explicit DCRDeadline( double seconds ) : seconds( seconds ),
  start( std::chrono::steady_clock::now() ) {}
 void check() const {
  if( std::isfinite( seconds ) &&
      std::chrono::duration<double>( std::chrono::steady_clock::now() - start ).count()
       >= seconds )
   throw Expired{};
 }
private:
 double seconds;
 std::chrono::steady_clock::time_point start;
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
