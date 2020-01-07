/// \file Timer.cpp
/// Implementation of the CppNoddy Timer object.

#include <ctime>
#include <iostream>

#include <Timer.h>
#include <Exceptions.h>

namespace CppNoddy {

  void Timer::start() {
    T_START = clock();
    STOPPED = false;
  }

  void Timer::stop() {
    if(!STOPPED) {
      // we add current time interval to that stored
      DELTA_T_STORE += clock() - T_START;
      STOPPED = true;
      m_deltaWall = get_wall_time() - m_wallStart;
    }
  }

  void Timer::reset() {
    stop();
    COUNTER = 0;
    STOPPED = true;
    DELTA_T_STORE = 0;
    m_deltaWall = 0.0;
  }

  double Timer::get_time() const {
    if(STOPPED) {
      // if stopped or paused return stored intervals
      return 1.e3 * double(DELTA_T_STORE) / CLOCKS_PER_SEC;
    } else {
      // running clock
      return 1.e3 * (double(DELTA_T_STORE) + double(clock()) - double(T_START)) / CLOCKS_PER_SEC;
    }
  }

  int& Timer::counter() {
    return COUNTER;
  }

  double Timer::time_per_count() const {
    if(STOPPED) {
      return 1.e3 * double(DELTA_T_STORE) / CLOCKS_PER_SEC / COUNTER;
    } else {
      std::string problem = HEADER;
      problem += "\n The Timer object can only return time_per_count after being stopped.\n";
      throw ExceptionRuntime(problem);
    }
    return 0; // dummy
  }

  void Timer::print() const {
    std::cout.precision(4);
    std::cout << "  " << HEADER << "\n";
    const double elapsed_time_in_ms(1.e3 * double(DELTA_T_STORE) / CLOCKS_PER_SEC);
    if(elapsed_time_in_ms > 1000) {
      std::cout << "  TOTAL CPU time taken  = " << elapsed_time_in_ms / 1000. << " s\n";
      std::cout << "  TOTAL wall time taken = " << m_deltaWall << " s\n";
    } else {
      std::cout << "  TOTAL CPU time taken = " << elapsed_time_in_ms << " ms\n";
    }
    if(COUNTER != 0) {
      std::cout << "  Number of loops during this time = " << COUNTER << "\n";
      std::cout << "  Throughput = " << 1.e3 * COUNTER / elapsed_time_in_ms << " runs/s \n";
    }
  }

} // end CppNoddy namespace
