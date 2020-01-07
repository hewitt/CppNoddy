/// \file Timer.h
/// A spec for the CppNoddy Timer object. This is used
/// to provide course timing of test codes if -DTIME compilation
/// flag is set.

#ifndef TIMER_H
#define TIMER_H

#include <ctime>
#include <string>
#include <time.h>
#include <sys/time.h>

namespace CppNoddy {

  /// A simple CPU-clock-tick timer for timing
  /// metods.

  class Timer {

   public:
    // CONSTRUCTORS
    Timer() : STOPPED(true), COUNTER(0),
	      T_START(clock()), DELTA_T_STORE(0), m_deltaWall(0.0) {
      HEADER = "";
      m_wallStart = get_wall_time();
    }

    Timer(std::string name) : STOPPED(true), COUNTER(0),
			      T_START(clock()), DELTA_T_STORE(0), m_deltaWall(0.0) {
      HEADER = name;
      m_wallStart = get_wall_time();
    }

    // wall timer
    double get_wall_time(){
      struct timeval time;
      if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
      }
      return (double)time.tv_sec + (double)time.tv_usec * .000001;      
    }

    // CLASS METHODS
    
    /// Start the timer & reset stored time to zero.
    void start();

    /// Stop the clock & add the current time interval to the
    /// previously stored values ready for printing.
    void stop();

    /// Pause the clock & add the time interval to the
    /// stored cumulative time.
    void reset();

    /// Write a string to cout stating the time taken.
    void print() const;

    /// Increment an internal discrete counter.
    int& counter();

    /// \return The average time taken per increment of the
    /// internal discrete counter.
    double time_per_count() const;

    /// Return the time.
    /// \return The time (in ms).
    double get_time() const;

   private:
    // ATTRIBUTES
    bool STOPPED;
    // a counter - useful for working out time per counter click
    int COUNTER;
    // the start time, stop time and time elapsed so far
    clock_t T_START, DELTA_T_STORE;
    // wall time
    double m_wallStart, m_deltaWall;
    // a string header for the object
    std::string HEADER;
  };

}

#endif // TIMER_H
