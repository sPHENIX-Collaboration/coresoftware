#ifndef PHOOL_PHTIMER_H
#define PHOOL_PHTIMER_H

/*!
\file		PHTimer.h
\brief	 high precision timer
\author	Sean Kelly, Hugo Pereira
\version $Revision: 1.6 $
\date		$Date: 2013/01/07 09:27:01 $
*/

#include <unistd.h>
#include <exception>
#include <iostream>
#include <string>
#include <sstream>

#define rdtsc(low, high)       \
  __asm__ __volatile__("rdtsc" \
                       : "=a"(low), "=d"(high))

/*! \ingroup classes */
//! high precision timer
/*! high precision timer */
class PHTimer
{
 public:
  //! enum for timer state
  enum State
  {
    STOP = 0,
    RUN = 1
  };

  //! access timer state
  State get_state(void) const
  {
    return _state;
  }

  //! Construct with a name
  explicit PHTimer(const std::string& name = "Generic Timer")
    : _name(name)
    , _state(STOP)
    , _start_time(get_clock_counts())
    , _stop_time(get_clock_counts())
    , _accumulated_time(0)
    , _ncycle(0)
  {
    _stop_time._low++;
  }

  //! stops the counter
  void stop()
  {
    if (_state == STOP) return;
    _stop_time = get_clock_counts();
    _state = STOP;
    _ncycle++;
    _accumulated_time += elapsed();
  }

  //! Restart timer
  void restart()
  {
    _start_time = get_clock_counts();
    _state = RUN;
  }

  //! Dump elapsed time to provided ostream
  void print(std::ostream& os = std::cout) const
  {
    double elapse(elapsed());
    PRINT(os, "Timing for " + _name);
    os << "time (ms): " << elapse << std::endl;
    PRINT();
  }

  //! Dump statistics
  void print_stat(std::ostream& os = std::cout) const
  {
    //   PRINT(os, std::string("Stats for " + _name));
    if (_ncycle)
    {
      os << _name << ": accumulated time (ms):  " << _accumulated_time << std::endl;
      os << _name << ": per event time (ms):    " << _accumulated_time / _ncycle << std::endl;
    }
    else
    {
      os <<  _name << ": timer never started.\n";
    }
    PRINT(os, "**");
  }

  //! Set timer name
  void set_name(const std::string& name)
  {
    _name = name;
  }

  //! get timer name
  std::string get_name(void) const
  {
    return _name;
  }

  //! get cumulated time
  double get_accumulated_time(void) const
  {
    return _accumulated_time;
  }

  //! get number of cycles
  unsigned int get_ncycle(void) const
  {
    return _ncycle;
  }

  //! get averaged time/cycle
  double get_time_per_cycle(void) const
  {
    return _accumulated_time / _ncycle;
  }

  //! retrieve elapsed value since last restart (in ms)
  double elapsed(void) const
  {
    return 1000.0 * get_difference(
                        (_state == RUN) ? get_clock_counts() : _stop_time,
                        _start_time);
  }

  //! test PHTimer for a given amount of time (in ms)
  void test(double time, std::ostream& os = std::cout)
  {
    std::ostringstream tmp;
    tmp << "Test for " << _name << " - " << time << "ms";
    PRINT(os, tmp.str());
    restart();
    usleep((unsigned int) (time * 1000));
    print(os);
  }

  //! print a message (formated) to a stream
  static void PRINT(std::ostream& os = std::cout, const std::string& message = "");

 private:
  //! internal frequency read from cpu information file
  class Frequency
  {
   public:
    //! constructor
    Frequency()
    {
      try
      {
        set_cpu_freq();
      }
      catch (std::exception& e)
      {
        std::cerr << e.what() << std::endl;
      }
      _period = 1.0 / _frequency;
    }

    //! frequency accessor
    operator double() const
    {
      return _frequency;
    }

    //! period accessor
    double period() const
    {
      return _period;
    }

   private:
    //! pc frequency
    double _frequency;

    //! pc period
    double _period;

    //! read pc frequency from cpuinfo place
    void set_cpu_freq(const std::string &cpuinfopath = "/proc/cpuinfo");
  };

  //! used to store high precision time using two integers
  struct time_struct
  {
    //! constructor
    time_struct(void)
      : _low(0)
      , _high(0)
    {
    }

    //! low wheight bits cpu count
    unsigned long _low;

    //! high wheight bits cpu count
    unsigned long _high;
  };

  //! gets time from cpu clock counts
  static time_struct get_clock_counts(void)
  {
    time_struct t;
    rdtsc(t._low, t._high);
    return t;
  }

  //! returns difference between to time
  static double get_difference(
      const time_struct&,
      const time_struct&);

  //! static frequency object
  static Frequency _frequency;

  //! to stores 2^32
  static const double _twopower32;

  //! timer name
  std::string _name;

  //! timer state
  State _state;

  //! start time structure
  time_struct _start_time;

  //! stop time structure
  time_struct _stop_time;

  //! cumulated time
  double _accumulated_time;

  //! number of restart/stop cycles
  unsigned int _ncycle;
};

#endif
