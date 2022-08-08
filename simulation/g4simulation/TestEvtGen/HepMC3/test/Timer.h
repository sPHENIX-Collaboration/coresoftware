// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2021 The HepMC collaboration (see AUTHORS for details)
//
#ifndef BENCHMARK_TIMER_H
#define BENCHMARK_TIMER_H

#include <iostream>

#if defined(_MSC_VER)|| defined(_WIN32)
/**
 *  @class Timer
 *  @brief Used to benchmark MC generators
 *
 */
class Timer {
public:
    /** Default constructor */
    Timer(const char* name):m_name(name) { reset(); }
    void start() {}
    int elapsed_time() {return 0;}
    int total_time() {return 0;}
    void stop() {}
    void reset() {}
    void print() { printf("<sys/times.h> header is not present in MS Visual Studio. Dummy implementation of Timer class is used.\n");}
private:
    const char *m_name;

};
#else
#include <sys/times.h>
/**
 *  @class Timer
 *  @brief Used to benchmark MC generators
 *
 */
class Timer {
public:
    /** Default constructor */
    Timer(const char* name):m_name(name) { reset(); }

    /** Start or restart the timer */
    void start() {
        times(&m_start);
    }

    /** Get time elapsed since timer started */
    int elapsed_time() {
        times(&m_stop);
        return 10*(m_stop.tms_utime - m_start.tms_utime + m_stop.tms_stime - m_start.tms_stime);
    }

    /** Get total time counted by the timer */
    int total_time() {
        return 10*(m_stored.tms_utime + m_stored.tms_stime);
    }

    /** Save end time and aggregate build-in clock */
    void stop() {
        // Do nothing if timer has not been started
        if(m_start.tms_utime == 0) return;

        times(&m_stop);

        m_stored.tms_utime += m_stop.tms_utime - m_start.tms_utime;
        m_stored.tms_stime += m_stop.tms_stime - m_start.tms_stime;

        m_start.tms_utime = 0;
        m_start.tms_stime = 0;
    }

    /** Reset the clock */
    void reset() {
        m_start.tms_utime = 0;
        m_start.tms_stime = 0;
        m_stored.tms_utime = 0;
        m_stored.tms_stime = 0;
    }

    /** Print time elapsed */
    void print() {
        std::cout << m_name << ":" << std::endl;
        std::cout << "  user:   " << m_stored.tms_utime*10 << " ms" << std::endl;
        std::cout << "  system: " << m_stored.tms_stime*10 << " ms" << std::endl;
    }

private:
    const char *m_name; ///< Name of the object
    struct tms  m_start; ///< Start
    struct tms  m_stop; ///< Stop
    struct tms  m_stored; ///< Stored
};
#endif
#endif
