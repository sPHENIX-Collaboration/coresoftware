// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2019 The HepMC collaboration (see AUTHORS for details)
//
#ifndef VALIDATION_CONTROL_H
#define VALIDATION_CONTROL_H

#ifdef HEPMC2
#include "HepMC/GenEvent.h"
#else
#include "HepMC3/GenEvent.h"
#include "HepMC3/Print.h"
#endif // ifdef HEPMC2

#include "ValidationTool.h"
#include "Timer.h"

#include <vector>
#include <string.h>
/// @class ValidationControl
/// @brief Runs multiple validation tools
class ValidationControl {
//
// Constructors
//
public:
    /** @brief Constructor */
    ValidationControl();
    /** @brief Destructor */
    ~ValidationControl();

//
// Functions
//
public:
    /** @brief Read file */
    void read_file(const std::string &filename);
    /** @brief New event */
    bool new_event();
    /** @brief Init function */
    void initialize();
    /** @brief Process event */
    void process(GenEvent &hepmc);
    /** @brief Finalize */
    void finalize();

//
// Accessors
//
public:
    /** @brief Toolchain */
    const std::vector<ValidationTool*>& toolchain() { return m_toolchain; }
    /** @brief Event limit */
    int   event_limit()                             { return m_events;    }
    /** @brief Set event limit */
    void  set_event_limit(int events)               { m_events = events;  }
    /** @brief N events to print*/
    void  print_events(int events)              { m_print_events          = events; }
    /** @brief N events to check momentum*/
    void  check_momentum_for_events(int events) { m_momentum_check_events = events; }

//
// Fields
//
private:
    std::vector<ValidationTool*> m_toolchain;  ///< Toolchain

    int    m_events;                     ///< events
    int    m_events_print_step;          ///< events print step
    int    m_momentum_check_events;      ///< mom check events
    double m_momentum_check_threshold;   ///< mom check threshold
    int    m_print_events;               ///< print events
    int    m_event_counter;               ///< counter of events
    int    m_status;                                     ///< status
    Timer  m_timer;              ///< Times

    bool m_has_input_source;     ///< Input source flag

    /** @brief parsing stutus */
    enum PARSING_STATUS {
        PARSING_OK,
        UNRECOGNIZED_COMMAND,
        UNRECOGNIZED_OPTION,
        UNRECOGNIZED_INPUT,
        UNRECOGNIZED_TOOL,
        UNAVAILABLE_TOOL,
        ADDITIONAL_INPUT,
        CANNOT_OPEN_FILE
    };
};

#endif
