// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2019 The HepMC collaboration (see AUTHORS for details)
//
#ifndef OUTPUT_VALIDATION_TOOL_H
#define OUTPUT_VALIDATION_TOOL_H

#ifdef HEPMC2
#include "HepMC/GenEvent.h"
#include "HepMC/IO_GenEvent.h"
#else
#include "HepMC3/GenEvent.h"
#include "HepMC3/WriterAsciiHepMC2.h"
#endif

#include "ValidationTool.h"
#include "Timer.h"

/// @class OutputValidationTool
/// @brief Interface for validatio to Pythia
class OutputValidationTool : public ValidationTool {
public:
    OutputValidationTool( const std::string &filename ); ///< Constructor

    const std::string name()      { return "OUTPUT"; }
    const std::string long_name() { return name() + " config file: " + m_filename; }

    bool   tool_modifies_event() { return false;      }
    Timer* timer()               { return &m_timer;  }

    void initialize();
    int  process(GenEvent &hepmc);
    void finalize();

private:
    std::string     m_filename; ///< Used file
    Timer           m_timer;  ///< Timer
    HEPMC2CODE( IO_GenEvent * m_file; )
    HEPMC3CODE( WriterAsciiHepMC2 * m_file; )
};

#endif
