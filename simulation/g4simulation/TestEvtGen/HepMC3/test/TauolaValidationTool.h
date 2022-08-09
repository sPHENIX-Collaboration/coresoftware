// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2019 The HepMC collaboration (see AUTHORS for details)
//
#ifndef TAUOLA_VALIDATION_TOOL_H
#define TAUOLA_VALIDATION_TOOL_H

#ifdef HEPMC2
#include "Tauola/TauolaHepMCEvent.h"
#include "HepMC/GenEvent.h"
#else
#include "Tauola/TauolaHepMC3Event.h"
#include "HepMC3/GenEvent.h"
#endif // ifdef HEPMC2

#include "ValidationTool.h"
#include "Timer.h"

#include "Tauola/Tauola.h"
#include "Tauola/Log.h"
/// @class TauolaValidationTool
/// @brief Interface for validatio to Tauola
class TauolaValidationTool : public ValidationTool {
public:
    TauolaValidationTool():m_timer("Tauola++ processing time") {}

public:
    const std::string name()     { return "Tauola++"; }
    bool   tool_modifies_event() { return true;       }
    Timer* timer()               { return &m_timer;   }

    void initialize();
    int  process(GenEvent &hepmc);
    void finalize();

private:
    Timer m_timer; ///< Timer
};

#endif
