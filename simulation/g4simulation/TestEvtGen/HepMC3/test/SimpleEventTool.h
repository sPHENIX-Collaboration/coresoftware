// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2019 The HepMC collaboration (see AUTHORS for details)
//
#ifndef SIMPLE_EVENT_TOOL_H
#define SIMPLE_EVENT_TOOL_H

#ifdef HEPMC2
#include "HepMC/IO_GenEvent.h"
#else
#include "HepMC3/GenVertex.h"
#include "HepMC3/GenParticle.h"
#endif // ifdef HEPMC2

#include "ValidationTool.h"
/// @class SimpleEventTool
/// @brief Simple validation
class SimpleEventTool : public ValidationTool {
public:
    SimpleEventTool() {}

    const std::string name()    { return "SimpleEventTool"; }
    bool  tool_modifies_event() { return true; }

    void initialize() {}
    int  process(GenEvent &hepmc);
    void finalize()   {}

};

#endif
