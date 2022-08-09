// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2019 The HepMC collaboration (see AUTHORS for details)
//
#ifndef MCTESTER_TOOL_H
#define MCTESTER_TOOL_H

#ifdef  HEPMC2
#include "HepMC/GenEvent.h"
#include "HepMCEvent.H"
#else
#include "HepMC3/GenEvent.h"
#include "HepMC3Event.h"
#endif // ifdef HEPMC2

#include "ValidationTool.h"

#include "Setup.H"
#include "Generate.h"
/// @class McTesterValidationTool
/// @brief Interface to MCTester
class McTesterValidationTool : public ValidationTool {
public:
    const std::string name()    { return "MC-TESTER"; }
    bool  tool_modifies_event() { return false; }

    void initialize();
    int  process(GenEvent &hepmc);
    void finalize();
};

#endif
