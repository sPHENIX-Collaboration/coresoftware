// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2019 The HepMC collaboration (see AUTHORS for details)
//
#include "PythiaValidationTool.h"

PythiaValidationTool::PythiaValidationTool( const std::string &filename ):m_filename(filename),m_timer("pythia8 conversion time") {
    m_pythia.readFile(m_filename);
}

void PythiaValidationTool::initialize() {
    m_pythia.init();
}

int PythiaValidationTool::process(GenEvent &hepmc) {
    if( !m_pythia.next() ) return 1;

    // Exclude generation time
    m_timer.start();

    m_tohepmc.fill_next_event(m_pythia.event, &hepmc, -1, &m_pythia.info);
    return 0;
}

void PythiaValidationTool::finalize() {
    /* The condition below is true at least for 8.209+. 8.209- will probably fail */
#if defined(PYTHIA_VERSION_INTEGER) || defined (PYTHIA_VERSION)
    m_pythia.stat();
#else
    m_pythia.statistics();
#endif
}
