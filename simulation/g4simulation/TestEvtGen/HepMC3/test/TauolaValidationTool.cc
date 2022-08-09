// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2019 The HepMC collaboration (see AUTHORS for details)
//
#include "TauolaValidationTool.h"

void TauolaValidationTool::initialize() {
    Tauolapp::Tauola::setSameParticleDecayMode(4);
    Tauolapp::Tauola::setOppositeParticleDecayMode(4);
    Tauolapp::Tauola::initialize();
}

int TauolaValidationTool::process(GenEvent &hepmc) {

    HEPMC2CODE( Tauolapp::TauolaHepMCEvent  t_event(&hepmc); )
    HEPMC3CODE( Tauolapp::TauolaHepMC3Event t_event(&hepmc); )

    //t_event.undecayTaus();
    t_event.decayTaus();

    return 0;
}

void TauolaValidationTool::finalize() {
    Tauolapp::Tauola::summary();
    Tauolapp::Log::Summary();
}
