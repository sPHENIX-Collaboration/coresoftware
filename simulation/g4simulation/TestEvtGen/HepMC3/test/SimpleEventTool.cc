// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2019 The HepMC collaboration (see AUTHORS for details)
//
#include "SimpleEventTool.h"

int SimpleEventTool::process(GenEvent &hepmc) {
    // Create some four vectors for the electrons
    double ele_mass_sqr = 0.000511*0.000511;
    FourVector momentum_e1;
    FourVector momentum_e2;

    momentum_e1.setPz( -2);
    momentum_e2.setPz(3.5);

    momentum_e1.setE(std::sqrt(momentum_e1.pz()*momentum_e1.pz() + ele_mass_sqr));
    momentum_e2.setE(std::sqrt(momentum_e2.pz()*momentum_e2.pz() + ele_mass_sqr));

    FourVector momentum_tau1(+1.38605041e+00,+1.38605041e+00,+7.50000000e-01,+2.75000005e+00);
    FourVector momentum_tau2(-1.38605041e+00,-1.38605041e+00,+7.50000000e-01,+2.75000005e+00);

    // Make particles
    HEPMC2CODE(
        GenParticle *e1     = new GenParticle( momentum_e1,   -11, 2 );
        GenParticle *e2     = new GenParticle( momentum_e2,    11, 2 );
        GenParticle *tau1   = new GenParticle( momentum_tau1, -15, 1 );
        GenParticle *tau2   = new GenParticle( momentum_tau2,  15, 1 );
        GenVertex   *vertex = new GenVertex();
    )
    HEPMC3CODE(
        // Although the code for HepMC2 would work (thanks to backward compatibility)
        // we don't want to use deprecated functions
        GenParticlePtr e1     = std::make_shared<GenParticle>( momentum_e1,   -11, 2 );
        GenParticlePtr e2     = std::make_shared<GenParticle>( momentum_e2,    11, 2 );
        GenParticlePtr tau1   = std::make_shared<GenParticle>( momentum_tau1, -15, 1 );
        GenParticlePtr tau2   = std::make_shared<GenParticle>( momentum_tau2,  15, 1 );
        GenVertexPtr   vertex = std::make_shared<GenVertex>();
    )

    // Set masses
    e1->  set_generated_mass(0.000511);
    e2->  set_generated_mass(0.000511);
    tau1->set_generated_mass(1.777);
    tau2->set_generated_mass(1.777);

    // Make vertex
    vertex->add_particle_in(e1);
    vertex->add_particle_in(e2);
    vertex->add_particle_out(tau1);
    vertex->add_particle_out(tau2);

    hepmc.add_vertex(vertex);
    return 0;
}
