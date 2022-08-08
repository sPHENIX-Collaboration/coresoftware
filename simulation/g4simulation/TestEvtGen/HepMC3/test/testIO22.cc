// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2022 The HepMC collaboration (see AUTHORS for details)
//
// -- Purpose: Test that an event can be written out to a protobuf file
//
// -- Note: The output of this test was used to generate inputs to tests 20 and
// 21
//

#include "HepMC3/Writerprotobuf.h"

#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenRunInfo.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/Print.h"

#include <cassert>

int main() {

  int const kISStatus = 4;
  int const kFSStatus = 1;
  HepMC3::GenParticlePtr ISNeutron = std::make_shared<HepMC3::GenParticle>(
      HepMC3::FourVector{1.5255172492130473e+02, 8.9392830847276528e+01,
                         6.4870597568257821e+01, 9.5825554558124941e+02},
      2112, kISStatus);
  HepMC3::GenParticlePtr ISnumu = std::make_shared<HepMC3::GenParticle>(
      HepMC3::FourVector{0.0000000000000000e+00, 0.0000000000000000e+00,
                         1.5000000000000000e+03, 1.5000000000000000e+03},
      14, kISStatus);
  HepMC3::GenParticlePtr FSmuon = std::make_shared<HepMC3::GenParticle>(
      HepMC3::FourVector{-6.8928697531845643e+01, 4.8219068401438176e+02,
                         1.2406574501351240e+03, 1.3370316161682497e+03},
      13, kFSStatus);
  HepMC3::GenParticlePtr FSProton = std::make_shared<HepMC3::GenParticle>(
      HepMC3::FourVector{2.2148042245314980e+02, -3.9279785316710411e+02,
                         3.2421314743313258e+02, 1.0903266675337304e+03},
      2212, kFSStatus);

  // Set masses
  ISNeutron->set_generated_mass(9.3956499999999994e+02);
  ISnumu->set_generated_mass(0.0000000000000000e+00);
  FSmuon->set_generated_mass(1.0565800000000023e+02);
  FSProton->set_generated_mass(9.3827200000000005e+02);

  // Make vertex
  HepMC3::GenVertexPtr vertex = std::make_shared<HepMC3::GenVertex>(
      HepMC3::FourVector{1E1, 2E2, 3E3, 4E4});
  vertex->add_particle_in(ISNeutron);
  vertex->add_particle_in(ISnumu);
  vertex->add_particle_out(FSmuon);
  vertex->add_particle_out(FSProton);

  HepMC3::GenEvent evt;

  evt.add_vertex(vertex);
  evt.weights() = std::vector<double>{1.23456789, 9.87654321};
  evt.set_event_number(1337);
  evt.set_units(HepMC3::Units::MEV, HepMC3::Units::CM);
  evt.add_attribute("HardScatterMode", std::make_shared<HepMC3::IntAttribute>(1));

  std::shared_ptr<HepMC3::GenRunInfo> gri =
      std::make_shared<HepMC3::GenRunInfo>();
  gri->tools().emplace_back(HepMC3::GenRunInfo::ToolInfo{
      "NuDum", "0.99.0", "A dummy neutrino event generator"});
  evt.set_run_info(gri);

  auto writer =
      std::make_shared<HepMC3::Writerprotobuf>("outputIO24.proto", gri);
  writer->write_event(evt);
  
  HepMC3::Print::listing(*evt.run_info());
  HepMC3::Print::listing(evt);
  writer->close();
  return 0;
}
