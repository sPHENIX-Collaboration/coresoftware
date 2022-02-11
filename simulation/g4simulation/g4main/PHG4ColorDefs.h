// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4COLORDEFS_H
#define G4MAIN_PHG4COLORDEFS_H

#include <Geant4/G4VisAttributes.hh>

namespace PHG4TpcColorDefs
{
  static G4Colour tpc_cu_color = G4Colour::Yellow();
  static G4Colour tpc_fr4_color = G4Colour::Grey();
  static G4Colour tpc_gas_color = G4Colour::Magenta();
  static G4Colour tpc_honeycomb_color = G4Colour::White();
  static G4Colour tpc_kapton_color = G4Colour::Green();
  static G4Colour tpc_pcb_color = G4Colour::Blue();
}

#endif
