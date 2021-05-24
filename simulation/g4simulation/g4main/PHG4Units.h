// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4UNITS_H
#define G4MAIN_PHG4UNITS_H

// We need these units for initializers in the ctors of 
// so we cannot use G4UnitDefinition for these
#include <Geant4/G4SystemOfUnits.hh>

static const double inch = 2.54*cm;

#endif
