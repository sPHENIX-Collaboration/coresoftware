// $Id: $

/*!
 * \file PHG4OHCalFieldSetup.cc
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#include "PHG4OHCalFieldSetup.h"

#include <g4main/PHG4MagneticField.h>

#include <phfield/PHFieldUtility.h>
#include <phfield/PHFieldConfigv1.h>

#include <Geant4/G4ChordFinder.hh>
#include <Geant4/G4ClassicalRK4.hh>
#include <Geant4/G4FieldManager.hh>
#include <Geant4/G4MagIntegratorDriver.hh>
#include <Geant4/G4MagIntegratorStepper.hh>
#include <Geant4/G4Mag_UsualEqRhs.hh>
#include <Geant4/G4MagneticField.hh>
#include <Geant4/G4SystemOfUnits.hh>

#include <cassert>
#include <iostream>

PHG4OHCalFieldSetup::PHG4OHCalFieldSetup(const std::string &iron_fieldmap_path, const double scale, const double inner_radius, const double outer_radius, const double size_z)
  : fMinStep(0.005 * mm)
{
  static const G4int nvar = 8;

  // the new HCal expect 3D magnetic field
  PHFieldConfigv1 field_config (PHFieldConfig::Field3DCartesian, iron_fieldmap_path, scale);

  fEMfieldIron = new PHG4MagneticField(PHFieldUtility::BuildFieldMap(&field_config, inner_radius, outer_radius, size_z));
  assert(fEMfieldIron);

  fEquationIron = new G4Mag_UsualEqRhs(fEMfieldIron);

  fStepperIron = new G4ClassicalRK4(fEquationIron, nvar);

  fChordFinderIron = new G4ChordFinder(
      new G4MagInt_Driver(fMinStep, fStepperIron,
                          fStepperIron->GetNumberOfVariables()));

  fFieldManagerIron = new G4FieldManager();
  fFieldManagerIron->SetDetectorField(fEMfieldIron);
  fFieldManagerIron->SetChordFinder(fChordFinderIron);
}

PHG4OHCalFieldSetup::~PHG4OHCalFieldSetup()
{
  delete fEMfieldIron;
}
