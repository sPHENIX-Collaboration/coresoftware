// $Id: $

/*!
 * \file PHG4OuterHcalFieldSetup.cc
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#include "PHG4OuterHcalFieldSetup.h"
#include "PHG4OuterHcalField.h"

#include <Geant4/G4ChordFinder.hh>
#include <Geant4/G4ClassicalRK4.hh>
#include <Geant4/G4FieldManager.hh>
#include <Geant4/G4MagIntegratorDriver.hh>
#include <Geant4/G4MagIntegratorStepper.hh>
#include <Geant4/G4Mag_UsualEqRhs.hh>
#include <Geant4/G4MagneticField.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4Types.hh>

PHG4OuterHcalFieldSetup::PHG4OuterHcalFieldSetup(G4int steelPlates,
                                                 G4double scintiGap, G4double tiltAngle)
  : fMinStep(0.005 * mm)
  , n_steel_plates(steelPlates) /*G4int*/
  , scinti_gap(scintiGap)       /*G4double*/
  , tilt_angle(tiltAngle)       /*G4double*/
{
  G4int nvar = 8;

  {
    fEMfieldIron = new PHG4OuterHcalField(true, n_steel_plates, scinti_gap,
                                          tilt_angle);

    fEquationIron = new G4Mag_UsualEqRhs(fEMfieldIron);

    fStepperIron = new G4ClassicalRK4(fEquationIron, nvar);

    fChordFinderIron = new G4ChordFinder(
        new G4MagInt_Driver(fMinStep, fStepperIron,
                            fStepperIron->GetNumberOfVariables()));

    fFieldManagerIron = new G4FieldManager();
    fFieldManagerIron->SetDetectorField(fEMfieldIron);
    fFieldManagerIron->SetChordFinder(fChordFinderIron);
  }

  {
    fEMfieldGap = new PHG4OuterHcalField(false, n_steel_plates, scinti_gap,
                                         tilt_angle);

    fEquationGap = new G4Mag_UsualEqRhs(fEMfieldGap);

    fStepperGap = new G4ClassicalRK4(fEquationGap, nvar);

    fChordFinderGap = new G4ChordFinder(
        new G4MagInt_Driver(fMinStep, fStepperGap,
                            fStepperGap->GetNumberOfVariables()));

    fFieldManagerGap = new G4FieldManager();
    fFieldManagerGap->SetDetectorField(fEMfieldGap);
    fFieldManagerGap->SetChordFinder(fChordFinderGap);
  }
}
