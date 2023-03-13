// Tell emacs that this is a C++ source
//  -*- C++ -*-.
// $Id: $

/*!
 * \file PHG4OuterHcalFieldSetup.h
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#ifndef G4DETECTORS_PHG4OUTERHCALFIELDSETUP_H
#define G4DETECTORS_PHG4OUTERHCALFIELDSETUP_H

#include <Geant4/G4Types.hh>  // for G4double, G4int

#include <cmath>

class G4ChordFinder;
class G4FieldManager;
class G4Mag_UsualEqRhs;
class G4MagIntegratorStepper;
class G4MagneticField;

/*!
 * \brief PHG4OuterHcalFieldSetup following Geant4 example F03FieldSetup
 */
class PHG4OuterHcalFieldSetup
{
 public:
  PHG4OuterHcalFieldSetup(G4int steelPlates, G4double scintiGap,
                          G4double tiltAngle);
  //! delete copy ctor and assignment opertor (cppcheck)
  explicit PHG4OuterHcalFieldSetup(const PHG4OuterHcalFieldSetup&) = delete;
  PHG4OuterHcalFieldSetup& operator=(const PHG4OuterHcalFieldSetup&) = delete;

  virtual ~PHG4OuterHcalFieldSetup() {}

  G4FieldManager*
  get_Field_Manager_Gap() const
  {
    return fFieldManagerGap;
  }

  void
  set_Field_Manager_Gap(G4FieldManager* fieldManagerGap)
  {
    fFieldManagerGap = fieldManagerGap;
  }

  G4FieldManager*
  get_Field_Manager_Iron() const
  {
    return fFieldManagerIron;
  }

  void
  set_Field_Manager_Iron(G4FieldManager* fieldManagerIron)
  {
    fFieldManagerIron = fieldManagerIron;
  }

  G4double
  get_Min_Step() const
  {
    return fMinStep;
  }

  void
  set_Min_Step(G4double minStep)
  {
    fMinStep = minStep;
  }

  G4int
  get_steel_plates() const
  {
    return n_steel_plates;
  }

  void
  set_steel_plates(G4int steelPlates)
  {
    n_steel_plates = steelPlates;
  }

  G4double
  get_scinti_gap() const
  {
    return scinti_gap;
  }

  void
  set_scinti_gap(G4double scintiGap)
  {
    scinti_gap = scintiGap;
  }

  G4double
  get_tilt_angle() const
  {
    return tilt_angle;
  }

  void
  set_tilt_angle(G4double tiltAngle)
  {
    tilt_angle = tiltAngle;
  }

 private:
  G4FieldManager* fFieldManagerIron = nullptr;
  G4FieldManager* fFieldManagerGap = nullptr;
  G4Mag_UsualEqRhs* fEquationIron = nullptr;
  G4Mag_UsualEqRhs* fEquationGap = nullptr;
  G4ChordFinder* fChordFinderIron = nullptr;
  G4ChordFinder* fChordFinderGap = nullptr;
  G4MagneticField* fEMfieldIron = nullptr;
  G4MagneticField* fEMfieldGap = nullptr;
  G4MagIntegratorStepper* fStepperIron = nullptr;
  G4MagIntegratorStepper* fStepperGap = nullptr;

  G4double fMinStep = NAN;

  G4int n_steel_plates = -1;
  G4double scinti_gap = NAN;
  G4double tilt_angle = NAN;
};

#endif /* PHG4OUTERHCALFIELDSETUP_H_ */
