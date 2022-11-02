// Tell emacs that this is a C++ source
//  -*- C++ -*-.
// $Id: $

/*!
 * \file PHG4OHCalFieldSetup.h
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#ifndef G4OHCAL_PHG4OHCALFIELDSETUP_H
#define G4OHCAL_PHG4OHCALFIELDSETUP_H

#include <Geant4/G4Types.hh>  // for G4double, G4int

#include <string>

class G4ChordFinder;
class G4FieldManager;
class G4Mag_UsualEqRhs;
class G4MagIntegratorStepper;
class G4MagneticField;

/*!
 * \brief PHG4OHCalFieldSetup following Geant4 example F03FieldSetup
 */
class PHG4OHCalFieldSetup
{
 public:
  PHG4OHCalFieldSetup(const std::string & iron_fieldmap_path, const double scale = 1., const double inner_radius = 0., const double outer_radius = 1.e10, const double size_z = 1.e10);

  // delete copy ctor and assignment opertor (cppcheck)
  explicit PHG4OHCalFieldSetup(const PHG4OHCalFieldSetup&) = delete;
  PHG4OHCalFieldSetup& operator=(const PHG4OHCalFieldSetup&) = delete;

  virtual ~PHG4OHCalFieldSetup();

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

 private:
  G4FieldManager* fFieldManagerIron = nullptr;
  G4Mag_UsualEqRhs* fEquationIron = nullptr;
  G4ChordFinder* fChordFinderIron = nullptr;
  G4MagneticField* fEMfieldIron = nullptr;
  G4MagIntegratorStepper* fStepperIron = nullptr;
  G4double fMinStep = NAN;
};

#endif /* G4OHCAL_PHG4OHCALFIELDSETUP_H */
