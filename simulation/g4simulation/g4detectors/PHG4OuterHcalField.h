// $Id: $

/*!
 * \file PHG4OuterHcalField.h
 * \brief
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#ifndef PHG4OUTERHCALFIELD_H_
#define PHG4OUTERHCALFIELD_H_

#include <Geant4/G4MagneticField.hh>
#include <Geant4/globals.hh>
#include <Geant4/G4ios.hh>
/*!
 * \brief PHG4OuterHcalField
 */
class PHG4OuterHcalField : public G4MagneticField
{
public:
  PHG4OuterHcalField();
  virtual
  ~PHG4OuterHcalField();

  void
  GetFieldValue(const double Point[4], double *Bfield) const;

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

  double
  get_relative_permeability_absorber() const
  {
    return relative_permeability_absorber;
  }

  void
  set_relative_permeability_absorber(double relativePermeabilityAbsorber)
  {
    relative_permeability_absorber = relativePermeabilityAbsorber;
  }

  double
  get_relative_permeability_gap() const
  {
    return relative_permeability_gap;
  }

  void
  set_relative_permeability_gap(double relativePermeabilityGap)
  {
    relative_permeability_gap = relativePermeabilityGap;
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

  double relative_permeability_absorber;
  double relative_permeability_gap;

  G4int n_steel_plates;
  G4double scinti_gap;
  G4double tilt_angle;
};

#endif /* PHG4OUTERHCALFIELD_H_ */
