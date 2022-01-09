// Tell emacs that this is a C++ source
//  -*- C++ -*-.
// $Id: $

/*!
 * \file PHG4OuterHcalField.h
 * \brief
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#ifndef G4DETECTORS_PHG4OUTERHCALFIELD_H
#define G4DETECTORS_PHG4OUTERHCALFIELD_H

#include <Geant4/G4MagneticField.hh>
#include <Geant4/G4Types.hh>  // for G4double, G4int

/*!
 * \brief PHG4OuterHcalField
 *
 * After burner to produce magnetic field in outer HCal on top of sPHENIX field map
 * In leading order, for the field within the plane of absorber plate,
 *  iron absorb almost all the flux, with field in the air-scintillator region
 * reduced to 1/relative_permeability_absorber of that in the iron.
 * The field strength perpendicular to the plate remain unchanged.
 *
 * relative_permeability_absorber = 1514, relative permeability for Steel 1006 @ B = 1.06T
 * http://www.fieldp.com/magneticproperties.html
 */
class PHG4OuterHcalField : public G4MagneticField
{
 public:
  PHG4OuterHcalField() = delete;

  PHG4OuterHcalField(bool isInIron, G4int steelPlates,
                     G4double scintiGap, G4double tiltAngle);

  ~PHG4OuterHcalField() override {}

  void
  GetFieldValue(const double Point[4], double *Bfield) const override;

  bool
  is_is_in_iron() const
  {
    return is_in_iron;
  }

  void
  set_is_in_iron(bool isInIron)
  {
    is_in_iron = isInIron;
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
  double relative_permeability_absorber = 1514.;
  // relative permeability for Steel 1006 @ B = 1.06T
  double relative_permeability_gap = 1;

  bool is_in_iron;
  G4int n_steel_plates;
  G4double scinti_gap;
  G4double tilt_angle;
};

#endif /* G4DETECTORS_PHG4OUTERHCALFIELD_H */
