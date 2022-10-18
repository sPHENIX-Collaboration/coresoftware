// Tell emacs that this is a C++ source
//  -*- C++ -*-.
// $Id: $

/*!
 * \file PHG4MagneticField.h
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#ifndef G4MAIN_PHG4MAGNETICFIELD_H
#define G4MAIN_PHG4MAGNETICFIELD_H

#include <Geant4/G4MagneticField.hh>

class PHField;

/*!
 * \brief PHG4MagneticField interfaces with Geant4
 */
class PHG4MagneticField : public G4MagneticField
{
 public:
  PHG4MagneticField(const PHField* field);
  ~PHG4MagneticField() override;

  const PHField* get_field() const
  {
    return field_;
  }

  void set_field(const PHField* field)
  {
    field_ = field;
  }

  void GetFieldValue( const double Point[4],    double *Bfield ) const override;

 private:
  const PHField* field_;
};

#endif /* SIMULATION_CORESOFTWARE_SIMULATION_G4SIMULATION_G4MAIN_PHG4MAGNETICFIELD_H_ */
