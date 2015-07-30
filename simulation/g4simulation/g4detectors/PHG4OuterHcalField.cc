// $Id: $

/*!
 * \file PHG4OuterHcalField.cc
 * \brief
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#include "PHG4OuterHcalField.h"

#include <Geant4/G4FieldManager.hh>
#include <Geant4/G4TransportationManager.hh>
#include <Geant4/G4EquationOfMotion.hh>
#include <iostream>
using namespace std;

PHG4OuterHcalField::PHG4OuterHcalField() :

    /*double*/relative_permeability_absorber(1),
    /*double*/relative_permeability_gap(1),
    /*G4int*/n_steel_plates(0),
    /*G4double*/scinti_gap(0),
    /*G4double*/tilt_angle(0)
{
}

PHG4OuterHcalField::~PHG4OuterHcalField()
{
}

void
PHG4OuterHcalField::GetFieldValue(const double Point[4], double *Bfield) const
{
  G4FieldManager* field_manager =
      G4TransportationManager::GetTransportationManager()->GetFieldManager();

  if (!field_manager)
    {
      static bool once = true;

      if (once)
        {
          once = false;
          cout << "PHG4OuterHcalField::GetFieldValue"
              << " - Error! can not find field manager in G4TransportationManager"
              << endl;
        }
    }

  const G4Field* default_field = field_manager->GetDetectorField();

  if (default_field)
    {
      default_field->GetFieldValue(Point, Bfield);

      // scale the field value;
    }
  else
    {
      static bool once = true;

      if (once)
        {
          once = false;
          cout << "PHG4OuterHcalField::GetFieldValue"
              << " - Error! can not find detecor field in field manager!"
              << endl;
        }

    }
}
