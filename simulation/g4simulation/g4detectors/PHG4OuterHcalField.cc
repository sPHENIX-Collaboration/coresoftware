// $Id: $

/*!
 * \file PHG4OuterHcalField.cc
 * \brief
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#include "PHG4OuterHcalField.h"

#include <Geant4/G4Vector3D.hh>
#include <Geant4/G4Field.hh>                   // for G4Field
#include <Geant4/G4FieldManager.hh>
#include <Geant4/G4TransportationManager.hh>
#include <Geant4/G4PhysicalConstants.hh>
#include <Geant4/G4Types.hh>                  // for G4double, G4int

#include <cassert>                            // for assert
#include <cmath>                              // for atan2, cos, sin, sqrt
#include <iostream>

using namespace std;

PHG4OuterHcalField::PHG4OuterHcalField(bool isInIron, G4int steelPlates,
    G4double scintiGap, G4double tiltAngle) :
    /*double*/relative_permeability_absorber(1514), // relative permeability for Steel 1006 @ B = 1.06T
    /*double*/relative_permeability_gap(1),
    /*bool*/is_in_iron(isInIron),
    /*G4int*/n_steel_plates(steelPlates),
    /*G4double*/scinti_gap(scintiGap),
    /*G4double*/tilt_angle(tiltAngle)
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

      // scale_factor for field component along the plate surface
      double x = Point[0];
      double y = Point[1];
//      double z = Point[2];

      assert(cos(tilt_angle)>0);
      assert(n_steel_plates>0);

      // input 2D magnetic field vector
      const G4Vector3D B0(Bfield[0], Bfield[1], Bfield[2]);
      const G4Vector3D B0XY(Bfield[0], Bfield[1], 0);
      const G4Vector3D B0Z(0, 0, Bfield[2]);

      const double R = sqrt(x * x + y * y);
      const double layer_RdPhi = R * twopi / n_steel_plates;
      const double layer_width = layer_RdPhi * cos(tilt_angle);
      const double gap_width = scinti_gap;

      assert(gap_width<layer_width);

      // sign definition of tilt_angle is rotation around the -z axis
      const G4Vector3D absorber_dir(cos(atan2(y, x) - tilt_angle),
          sin(atan2(y, x) - tilt_angle), 0);
      const G4Vector3D radial_dir(cos(atan2(y, x)), sin(atan2(y, x)), 0);

      const double radial_flux_per_layer = layer_RdPhi
          * (B0XY.dot(radial_dir));
      double B_XY_mag = radial_flux_per_layer
          / (relative_permeability_absorber * (layer_width - gap_width)
              + relative_permeability_gap * gap_width);
      B_XY_mag *=
          is_in_iron ?
              relative_permeability_absorber : relative_permeability_gap;

      const G4Vector3D B_New_XY = B_XY_mag * absorber_dir;

      const double z_flux_per_layer = layer_width * B0Z.z();
      double B_Z_mag = z_flux_per_layer
          / (relative_permeability_absorber * (layer_width - gap_width)
              + relative_permeability_gap * gap_width);
      B_Z_mag *=
          is_in_iron ?
              relative_permeability_absorber : relative_permeability_gap;
      const G4Vector3D B_New_Z(0, 0, B_Z_mag);

      // scale B_T component
      G4Vector3D B_New = B_New_Z + B_New_XY;
      Bfield[0] = B_New.x();
      Bfield[1] = B_New.y();
      Bfield[2] = B_New.z();
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
