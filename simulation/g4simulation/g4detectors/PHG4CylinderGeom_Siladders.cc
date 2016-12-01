#include "PHG4CylinderGeom_Siladders.h"
#include <cmath>

#include "g4main/PHG4Detector.h"

#include <G4RotationMatrix.hh>
#include <G4Transform3D.hh>

#include <boost/format.hpp>
#include <TMath.h>
#include <TVector3.h>
#include <TRotation.h>

ClassImp(PHG4CylinderGeom_Siladders)

PHG4CylinderGeom_Siladders::~PHG4CylinderGeom_Siladders()
{}

void PHG4CylinderGeom_Siladders::identify(std::ostream& os) const
  {}

bool PHG4CylinderGeom_Siladders::load_geometry()
{
  return true;
}

void PHG4CylinderGeom_Siladders::find_segment_center(const int segment_z_bin, const int segment_phi_bin, double location[])
{
  const int itype = (segment_z_bin==1 || segment_z_bin==2) ? 0 : 1;

  // Ladder
  const double phi  = offsetphi + dphi_ * (double)segment_phi_bin;
  location[0] = eff_radius * TMath::Cos(phi);
  location[1] = eff_radius * TMath::Sin(phi);
  location[2] = ladder_z_[itype];
}

void PHG4CylinderGeom_Siladders::find_strip_center(const int segment_z_bin, const int segment_phi_bin, const int strip_column, const int strip_index, double location[])
{
  // Ladder
  find_segment_center(segment_z_bin, segment_phi_bin, location);
  TVector3 ladder(location[0], location[1], location[2]);

  // Strip
  const int itype = (segment_z_bin==1 || segment_z_bin==2) ? 0 : 1;
  const double strip_z       = strip_z_[itype];
  const int nstrips_z_sensor = nstrips_z_sensor_[itype];

  const double strip_localpos_z = 2.*strip_z*(double)(strip_column%nstrips_z_sensor) -    strip_z*(double)nstrips_z_sensor + strip_z;
  const double strip_localpos_y = 2.*strip_y*(double)strip_index                     - 2.*strip_y*(double)nstrips_phi_cell + strip_y;

  TVector3 strip_localpos(strip_x_offset, strip_localpos_y, strip_localpos_z);

  // Strip rotation
  const double phi    = offsetphi + dphi_ * (double)segment_phi_bin;
  const double rotate = phi + offsetrot + TMath::Pi();

  TRotation rot;
  rot.RotateZ(rotate);
  strip_localpos = rot*strip_localpos + ladder;

  location[0] = strip_localpos.x();
  location[1] = strip_localpos.y();
  location[2] = strip_localpos.z();
}
