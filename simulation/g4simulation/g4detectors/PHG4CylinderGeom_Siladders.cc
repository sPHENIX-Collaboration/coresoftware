#include "PHG4CylinderGeom_Siladders.h"

#include <g4main/PHG4Detector.h>

#include <Geant4/G4RotationMatrix.hh>
#include <Geant4/G4Transform3D.hh>

#include <CLHEP/Vector/ThreeVector.h>
#include <CLHEP/Vector/Rotation.h>

#include <cmath>

using namespace std;

PHG4CylinderGeom_Siladders::PHG4CylinderGeom_Siladders():
  layer(-1),
  strip_x(NAN),
  strip_y(NAN),
  strip_z0(NAN),
  strip_z1(NAN),
  nstrips_z_sensor0(-1),
  nstrips_z_sensor1(-1),
  nstrips_phi_cell(-1),
  nladders_layer(-1),
  ladder_z0(NAN),
  ladder_z1(NAN),
  eff_radius(NAN),
  eff_radius_alternate(NAN),
  strip_x_offset(NAN),
  offsetphi(NAN),
  offsetrot(NAN),
  dphi_(NAN)
{
  strip_z_[0]  = NAN;
  strip_z_[1]  = NAN;
  ladder_z_[0] = NAN;
  ladder_z_[1] = NAN;
  nstrips_z_sensor_[0] = -1;
  nstrips_z_sensor_[1] = -1;
  return;
}

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
  const double signz = (segment_z_bin > 1) ? 1. : -1.;
  const int itype    = segment_z_bin % 2;

  // Ladder
  const double phi  = offsetphi + dphi_ * segment_phi_bin;
  location[0] = eff_radius  * cos(phi);
  location[1] = eff_radius  * sin(phi);
  location[2] = signz * ladder_z_[itype];
}

void PHG4CylinderGeom_Siladders::find_strip_center(const int segment_z_bin, const int segment_phi_bin, const int strip_column, const int strip_index, double location[])
{
  cout << endl << "Strip geom:  layer " << layer << " nstrips_z_sensor0 " << nstrips_z_sensor0 << " nstrips_z_sensor1 " 
       << nstrips_z_sensor1 << " nstrips_phi_cell " << nstrips_phi_cell
       << " ladderz0 " << ladder_z0 << " ladder_z1 " << ladder_z1 << " eff_radius " << eff_radius << " eff_radius_alternate " << eff_radius_alternate 
       << " strip_x_offset " << strip_x_offset << " offsetphi " << offsetphi 
       << endl;

    cout << "    Input :  ladder_z_index " << segment_z_bin << " ladder_phi_index " << segment_phi_bin 
       << " strip_z_index " << strip_column << " strip_y_index " << strip_index << endl;

  // Ladder
  find_segment_center(segment_z_bin, segment_phi_bin, location);
  CLHEP::Hep3Vector ladder(location[0], location[1], location[2]);

  cout << "   segment center is at " << location[0] << "  " << location[1] << "  " << location[2] << endl;

  // Strip
  const int itype = segment_z_bin % 2;
  const double strip_z       = strip_z_[itype];
  const int nstrips_z_sensor = nstrips_z_sensor_[itype];

  const double strip_localpos_z = strip_z*(strip_column%nstrips_z_sensor) -    strip_z/2.*nstrips_z_sensor + strip_z/2.;
  const double strip_localpos_y = strip_y*strip_index                     - strip_y*nstrips_phi_cell + strip_y/2.;
  cout << "   Sensor parameters: itype " << itype << " nstrips_z_sensor " << nstrips_z_sensor << " strip_z " << strip_z   
       << " nstrips_phi " << nstrips_phi_cell << " strip_y " << strip_y << endl;

  CLHEP::Hep3Vector strip_localpos(strip_x_offset, strip_localpos_y, strip_localpos_z);

  // Strip rotation
  const double phi    = offsetphi + dphi_ * segment_phi_bin;
  //const double rotate = phi + offsetrot + M_PI;
  const double rotate = phi + offsetrot;
  cout << "       phi " << phi << " offsetphi " << offsetphi << " dphi_ " << dphi_ << " segment_phi_bin " << segment_phi_bin << " rotate " << rotate << endl;

  cout << "       strip_localpos " << strip_localpos[0] << "  " << strip_localpos[1] << "  " << strip_localpos[2]  << endl;

  CLHEP::HepRotation rot;
  rot.rotateZ(rotate);
  strip_localpos = rot*strip_localpos;
  cout << "       rot*strip_localpos " << strip_localpos[0] << "  " << strip_localpos[1] << "  " << strip_localpos[2]  << endl;
  strip_localpos += ladder;

  location[0] = strip_localpos.x();
  location[1] = strip_localpos.y();
  location[2] = strip_localpos.z();

  cout << "    x " << location[0] << " y " << location[1] << " z " << location[2] << endl;
}
