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
  sensor_radius_inner(NAN),
  sensor_radius_outer(NAN),
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
  location[0] = radius  * cos(phi);
  location[1] = radius  * sin(phi);
  location[2] = signz * ladder_z_[itype];
}

void PHG4CylinderGeom_Siladders::find_strip_center(const int segment_z_bin, const int segment_phi_bin, const int strip_column, const int strip_index, double location[])
{
  // The sub layer radius is determined from the ladder phi index, called segment_phi_bin
  
  radius = sensor_radius_inner;
  if(segment_phi_bin % 2)
    radius = sensor_radius_outer;
  //cout << "      setting working sensor radius to " << radius << endl;
  
  // Ladder
  find_segment_center(segment_z_bin, segment_phi_bin, location);
  CLHEP::Hep3Vector ladder(location[0], location[1], location[2]);

  // Strip
  const int itype = segment_z_bin % 2;
  const double strip_z       = strip_z_[itype];
  const int nstrips_z_sensor = nstrips_z_sensor_[itype];

  const double strip_localpos_z = strip_z*(strip_column%nstrips_z_sensor) -    strip_z/2.*nstrips_z_sensor + strip_z/2.;
  // distance from bottom of sensor = strip_y*strip_index +strip_y/2.0, then subtract nstrips_phi_cell * strip_y / 2.0
  const double strip_localpos_y = strip_y*strip_index + strip_y/2. - nstrips_phi_cell * strip_y / 2.0;

  CLHEP::Hep3Vector strip_localpos(strip_x_offset, strip_localpos_y, strip_localpos_z);

  // Strip rotation
  const double phi    = offsetphi + dphi_ * segment_phi_bin;
  //const double rotate = phi + offsetrot + M_PI;
  const double rotate = phi + offsetrot;

  CLHEP::HepRotation rot;
  rot.rotateZ(rotate);
  strip_localpos = rot*strip_localpos;
  strip_localpos += ladder;

  location[0] = strip_localpos.x();
  location[1] = strip_localpos.y();
  location[2] = strip_localpos.z();
}

void PHG4CylinderGeom_Siladders::find_strip_index_values(const int segment_z_bin, const double yin, const double zin,  int &strip_y_index, int &strip_z_index)
{
  // Given the location in y and z in sensor local coordinates, find the strip y and z index values

  // find the sensor type (inner or outer) from the segment_z_bin (location of sensor on ladder)
  const int itype = segment_z_bin % 2;  
  if(itype != 0 && itype != 1)
    {
      cout  << "Problem: itype = " << itype << endl; 
      return;
    }

  // expect cm
  double zpos = zin;
  double ypos = yin;
  
  const double strip_z  = strip_z_[itype];
  const int nstrips_z_sensor = nstrips_z_sensor_[itype];
  const int nstrips_y_sensor = nstrips_phi_cell;
  
  // get the strip z index
  double zup = (double) nstrips_z_sensor * strip_z / 2.0 + zpos;  
  strip_z_index = (int) (zup / strip_z);

  // get the strip y index
  double yup = (double) nstrips_y_sensor * strip_y / 2.0 + ypos;
  strip_y_index = (int) (yup / strip_y);

  /*
  cout << "segment_z_bin " << segment_z_bin << " ypos " << ypos << " zpos " << zpos << " zup " << zup << " yup " << yup << endl;
  cout << "      -- itype " << itype << " strip_y " << strip_y << " strip_z " << strip_z << " nstrips_z_sensor " << nstrips_z_sensor 
       << " nstrips_y_sensor " << nstrips_y_sensor << endl;
  cout << "      --  strip_z_index " << strip_z_index << " strip_y_index " << strip_y_index << endl;
  */  
}

void PHG4CylinderGeom_Siladders::find_strip_center_localcoords(const int segment_z_bin, const int strip_y_index, const int strip_z_index, double location[])
{
  // find the sensor type (inner or outer) from the segment_z_bin (location of sensor on ladder)
  const int itype = segment_z_bin % 2;  
  if(itype != 0 && itype != 1)
    {
      cout  << "Problem: itype = " << itype << endl; 
      return;
    }

  const double strip_z  = strip_z_[itype];
  const int nstrips_z_sensor = nstrips_z_sensor_[itype];
  const int nstrips_y_sensor = nstrips_phi_cell;

  // center of strip in y
  double ypos = (double) strip_y_index * strip_y + strip_y/2.0 - (double) nstrips_y_sensor * strip_y/2.0;

  // center of strip in z
  double zpos = (double) strip_z_index * strip_z + strip_z/2.0 - (double) nstrips_z_sensor * strip_z/2.0;

  location[0] = 0.0;
  location[1] = ypos;
  location[2] = zpos;


}

