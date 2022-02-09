#include "PHG4CylinderGeomv4.h"

#include <cmath>

using namespace std;

PHG4CylinderGeomv4::PHG4CylinderGeomv4()
  : N_sensors_in_layer(-1)
  , layer(-1)
  , layer_radius(NAN)
  , radius_stagger(NAN)
  , layer_NZ(-1)
  , segment_z_step(NAN)
  , segment_phi_step(NAN)
  , sensor_x_offset(NAN)
  , sensor_y_offset(NAN)
  , N_strip_columns(-1)
  , N_strips_per_column(-1)
  , N_staggers(-1)
  , strip_z_spacing(NAN)
  , strip_y_spacing(NAN)
  , thickness(NAN)
  , strip_tilt(NAN)

{
  return;
}

void PHG4CylinderGeomv4::identify(std::ostream& os) const
{
  os << "PHG4CylinderGeomv4: layer: " << layer
     << ", layer_radius: " << layer_radius
     << ", radius_stagger: " << radius_stagger
     << ", N_sensors_in_layer: " << N_sensors_in_layer
     << ", layer_NZ: " << layer_NZ
     << ", segment_z_step: " << segment_z_step
     << ", segment_phi_step: " << segment_phi_step
     << ", sensor_x_offset: " << sensor_x_offset
     << ", sensor_y_offset: " << sensor_y_offset
     << ", N_strip_columns: " << N_strip_columns
     << ", N_strips_per_column: " << N_strips_per_column
     << ", N_staggers " << N_staggers
     << ", strip_z_spacing: " << strip_z_spacing
     << ", strip_y_spacing: " << strip_y_spacing
     << ", strip_tilt: " << strip_tilt
     << endl;
  return;
}

void PHG4CylinderGeomv4::find_segment_center(int segment_z_bin, int segment_phi_bin, double location[])
{
  double z_location = (double) (segment_z_bin - layer_NZ / 2) * segment_z_step;

  // this determines the stggered layer radius
  int istagger = segment_phi_bin % N_staggers;

  // We need to stagger the radii at alternate phi values by radius_stagger, since the ladders overlap in phi
  // The number of staggers is an input number, since it has to be the same for both parts of a double layer!
  double R_layer = layer_radius + (double) istagger * radius_stagger;

  // Place the ladder segment envelopes at the correct z and phi
  double phi = (double) segment_phi_bin * segment_phi_step;

  double x_location = R_layer * cos(phi);
  double y_location = R_layer * sin(phi);

  location[0] = x_location;
  location[1] = y_location;
  location[2] = z_location;
}

void PHG4CylinderGeomv4::find_strip_center(int segment_z_bin, int segment_phi_bin, int strip_column, int strip_index, double location[])
{
  // Start  by getting the ladder segment center location in the sPHENIX frame
  find_segment_center(segment_z_bin, segment_phi_bin, location);

  // Now calculate the strip x, y and z position in the frame of the ladder segment
  // if N_strip_columns is even, the center of the sensor is a boundary between strip columns, the first sensor is 1/2 strip_z_spacing from zero
  // if it is odd, the center of the sensor is in the middle of a strip column, one strip is centered at zero

  double strip_sensor_z = 0.0;
  if (N_strip_columns % 2)
    strip_sensor_z = ((double) (strip_column - N_strip_columns / 2)) * strip_z_spacing;
  else
    strip_sensor_z = ((double) (strip_column - N_strip_columns / 2) + 0.5) * strip_z_spacing;

  double strip_sensor_y = 0.0;
  if (N_strips_per_column % 2)
    strip_sensor_y = (double) (strip_index - N_strips_per_column / 2) * strip_y_spacing;
  else
    strip_sensor_y = ((double) (strip_index - N_strips_per_column / 2) + 0.5) * strip_y_spacing;

  // The sensor is set forward of the center in the ladder segment volume
  double strip_sensor_x = sensor_x_offset;

  // If there is only an upper ROC, the sensor is not centered in the ladder segment
  // Add the sensor offset to the strip_sensor_y value

  strip_sensor_y += sensor_y_offset;

  // Now we need to transform the position in the ladder segment frame to that in the sPHENIX frame
  // this is just a rotation around the z axis in phi
  double phi = (double) segment_phi_bin * segment_phi_step;
  double x = strip_sensor_x * cos(phi) - strip_sensor_y * sin(phi);
  double y = strip_sensor_y * cos(phi) + strip_sensor_x * sin(phi);

  // now add these to the location of the sensor center in the sPHENIX frame
  location[0] += x;
  location[1] += y;
  location[2] += strip_sensor_z;

  return;
}
