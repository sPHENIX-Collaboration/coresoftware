#include "CylinderGeom_Mvtx.h"

#include "SegmentationAlpide.h"

#include <TRotation.h>
#include <TVector3.h>

#include <cmath>
#include <ostream>  // for operator<<, basic_ostream::operator<<, basic_...

using namespace std;
using Segmentation = SegmentationAlpide;

CylinderGeom_Mvtx::CylinderGeom_Mvtx(
    int in_layer,
    int in_N_staves,
    double in_layer_nominal_radius,
    double in_phistep,
    double in_phitilt,
    double in_phi0)
  : layer(in_layer)
  , N_staves(in_N_staves)
  , N_half_staves(0)
  , layer_radius(in_layer_nominal_radius)
  , stave_phi_step(in_phistep)
  , stave_phi_tilt(in_phitilt)
  , stave_phi_0(in_phi0)
  , pixel_x(Segmentation::PitchRow)
  , pixel_z(Segmentation::PitchCol)
  , pixel_thickness(Segmentation::SensorLayerThickness)
{
  // Note that stave is centered at origin with normal to face of sensor pointing in +y direction
  // Units here are cm, same as in the gdml file

  // for all layers
  double loc_sensor_in_chip_data[3] = {0.058128, -0.0005, 0.0};  // mvtx_stave_v1.gdml

  for (int i = 0; i < 3; i++)
  {
    loc_sensor_in_chip[i] = loc_sensor_in_chip_data[i];
  }

  // inner barrel layers stave construction
  //==========================
  // from mvtx_stave_v1.gdml
  // chip 0 is the closet to connectors (-Z)
  double inner_loc_chip_in_module_data[9][3] = {
      {0.0275, -0.02075, -12.060},
      {0.0275, -0.02075, -9.0450},
      {0.0275, -0.02075, -6.0300},
      {0.0275, -0.02075, -3.0150},
      {0.0275, -0.02075, 0.0},
      {0.0275, -0.02075, 3.0150},
      {0.0275, -0.02075, 6.0300},
      {0.0275, -0.02075, 9.0450},
      {0.0275, -0.02075, 12.060}};

  double inner_loc_module_in_halfstave_data[3] = {0.0, 0.0, 0.0};  // only one module

  double inner_loc_halfstave_in_stave_data[3] = {-0.0275, 0.01825, 0.0};

  for (int i = 0; i < 3; i++)
  {
    inner_loc_module_in_halfstave[i] = inner_loc_module_in_halfstave_data[i];
    inner_loc_halfstave_in_stave[i] = inner_loc_halfstave_in_stave_data[i];
    for (int j = 0; j < 9; j++)
    {
      inner_loc_chip_in_module[j][i] = inner_loc_chip_in_module_data[j][i];
    }
  }

  return;
}

void CylinderGeom_Mvtx::get_sensor_indices_from_world_coords(std::vector<double>& world, unsigned int& stave_index, unsigned int& chip_index)
{
  // stave number is fom phi
  double phi = atan2(world[1], world[0]);
  if (phi < 0)
  {
    phi += 2.0 * M_PI;
  }

  // int stave_tmp = (int) ( (phi - stave_phi_0) / stave_phi_step );
  int stave_tmp = round((phi - stave_phi_0) / stave_phi_step);
  // std::cout << "  phi " << phi << " stave_phi_0 " << stave_phi_0 << " stave_phi_step " << stave_phi_step << " stave_tmp " << stave_tmp << std::endl;

  // sensor is from z
  double chip_delta_z = (inner_loc_chip_in_module[8][2] - inner_loc_chip_in_module[0][2]) / 8.0;
  // int chip_tmp = (int) (world[2]/chip_delta_z) + 4;  // 0-9
  int chip_tmp = round(world[2] / chip_delta_z) + 4;  // 0-9
  // std::cout << "  z " << world[2] << " chip_delta_z " << chip_delta_z << " chip_tmp " << chip_tmp << endl;

  stave_index = stave_tmp;
  chip_index = chip_tmp;
}

bool CylinderGeom_Mvtx::get_pixel_from_local_coords(TVector3 sensor_local, int& iRow, int& iCol)
{
  // YCM (2020-01-02): It seems that due some round issues, local coords of hits at the edge of the sensor volume
  //                   are out by some fraction of microns from the ActiveMatrix. Making a safety check inside 0.1 um
  double EPS = 5e-6;
  if (fabs(fabs(sensor_local.X()) - SegmentationAlpide::ActiveMatrixSizeRows / 2.F) < EPS)
  {
    // cout << " Adjusting X,  before X= " << sensor_local.X() << endl;
    sensor_local.SetX(((sensor_local.X() < 0) ? -1 : 1) * (SegmentationAlpide::ActiveMatrixSizeRows / 2.F - EPS));
    // cout << " Adjusting X,  after X= " << sensor_local.X() << endl;
  }
  if (fabs(fabs(sensor_local.Z()) - SegmentationAlpide::ActiveMatrixSizeCols / 2.F) < EPS)
  {
    // cout << " Adjusting Z,  before Z= " << sensor_local.Z() << endl;
    sensor_local.SetZ(((sensor_local.Z() < 0) ? -1 : 1) * (SegmentationAlpide::ActiveMatrixSizeCols / 2.F - EPS));
    // cout << " Adjusting Z,  after Z= " << sensor_local.Z() << endl;
  }
  // YCM (2020-01-02): go from sensor to chip local coords
  TVector3 in_chip = sensor_local;
  TVector3 tr(loc_sensor_in_chip[0], loc_sensor_in_chip[1], loc_sensor_in_chip[2]);
  in_chip += tr;

  return SegmentationAlpide::localToDetector(in_chip.X(), in_chip.Z(), iRow, iCol);
}

int CylinderGeom_Mvtx::get_pixel_from_local_coords(const TVector3& sensor_local)
{
  int Ngridx, Ngridz;
  bool px_in = get_pixel_from_local_coords(sensor_local, Ngridx, Ngridz);
  if (!px_in)
  {
    cout << PHWHERE
         << " Pixel is out sensor. ("
         << sensor_local.X() << ", "
         << sensor_local.Y() << ", "
         << sensor_local.Z() << ")."
         << endl;
  }

  if (Ngridx < 0 || Ngridx >= get_NX() || Ngridz < 0 || Ngridz >= get_NZ())
  {
    cout << PHWHERE << "Wrong pixel value X= " << Ngridx << " and Z= " << Ngridz << endl;
  }

  // numbering starts at zero
  return Ngridx + Ngridz * get_NX();
}

TVector3 CylinderGeom_Mvtx::get_local_coords_from_pixel(int NXZ)
{
  int Ngridz = NXZ / get_NX();
  int Ngridx = NXZ % get_NX();

  return get_local_coords_from_pixel(Ngridx, Ngridz);
}

TVector3 CylinderGeom_Mvtx::get_local_coords_from_pixel(int iRow, int iCol)
{
  TVector3 local;
  bool check = SegmentationAlpide::detectorToLocal((float) iRow, (float) iCol, local);
  if (!check)
  {
    cout << PHWHERE << "Pixel coord ( " << iRow << ", " << iCol << " )"
         << "out of range" << endl;
  }
  // Transform location in chip to location in sensors
  TVector3 trChipToSens(loc_sensor_in_chip[0],
                        loc_sensor_in_chip[1],
                        loc_sensor_in_chip[2]);
  local -= trChipToSens;
  return local;
}

void CylinderGeom_Mvtx::identify(std::ostream& os) const
{
  os << "CylinderGeom_Mvtx: layer: " << layer
     << ", layer_radius: " << layer_radius
     << ", N_staves in layer: " << N_staves
     << ", N_half_staves in layer: " << N_half_staves
     << ", pixel_x: " << pixel_x
     << ", pixel_z: " << pixel_z
     << ", pixel_thickness: " << pixel_thickness
     << endl;
  return;
}

// hide SegmentationAlpide include from root5 rootcint (rootcling is fine)
int CylinderGeom_Mvtx::get_NZ() const
{
  return SegmentationAlpide::NCols;
}

int CylinderGeom_Mvtx::get_NX() const
{
  return SegmentationAlpide::NRows;
}

int CylinderGeom_Mvtx::get_pixel_X_from_pixel_number(int NXZ)
{
  return NXZ % get_NX();
}

int CylinderGeom_Mvtx::get_pixel_Z_from_pixel_number(int NXZ)
{
  return NXZ / get_NX();
}

int CylinderGeom_Mvtx::get_pixel_number_from_xbin_zbin(int xbin, int zbin)  // obsolete
{
  return xbin + zbin * get_NX();
}
