#include "PHG4CylinderGeom_MVTX.h"


#include <Math/GenVector/Translation3D.h>
#include <Math/Rotation3D.h>
#include <Math/Translation3D.h>
#include <Math/Vector3D.h>
#include <TRotation.h>
#include <TVector3.h>

#include <cmath>

using namespace ROOT::Math;
using namespace std;

PHG4CylinderGeom_MVTX::PHG4CylinderGeom_MVTX(int in_layer, int in_stave_type, int in_N_staves, double in_layer_nominal_radius, double in_phistep, double in_phitilt, double in_pixel_x, double in_pixel_z, double in_pixel_thickness)
  : layer(in_layer)
  , stave_type(in_stave_type)
  , N_staves(in_N_staves)
  , layer_radius(in_layer_nominal_radius)
  , stave_phi_step(in_phistep)
  , stave_phi_tilt(in_phitilt)
  , pixel_x(in_pixel_x)
  , pixel_z(in_pixel_z)
  , pixel_thickness(in_pixel_thickness)
{
  if (stave_type < 3)
    N_half_staves = 2;
  else
    N_half_staves = 1;

  // There are two sensor sizes, one for the inner layers and one for the middle/outer layers
  // Here the sensor is at 0,0,0 with the normal to the face pointing in +y direction
  // These are half-dimensions, double them to get the full dimensions
  //    For inner layer (stave type 0):                            0.7525 x 0.0009 x 1.5050
  //    For mid and outer layer (stave types 1 and 2):   0.7500 x 0.0009 x 1.5000
  if (stave_type == 0)
  {
    Zsensor = 2.994;  //3.01;   // cm
    Xsensor = 1.376;  //1.505;   // cm
  }
  else
  {
    Zsensor = 3.0;  // cm
    Xsensor = 1.5;  // cm
  }

  /*      
  // In the ITS we should have these numbers for how staves are built (from ITS.gdml file)
   lyr rad   L   staves     modules                                                                                                     chips/module
   0  23    290    12   1 (x=0, y=0, z=0)                                                                                     9 (x=0, y=-0.00875, z=- 12.04, -9,03, -6.02, -3.01, 0, 3.01, 6.02, 9.03, 12.04) 
   1  31    290    16   1 (x=0, y=0, z=0)                                                                                     9 (x=0, y=-0.00875, z=- 12.04, -9,03, -6.02, -3.01, 0, 3.01, 6.02, 9.03, 12.04) 
   2  39    290    20   1 (x=0, y=0, z=0)                                                                                     9 (x=0, y=-0.00875, z=- 12.04, -9,03, -6.02, -3.01, 0, 3.01, 6.02, 9.03, 12.04) 
   3  194  900    24   4 (x=0, y=-0.06075, z= -31.605, -10.535, 10.535, 31.605)                    14  (x = -0.755 or +0.755, y= -0.00825, z = -9.03, -9.03, -6.02, -6.02, -3.01, -3.01, 0, 0, 3.01, 3.01, 6.02, 6.02, 9.03, 9.03)      
   4  247  900    30   4 (x=0, y=-0.06075, z= -31.605, -10.535, 10.535, 31.605)                    14  (x = -0.755 or +0.755, y= -0.00825, z = -9.03, -9.03, -6.02, -6.02, -3.01, -3.01, 0, 0, 3.01, 3.01, 6.02, 6.02, 9.03, 9.03)      
   5  253 1500   42   7 (x=0, y=-0.06075, z = -63.21, -42.14, -21.07, 0.0, 21.07, 42.14, 63.21)   14  (x = -0.755 or +0.755, y= -0.00825, z = -9.03, -9.03, -6.02, -6.02, -3.01, -3.01, 0, 0, 3.01, 3.01, 6.02, 6.02, 9.03, 9.03) 
   6  405 1500   48  7  (x=0, y=-0.06075, z = -63.21, -42.14, -21.07, 0.0, 21.07, 42.14, 63.21)   14  (x = -0.755 or +0.755, y= -0.00825, z = -9.03, -9.03, -6.02, -6.02, -3.01, -3.01, 0, 0, 3.01, 3.01, 6.02, 6.02, 9.03, 9.03) 
   sensor is in chip at (x=0, y=-0.0016, z=0)
   3-6 half-staves are in 3-6 staves at:  (x = -1.29, y = +2.067, z = 0)  or (x = +1.29 cm, y = 2.243, z = 0)
   // layers 0,1,2 have one stave with 1 module and 7 chips in that module
   // where layers 3 and 4  have two half-staves with 4 modules and 7 chips/module
   // layers 5 and 6 have two half staves with 7 modules and 7 chips/module 
   */

  // Note that stave is centered at origin with normal to face of sensor pointing in +y direction
  // Units here are cm, same as in the gdml file

  // for all layers
  //double loc_sensor_in_chip_data[3] = {0.0, -0.0016, 0.0};
  double loc_sensor_in_chip_data[3] = {0.0620, -0.0016, 0.0};  // mvtx_stave_v01.gdml

  for (int i = 0; i < 3; i++)
    loc_sensor_in_chip[i] = loc_sensor_in_chip_data[i];

  // inner barrel layers stave construction
  //==========================
  // (stave_type == 0)
  /*
  double inner_loc_chip_in_module_data[9][3] = {
    0.0, -0.00875, -12.04,
    0.0, -0.00875, -9.03,
    0.0, -0.00875, -6.02,
    0.0, -0.00875, -3.01,	   
    0.0, -0.00875, 0.0,
    0.0, -0.00875, 3.01,	   
    0.0, -0.00875, 6.02,	   
    0.0, -0.00875, 9.03,	   
    0.0, -0.00875, 12.04};	   
  double inner_loc_module_in_halfstave_data[3] = {0.0, 0.0, 0.0};   // only one module
  double inner_loc_halfstave_in_stave_data[3] = {0.0, 0.00625, 0.0}; 
  */

  // from mvtx_stave_v01.gdml
  double inner_loc_chip_in_module_data[9][3] = {
      0.0, -0.00875, -12.060,
      0.0, -0.00875, -9.0450,
      0.0, -0.00875, -6.0300,
      0.0, -0.00875, -3.0150,
      0.0, -0.00875, 0.0,
      0.0, -0.00875, 3.0150,
      0.0, -0.00875, 6.0300,
      0.0, -0.00875, 9.0450,
      0.0, -0.00875, 12.060};
  double inner_loc_module_in_halfstave_data[3] = {0.0, 0.0, 0.0};  // only one module
  double inner_loc_halfstave_in_stave_data[3] = {0.0, 0.00625, 0.0};

  for (int i = 0; i < 3; i++)
  {
    inner_loc_module_in_halfstave[i] = inner_loc_module_in_halfstave_data[i];
    inner_loc_halfstave_in_stave[i] = inner_loc_halfstave_in_stave_data[i];
    for (int j = 0; j < 9; j++)
    {
      inner_loc_chip_in_module[j][i] = inner_loc_chip_in_module_data[j][i];
    }
  }

  // middle barrel layers stave construction
  //=============================
  // (stave_type == 1)
  double middle_loc_chip_in_module_data[14][3] = {
      -0.755, -0.00825, -9.03,
      0.755, -0.00825, -9.03,
      -0.755, -0.00825, -6.02,
      0.755, -0.00825, -6.02,
      -0.755, -0.00825, -3.01,
      0.755, -0.00825, -3.01,
      -0.755, -0.00825, 0.0,
      0.755, -0.00825, 0.0,
      -0.755, -0.00825, 3.01,
      0.755, -0.00825, 3.01,
      -0.755, -0.00825, 6.02,
      0.755, -0.00825, 6.02,
      -0.755, -0.00825, 9.03,
      0.755, -0.00825, 9.03};
  double middle_loc_module_in_halfstave_data[4][3] = {
      0.0, -0.06075, -31.605,
      0.0, -0.06075, -10.535,
      0.0, -0.06075, 10.535,
      0.0, -0.06075, 31.605};
  double middle_loc_halfstave_in_stave_data[2][3] = {
      -1.29, 2.067, 0.0,
      1.29, 2.243, 0.0};

  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 14; j++)
      middle_loc_chip_in_module[j][i] = middle_loc_chip_in_module_data[j][i];

    for (int j = 0; j < 4; j++)
      middle_loc_module_in_halfstave[j][i] = middle_loc_module_in_halfstave_data[j][i];

    for (int j = 0; j < 2; j++)
      middle_loc_halfstave_in_stave[j][i] = middle_loc_halfstave_in_stave_data[j][i];
  }

  // outer barrel layers stave construction
  //===========================
  // (stave_type == 2)
  double outer_loc_chip_in_module_data[14][3] = {
      -0.755, -0.00825, -9.03,
      0.755, -0.00825, -9.03,
      -0.755, -0.00825, -6.02,
      0.755, -0.00825, -6.02,
      -0.755, -0.00825, -3.01,
      0.755, -0.00825, -3.01,
      -0.755, -0.00825, 0.0,
      0.755, -0.00825, 0.0,
      -0.755, -0.00825, 3.01,
      0.755, -0.00825, 3.01,
      -0.755, -0.00825, 6.02,
      0.755, -0.00825, 6.02,
      -0.755, -0.00825, 9.03,
      0.755, -0.00825, 9.03};
  double outer_loc_module_in_halfstave_data[7][3] = {
      0.0, -0.06075, -63.21,
      0.0, -0.06075, -42.14,
      0.0, -0.06075, -21.07,
      0.0, -0.06075, 0.0,
      0.0, -0.06075, 21.07,
      0.0, -0.06075, 42.14,
      0.0, -0.06075, 63.21};
  double outer_loc_halfstave_in_stave_data[2][3] = {
      -1.29, 2.067, 0.0,
      1.29, 2.243, 0.0};

  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 14; j++)
      outer_loc_chip_in_module[j][i] = outer_loc_chip_in_module_data[j][i];

    for (int j = 0; j < 7; j++)
      outer_loc_module_in_halfstave[j][i] = outer_loc_module_in_halfstave_data[j][i];

    for (int j = 0; j < 2; j++)
      outer_loc_halfstave_in_stave[j][i] = outer_loc_halfstave_in_stave_data[j][i];
  }

  return;
}

TVector3
PHG4CylinderGeom_MVTX::get_local_from_world_coords(int stave, int half_stave, int module, int chip, TVector3 world_location)
{
  double stave_phi = stave_phi_step * (double) stave;
  double stave_phi_offset = M_PI / 2.0;  // stave initially points so that sensor faces upward in y

  /*
    cout << endl << "PHG4CylinderGeom_MVTX::get_local_from_world_coords: " << " Stave type " << stave_type 
	 << " chip " << chip 
	 << " world coords " << world_location.X() << " " << world_location.Y() << " " << world_location.Z() << endl;
  */

  TVector3 res;

  // Is this is an inner, middle or outer stave?
  if (stave_type == 0)
  {
    // Inner stave

    // transform location of stave from its location in the world - this is just a translation
    TVector3 tr4(layer_radius * cos(stave_phi), layer_radius * sin(stave_phi), 0.0);
    res = world_location - tr4;

    // Rotate stave from its angle in the world
    // This requires rotating it by:
    //   removing the tilt (for layers 0-2)
    //   removing the angle that makes it point at the origin when it was at it's location in the world
    //   rotating it by -90 degrees to make the face point vertically up in y
    TRotation R;
    R.RotateZ(-stave_phi - stave_phi_offset - stave_phi_tilt);
    res = R * res;  // rotates res using R

    // transform location in stave to location in half stave
    TVector3 tr3(inner_loc_halfstave_in_stave[0], inner_loc_halfstave_in_stave[1], inner_loc_halfstave_in_stave[2]);
    res = res - tr3;

    // transfor from half stave to module
    TVector3 tr2a(inner_loc_module_in_halfstave[0], inner_loc_module_in_halfstave[1], inner_loc_module_in_halfstave[2]);
    res = res - tr2a;

    // transform location in module to location in chip
    TVector3 tr2(inner_loc_chip_in_module[chip][0], inner_loc_chip_in_module[chip][1], inner_loc_chip_in_module[chip][2]);
    res = res - tr2;

    // transform location in chip to location in sensor
    TVector3 tr1(loc_sensor_in_chip[0], loc_sensor_in_chip[1], loc_sensor_in_chip[2]);
    res = res - tr1;
  }
  else if (stave_type == 1)
  {
    // transform location of stave from its location in the world - this is just a translation
    TVector3 tr4(layer_radius * cos(stave_phi), layer_radius * sin(stave_phi), 0.0);
    res = world_location - tr4;

    // Rotate stave from its angle in the world
    // This requires rotating it by:
    //     removing the tilt (non-zero tilt for layers 0-2)
    //     removing the angle that makes it point at the origin when it was at it's location in the world
    //     rotating it by -90 degrees to make the face point vertically up in y
    TRotation R;
    R.RotateZ(-stave_phi - stave_phi_offset - stave_phi_tilt);
    res = R * res;  // rotates res using R

    // transform location in stave to location in half stave
    TVector3 tr3(middle_loc_halfstave_in_stave[half_stave][0], middle_loc_halfstave_in_stave[half_stave][1], middle_loc_halfstave_in_stave[half_stave][2]);
    res = res - tr3;

    // transform from half stave to module
    TVector3 tr2a(middle_loc_module_in_halfstave[module][0], middle_loc_module_in_halfstave[module][1], middle_loc_module_in_halfstave[module][2]);
    res = res - tr2a;

    // transform location in module to location in chip
    TVector3 tr2(middle_loc_chip_in_module[chip][0], middle_loc_chip_in_module[chip][1], middle_loc_chip_in_module[chip][2]);
    res = res - tr2;

    // for odd numbered chips, the chip is flipped by 180 degrees around the x and z axes
    TRotation Rchip, Rchip_inv;
    if (chip % 2)
    {
      Rchip.RotateX(M_PI);
      Rchip.RotateZ(M_PI);
      res = Rchip * res;
    }

    // transform location in chip to location in sensor
    TVector3 tr1(loc_sensor_in_chip[0], loc_sensor_in_chip[1], loc_sensor_in_chip[2]);
    res = res - tr1;
  }
  else if (stave_type == 2)
  {
    // transform location of stave from its location in the world - this is just a translation
    TVector3 tr4(layer_radius * cos(stave_phi), layer_radius * sin(stave_phi), 0.0);
    res = world_location - tr4;

    // Rotate stave from its angle in the world
    // This requires rotating it by:
    //   removing the tilt (for layers 0-2)
    //   removing the angle that makes it point at the origin when it was at it's location in the world
    //   rotating it by -90 degrees to make the face point vertically up in y
    TRotation R;
    R.RotateZ(-stave_phi - stave_phi_offset - stave_phi_tilt);
    res = R * res;  // rotates res using R

    // transform location in stave to location in half stave
    TVector3 tr3(outer_loc_halfstave_in_stave[half_stave][0], outer_loc_halfstave_in_stave[half_stave][1], outer_loc_halfstave_in_stave[half_stave][2]);
    res = res - tr3;

    // transform from half stave to module
    TVector3 tr2a(outer_loc_module_in_halfstave[module][0], outer_loc_module_in_halfstave[module][1], outer_loc_module_in_halfstave[module][2]);
    res = res - tr2a;

    // transform location in module to location in chip
    TVector3 tr2(outer_loc_chip_in_module[chip][0], outer_loc_chip_in_module[chip][1], outer_loc_chip_in_module[chip][2]);
    res = res - tr2;

    // for odd numbered chips, the chip is flipped by 180 degrees around the x and z axes
    TRotation Rchip, Rchip_inv;
    if (chip % 2)
    {
      Rchip.RotateX(M_PI);
      Rchip.RotateZ(M_PI);
      res = Rchip * res;
    }

    // transform location in chip to location in sensor
    TVector3 tr1(loc_sensor_in_chip[0], loc_sensor_in_chip[1], loc_sensor_in_chip[2]);
    res = res - tr1;
  }

  return res;
}

TVector3
PHG4CylinderGeom_MVTX::get_world_from_local_coords(int stave, int half_stave, int module, int chip, TVector3 sensor_local)
{
  double stave_phi = stave_phi_step * (double) stave;
  double stave_phi_offset = M_PI / 2.0;  // stave initially points so that sensor faces upward in y

  /*  
    cout << endl << "PHG4CylinderGeom_MVTX::get_world_from_local_coords: " << " stave type " << stave_type
	 << " chip " << chip
	 << " local coords " << sensor_local.X() << " " << sensor_local.Y() << " " << sensor_local.Z() << endl;
  */

  // Is this is an inner, middle or outer stave?
  if (stave_type == 0)
  {
    // Inner stave
    // Start with the point in sensor local coords
    TVector3 pos1 = sensor_local;

    // transform sensor location to location in chip
    TVector3 tr1(loc_sensor_in_chip[0], loc_sensor_in_chip[1], loc_sensor_in_chip[2]);
    TVector3 res = pos1 + tr1;

    // transform location in chip to location in module
    TVector3 tr2(inner_loc_chip_in_module[chip][0], inner_loc_chip_in_module[chip][1], inner_loc_chip_in_module[chip][2]);
    res = res + tr2;

    // module to half stave
    TVector3 tr2a(inner_loc_module_in_halfstave[0], inner_loc_module_in_halfstave[1], inner_loc_module_in_halfstave[2]);
    res = res + tr2a;

    // transform location in half stave to location in stave
    TVector3 tr3(inner_loc_halfstave_in_stave[0], inner_loc_halfstave_in_stave[1], inner_loc_halfstave_in_stave[2]);
    res = res + tr3;

    // Rotate stave to its angle in the world
    // This requires rotating it by
    //    90 degrees to make the face point to the origin instead of vertically up in y when it is at phi = 0 - stave_phi_offset is 90 degrees in CCW direction
    //    Rotating it fiurther so that it points at the origin after being translated to the x and y coords of the stave phi location - stave_phi derived from
    //    stave_phi_step  and stave (number), both constructor parameters
    //    Adding the tilt (for layers 0-2) - stave_phi_tilt is a constructor parameter provided by PHG4MVTXDetector
    // for a rotation
    TRotation R;
    R.RotateZ(stave_phi + stave_phi_offset + stave_phi_tilt);
    res = R * res;  // rotates res using R

    // transform location of stave to its location in the world
    TVector3 tr4(layer_radius * cos(stave_phi), layer_radius * sin(stave_phi), 0.0);
    res = res + tr4;

    return res;
  }
  else if (stave_type == 1)
  {
    // Middle stave
    // Start with the point in sensor local coords
    TVector3 pos1 = sensor_local;

    // transform sensor location to location in chip
    TVector3 tr1(loc_sensor_in_chip[0], loc_sensor_in_chip[1], loc_sensor_in_chip[2]);
    TVector3 res = pos1 + tr1;

    // transform location in chip to location in module
    // for odd numbered chips, the chip is flipped by 180 degrees around the x and z axes
    TRotation Rchip;
    if (chip % 2)
    {
      Rchip.RotateX(-M_PI);
      Rchip.RotateZ(-M_PI);
      res = Rchip * res;
    }

    TVector3 tr2(middle_loc_chip_in_module[chip][0], middle_loc_chip_in_module[chip][1], middle_loc_chip_in_module[chip][2]);
    res = res + tr2;

    // module to half stave
    TVector3 tr2a(middle_loc_module_in_halfstave[module][0], middle_loc_module_in_halfstave[module][1], middle_loc_module_in_halfstave[module][2]);
    res = res + tr2a;

    // transform location in half stave to location in stave
    TVector3 tr3(middle_loc_halfstave_in_stave[half_stave][0], middle_loc_halfstave_in_stave[half_stave][1], middle_loc_halfstave_in_stave[half_stave][2]);
    res = res + tr3;

    // Rotate stave to its angle in the world
    //    90 degrees to make it point to the origin instead of vertically up in y when it is at phi = 0
    //    Rotating it fiurther so that it points at the origin after being translated to the x and y coords of the stave phi location
    //    adding the tilt (for layers 0-2)
    TRotation R;
    R.RotateZ(stave_phi + stave_phi_offset + stave_phi_tilt);
    res = R * res;  // rotates res using R

    // transform location of stave to its location in the world
    TVector3 tr4(layer_radius * cos(stave_phi), layer_radius * sin(stave_phi), 0.0);
    res = res + tr4;

    return res;
  }
  else
  {
    // stave_type = 2
    // Outer stave
    //============

    // Start with the point in sensor local coords
    TVector3 pos1 = sensor_local;

    // transform sensor location to location in chip
    TVector3 tr1(loc_sensor_in_chip[0], loc_sensor_in_chip[1], loc_sensor_in_chip[2]);
    TVector3 res = pos1 + tr1;

    // transform location in chip to location in module
    // for odd numbered chips, the chip is flipped by 180 degrees around the x and z axes
    TRotation Rchip;
    if (chip % 2)
    {
      Rchip.RotateX(-M_PI);
      Rchip.RotateZ(-M_PI);
      res = Rchip * res;
    }

    TVector3 tr2(outer_loc_chip_in_module[chip][0], outer_loc_chip_in_module[chip][1], outer_loc_chip_in_module[chip][2]);
    res = res + tr2;

    // module to half stave
    TVector3 tr2a(outer_loc_module_in_halfstave[module][0], outer_loc_module_in_halfstave[module][1], outer_loc_module_in_halfstave[module][2]);
    res = res + tr2a;

    // transform location in half stave to location in stave
    TVector3 tr3(outer_loc_halfstave_in_stave[half_stave][0], outer_loc_halfstave_in_stave[half_stave][1], outer_loc_halfstave_in_stave[half_stave][2]);
    res = res + tr3;

    // Rotate stave to its angle in the world. Rotate by:
    //    90 degrees to make it point to the origin instead of vertically up in y when it is at phi = 0
    //    Rotating it fiurther so that it points at the origin after being translated to the x and y coords of the stave phi location
    //    adding the tilt (for layers 0-2)
    TRotation R;
    R.RotateZ(stave_phi + stave_phi_offset + stave_phi_tilt);
    res = R * res;  // rotates res using R

    // transform location of stave to its location in the world
    TVector3 tr4(layer_radius * cos(stave_phi), layer_radius * sin(stave_phi), 0.0);
    res = res + tr4;

    return res;
  }
}

int PHG4CylinderGeom_MVTX::get_pixel_number_from_xbin_zbin(int xbin, int ybin)
{
  //NZ = (int)  ( Zsensor / (pixel_z) );
  //  NX = (int)  ( Xsensor / (pixel_x) );

  int NXZ = xbin + ybin * get_NX();

  return NXZ;
}

int PHG4CylinderGeom_MVTX::get_pixel_X_from_pixel_number(int NXZ)
{
  //  NZ = (int)  ( Zsensor / (pixel_z) );
  //  NX = (int)  ( Xsensor / (pixel_x) );

  int Ngridx = NXZ % get_NX();

  return Ngridx;
}

int PHG4CylinderGeom_MVTX::get_pixel_Z_from_pixel_number(int NXZ)
{
  //  NZ = (int)  ( Zsensor / (pixel_z) );
  //  NX = (int)  ( Xsensor / (pixel_x) );

  int Ngridz = NXZ / get_NX();

  return Ngridz;
}

int PHG4CylinderGeom_MVTX::get_pixel_from_local_coords(TVector3 sensor_local)
{
  //  NZ = (int)  ( Zsensor / (pixel_z) );
  //  NX = (int)  ( Xsensor / (pixel_x) );

  //cout  << " Pixels in X: NX  " << NX  << " pixels in Z: NZ  " << NZ << endl;

  // start pixel numbering from the middle of the sensor
  // find the pixel grid point

  double npix_x = sensor_local.X() / pixel_x;
  int Ngridx = int(npix_x);

  double npix_z = sensor_local.Z() / pixel_z;
  int Ngridz = int(npix_z);

  //  Combine the grid locations into a single integer
  // transform to the grid location referenced to top left corner of the chip as (0,0)
  Ngridx += get_NX() / 2;
  Ngridz += get_NZ() / 2;

  /*
  cout << "Transformed grid locations: " 
       << " Ngridx (ref to neg x, neg y corner) " << Ngridx
       << " Ngridz (ref to neg x, neg y corner) " << Ngridz
       << endl;
  */

  // numbering starts at zero
  int NXZ = Ngridx + Ngridz * get_NX();
  //cout << " pixel number is " << NXZ << endl;

  return NXZ;
}

TVector3 PHG4CylinderGeom_MVTX::get_local_coords_from_pixel(int NXZ)
{
  //  NZ = (int)  ( Zsensor / (pixel_z) );
  //  NX = (int)  ( Xsensor / (pixel_x) );

  //cout  << " Pixels in X: NX  " << NX  << " pixels in Z: NZ  " << NZ << endl;

  int Ngridz = NXZ / get_NX();
  int Ngridx = NXZ % get_NX();

  // change to a grid centered on the sensor
  Ngridx -= get_NX() / 2;
  Ngridz -= get_NZ() / 2;

  double sensor_local_x = (double) Ngridx * pixel_x;
  if (sensor_local_x < 0)
    sensor_local_x -= pixel_x / 2.0;
  else
    sensor_local_x += pixel_x / 2.0;

  double sensor_local_z = (double) Ngridz * pixel_z;
  if (sensor_local_z < 0)
    sensor_local_z -= pixel_z / 2.0;
  else
    sensor_local_z += pixel_z / 2.0;

  // The front of the sensor is at y = 0.0009, the back is at y = -0.0009
  // if we wanted the coordinates of the entrance or exit point, we would make y = to 0.0009 or -0.0009
  // Because we want the coords of the center of the pixel, we take y = 0
  TVector3 sensor_local_coords(sensor_local_x, 0.0, sensor_local_z);

  return sensor_local_coords;
}

void PHG4CylinderGeom_MVTX::identify(std::ostream& os) const
{
  os << "PHG4CylinderGeom_MVTX: layer: " << layer
     << ", layer_radius: " << layer_radius
     << " , stave_type " << stave_type
     << ", N_staves in layer: " << N_staves
     << ", N_half_staves in layer: " << N_half_staves
     << ", pixel_x: " << pixel_x
     << ", pixel_z: " << pixel_z
     << ", pixel_thickness: " << pixel_thickness
     << endl;
  return;
}

int PHG4CylinderGeom_MVTX::get_ladder_z_index(int module, int chip)
{
  // return the z index of the sensor on the stave
  // this depends on the stave type
  // stave_type 0 has 1 module with 9 chips / module in 9 z locations (no half staves)
  // stave_type 1 has 2 half-staves separted in azimuth, with 4 modules each lined up in z, with 14 chips / module paired in 7 z locations,  the pairs separated in azimuth
  // stave_type 2 has 2 half staves separated in azimuth, each with 7 modules lined up in z, with 14 chips / module paired in 7 z locations, the pairs separated in azimuth

  int ladder_z_index;

  if (stave_type == 0)
  {
    // only 1 module
    ladder_z_index = chip;
  }
  else
  {
    // stave_type 1 or 2
    // half stave number does not affect z location, it affects only phi location
    ladder_z_index = module * 7 + chip;
  }

  return ladder_z_index;
}

int PHG4CylinderGeom_MVTX::get_ladder_phi_index(int stave, int half_stave, int chip)
{
  // return the phi index of the sensor on the stave
  // this depends on the stave type
  // stave_type 0 has 1 module with 9 chips / module in 9 z locations (no half staves)
  // stave_type 1 has 2 half-staves separted in azimuth, with 4 modules each lined up in z, with 14 chips / module paired in 7 z locations,  the pairs separated in azimuth
  // stave_type 2 has 2 half staves separated in azimuth, each with 7 modules lined up in z, with 14 chips / module paired in 7 z locations, the pairs separated in azimuth

  int ladder_phi_index;

  if (stave_type == 0)
  {
    // no half staves
    ladder_phi_index = stave;
  }
  else
  {
    // stave_type 1 or 2
    // each stave has two half staves separated in azimuth, with two rows of chips each also separated in azimuth
    ladder_phi_index = stave * 4 + half_stave * 2 + chip % 2;
  }

  return ladder_phi_index;
}

void PHG4CylinderGeom_MVTX::find_sensor_center(int stave, int half_stave, int module, int chip, double location[])
{
  TVector3 sensor_local(0.0, 0.0, 0.0);

  TVector3 sensor_world = get_world_from_local_coords(stave, half_stave, module, chip, sensor_local);

  location[0] = sensor_world.X();
  location[1] = sensor_world.Y();
  location[2] = sensor_world.Z();

  return;
}
