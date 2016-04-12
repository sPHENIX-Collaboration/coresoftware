#include "PHG4CylinderGeom_MAPS.h"
#include "Math/Vector3D.h"
#include "TVector3.h"
#include "TRotation.h"
#include "Math/Translation3D.h"
#include "Math/GenVector/Translation3D.h"
#include "Math/Rotation3D.h"
#include <cmath>

ClassImp(PHG4CylinderGeom_MAPS)

using namespace ROOT::Math;
using namespace std;

PHG4CylinderGeom_MAPS::PHG4CylinderGeom_MAPS(int in_layer, int in_stave_type, int in_N_staves, double in_layer_nominal_radius, double in_phistep, double in_phitilt):
  layer(in_layer),
  stave_type(in_stave_type),
  N_staves(in_N_staves),
  layer_radius(in_layer_nominal_radius),
  stave_phi_step(in_phistep),
  stave_phi_tilt(in_phitilt)
{
  // construct the geometry
  //construct_geometry();

  return;
}

TVector3 
PHG4CylinderGeom_MAPS::get_world_from_local_coords(int stave, int half_stave, int module, int chip, TVector3 sensor_local)
{
  // ITS inner barrel layer geometry

 /*      
  // In the ITS we should have these numbers for how staves are built (from ITS.gdml file)
   lyr rad   L   staves     modules                                                                                                     chips/module
   0  23    290    12    1 (x=0, y=0, z=0)                                                                                     9 (x=0, y=-0.00875, z=- 12.04, -9,03, -6.02, -3.01, 0, 3.01, 6.02, 9.03, 12.04) 
   1  31    290    16   1 (x=0, y=0, z=0)                                                                                     9 (x=0, y=-0.00875, z=- 12.04, -9,03, -6.02, -3.01, 0, 3.01, 6.02, 9.03, 12.04) 
   2  39    290    20   1 (x=0, y=0, z=0)                                                                                     9 (x=0, y=-0.00875, z=- 12.04, -9,03, -6.02, -3.01, 0, 3.01, 6.02, 9.03, 12.04) 
   3  194  900    24   4 (x=0, y=-0.06075, z= -31.605, -10.535, 10.535, 31.605)                    14  (x = -0.755 or +0.755, y= -0.00825, z = -9.03, -9.03, -6.02, -6.02, -3.01, -3.01, 0, 0, 3.01, 3.01, 6.02, 6.02, 9.03, 9.03)      
   4  247  900    30   4 (x=0, y=-0.06075, z= -31.605, -10.535, 10.535, 31.605)                    14  (x = -0.755 or +0.755, y= -0.00825, z = -9.03, -9.03, -6.02, -6.02, -3.01, -3.01, 0, 0, 3.01, 3.01, 6.02, 6.02, 9.03, 9.03)      
   5  253 1500   42   7 (x=0, y=-0.06075, z = -63.21, -42.14, -21.07, 0.0, 21.07, 42.14, 63.21)   14  (x = -0.755 or +0.755, y= -0.00825, z = -9.03, -9.03, -6.02, -6.02, -3.01, -3.01, 0, 0, 3.01, 3.01, 6.02, 6.02, 9.03, 9.03) 
   6  405 1500  48  7 (x=0, y=-0.06075, z = -63.21, -42.14, -21.07, 0.0, 21.07, 42.14, 63.21)   14  (x = -0.755 or +0.755, y= -0.00825, z = -9.03, -9.03, -6.02, -6.02, -3.01, -3.01, 0, 0, 3.01, 3.01, 6.02, 6.02, 9.03, 9.03) 
   sensor is in chip at (x=0, y=-0.0016, z=0)
   3-6 half-staves are in 3-6 staves at:  (x = -1.29, y = +2.067, z = 0)  or (x = +1.29 cm, y = 2.243, z = 0)
   // layers 0,1,2 have one stave with 1 module and 7 chips in that module
   // where layers 3 and 4  have two half-staves with 4 modules and 7 chips/module
   // layers 5 and 6 have two half staves with 7 modules and 7 chips/module 
   */

  // Note that stave is centered at origin with normal to face of sensor pointing in +y direction

  // UNITS? The required units here are cm, I think, and that is what the placement instructions in the gdml file use

  // for all layers
  double loc_sensor_in_chip[3] = {0.0, -0.0016, 0.0};

  // inner barrel layers stave construction 
  //========================== 
  // (stave_type == 0)
  double inner_loc_chip_in_module[9][3] = {
    0.0, -0.00875, -12.04,
    0.0, -0.00875, -9.03,
    0.0, -0.0875, -6.02,
    0.0, -0.0875, -3.01,	   
    0.0, -0.0875, 0.0,
    0.0, -0.0875, 3.01,	   
    0.0, -0.0875, 6.02,	   
    0.0, -0.0875, 9.03,	   
    0.0, -0.0875, 12.04};	   
  double inner_loc_module_in_halfstave[3] = {0.0, 0.0, 0.0};   // only one module
  double inner_loc_halfstave_in_stave[3] = {0.0, 0.00625, 0.0}; 

  // middle barrel layers stave construction 
  //=============================
  // (stave_type == 1)
  // indices are half_stave, chip
  double middle_loc_chip_in_module[14][3] = {
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
    0.755, -0.00825, 9.03  };
  double middle_loc_module_in_halfstave[4][3] = {
    0.0, -0.06075, -31.605,
    0.0, -0.06075, -10.535,
    0.0, -0.06075, 10.535,
    0.0, -0.06075, 31.605};
  double middle_loc_halfstave_in_stave[2][3] = {
    -1.29, 2.067, 0.0,
    1.29, 2.243, 0.0 };

  // outer barrel layers stave construction 
  //===========================
  // (stave_type == 2)
  double outer_loc_chip_in_module[14][3] = {
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
    0.755, -0.00825, 9.03 };
  double outer_loc_module_in_halfstave[7][3] = {
    0.0, -0.06075, -63.21,
    0.0, -0.06075, -42.14,
    0.0, -0.06075, -21.07,
    0.0, -0.06075, 0.0,
    0.0, -0.06075, 21.07,
    0.0, -0.06075, 42.14,
    0.0, -0.06075, 63.21 };
  double outer_loc_halfstave_in_stave[2][3] = {
    -1.29, 2.067, 0.0,
    1.29, 2.243, 0.0 };

  // This method will be called only for the layer and stave type for which it was constructed
  // Assume we are given the stave number, half stave number, module number and chip number
  // from that we can find the sensor center by the following procedure:
  //   Find the location of the sensor in the chip using "loc_sensor_in_chip"
  //   transform the coord in the chip to the coord in the module using "loc_chip_in_module"
  //   transform the coord in the module to the coord in the halfstave using "loc_module_in_halfstave"
  //   transform the coord in the halfstave to the coord in the stave using "loc_halfstave_in_stave"
  //   transform the coord in the stave to the coord in the world using  layer_radius, stave_number, stave_phi_step, stave_phi_tilt
  
  // Assume we are given:
  // chip number, module number, half stave number, stave number, local position in sensor (xsensor, ysensor, zsensor)
  // this geometry module knows the layer radius, stave phistep, stave phi offset, stave phitilt
  // we want the world coords of this hit

  double stave_phi = stave_phi_step * stave;
  double stave_phi_offset =  M_PI /2.0;    // stave initially points so that sensor faces upward in y

  // Is this is an inner, middle or outer stave?
  if( stave_type == 0 )
    {
      // Start with the point in sensor local coords 
      TVector3   pos1 = sensor_local;
      cout << " Start with local coords in sensor:  pos1 = " << pos1.X() << " " << pos1.Y() << " " << pos1.Z() << endl;

      // transform sensor location to location in chip
      TVector3 tr1(loc_sensor_in_chip[0], loc_sensor_in_chip[1], loc_sensor_in_chip[2]);
      TVector3  res = pos1 + tr1;

      cout << " tr11 = " << tr1.X() << " " << tr1.Y() << " " << tr1.Z() << endl;
      cout << " Translated  to local coords in chip: res = " << res.X() << " " << res.Y() << " " << res.Z() << endl;
      
      // transform location in chip to location in module
      TVector3 tr2(inner_loc_chip_in_module[chip][0], inner_loc_chip_in_module[chip][1], inner_loc_chip_in_module[chip][2]);
      res = res + tr2;
      cout << " tr2 = " << tr2.X() << " " << tr2.Y() << " " << tr2.Z() << endl;
      cout << " Translated to local coords in module: res = " << res.X() << " " << res.Y() << " " << res.Z() << endl;

      // module to half stave
      TVector3 tr2a(inner_loc_module_in_halfstave[0], inner_loc_module_in_halfstave[1], inner_loc_module_in_halfstave[2] );
      res = res + tr2a;
      cout << " Translated to local coords in half stave: res = " << res.X() << " " << res.Y() << " " << res.Z() << endl;

      // transform location in half stave to location in stave 
      TVector3 tr3(inner_loc_halfstave_in_stave[0],  inner_loc_halfstave_in_stave[1],inner_loc_halfstave_in_stave[2]);
      res = res + tr3;
      cout << " tr3 = " << tr3.X() << " " << tr3.Y() << " " << tr3.Z() << endl;
      cout << " Translated to local coords in stave: res = " << res.X() << " " << res.Y() << " " << res.Z() << endl;

      // Rotate stave to its angle in the world
      // This requires rotating it by
      //    90 degrees to make it point to the origin instead of vertically up in y when it is at phi = 0
      //    Rotating it fiurther so that it points at the origin after being translated to the x and y coords of the stave phi location
      //    adding the tilt (for layers 0-2)
      // for a rotation
      TRotation R;
      R.RotateZ(stave_phi + stave_phi_offset + stave_phi_tilt);
      res = R * res;    // rotates res using R
      cout << "Rotate through phi = " << stave_phi << " + phi_offset = " << stave_phi_offset << " + phitilt = " << stave_phi_tilt << endl;
      cout << " Rotated stave to point at origin, then tilted it: res = " << res.X() << " " << res.Y() << " " << res.Z() << endl;

      // transform location of stave to its location in the world 
      TVector3 tr4(layer_radius * cos(stave_phi), layer_radius * sin(stave_phi), 0.0);
      res = res + tr4;
      cout << " tr4 = " << tr4.X() << " " << tr4.Y() << " " << tr4.Z() << endl;
      cout << " Translated stave to location in world: res = " << res.X() << " " << res.Y() << " " << res.Z() << endl;

      return res;
    }
  else if(stave_type == 1)
    {
      // Start with the point in sensor local coords 
      TVector3   pos1 = sensor_local;
      //cout << " Start with local coords in sensor:  pos1 = " << pos1.X() << " " << pos1.Y() << " " << pos1.Z() << endl;

      // transform sensor location to location in chip
      TVector3 tr1(loc_sensor_in_chip[0], loc_sensor_in_chip[1], loc_sensor_in_chip[2]);
      TVector3  res = pos1 + tr1;
      //cout << " tr11 = " << tr1.X() << " " << tr1.Y() << " " << tr1.Z() << endl;
      //cout << " Translated  to local coords in chip: res = " << res.X() << " " << res.Y() << " " << res.Z() << endl;
      
      // transform location in chip to location in module
      // for odd numbered chips, the chip is flipped by 180 degrees around the x and z axes
      TRotation Rchip;
      if(chip % 2)
	{
	  Rchip.RotateZ(M_PI);
	  Rchip.RotateX(M_PI);
	  res = Rchip * res;
	  //cout << "Rotaed 180- deg around Z and X axes, now res = "  << res.X() << " " << res.Y() << " " << res.Z() << endl;
	}

      TVector3 tr2(middle_loc_chip_in_module[chip][0], middle_loc_chip_in_module[chip][1], middle_loc_chip_in_module[chip][2]);
      res = res + tr2;
      //cout << " tr2 = " << tr2.X() << " " << tr2.Y() << " " << tr2.Z() << endl;
      //cout << " Translated to local coords in module: res = " << res.X() << " " << res.Y() << " " << res.Z() << endl;

      // module to half stave
      TVector3 tr2a(middle_loc_module_in_halfstave[module][0], middle_loc_module_in_halfstave[module][1], middle_loc_module_in_halfstave[module][2] );
      res = res + tr2a;
      //cout << " tr2a = " << tr2a.X() << " " << tr2a.Y() << " " << tr2a.Z() << endl;
      //cout << " Translated to local coords in half stave: res = " << res.X() << " " << res.Y() << " " << res.Z() << endl;

      // transform location in half stave to location in stave 
      TVector3 tr3(middle_loc_halfstave_in_stave[half_stave][0],  middle_loc_halfstave_in_stave[half_stave][1],middle_loc_halfstave_in_stave[half_stave][2]);
      res = res + tr3;
      //cout << " tr3 = " << tr3.X() << " " << tr3.Y() << " " << tr3.Z() << endl;
      //cout << " Translated to local coords in stave: res = " << res.X() << " " << res.Y() << " " << res.Z() << endl;

      // Rotate stave to its angle in the world
      // This requires rotating it by
      //    90 degrees to make it point to the origin instead of vertically up in y when it is at phi = 0
      //    Rotating it fiurther so that it points at the origin after being translated to the x and y coords of the stave phi location
      //    adding the tilt (for layers 0-2)
      // for a rotation
      TRotation R;
      R.RotateZ(stave_phi + stave_phi_offset + stave_phi_tilt);
      res = R * res;    // rotates res using R
      //cout << "Rotate through phi = " << stave_phi << " + phi_offset = " << stave_phi_offset << " + phitilt = " << stave_phi_tilt << endl;
      //cout << " Rotated stave to point at origin, then tilted it: res = " << res.X() << " " << res.Y() << " " << res.Z() << endl;

      // transform location of stave to its location in the world 
      TVector3 tr4(layer_radius * cos(stave_phi), layer_radius * sin(stave_phi), 0.0);
      res = res + tr4;
      //cout << " tr4 = " << tr4.X() << " " << tr4.Y() << " " << tr4.Z() << endl;
      //cout << " Translated stave to location in world: res = " << res.X() << " " << res.Y() << " " << res.Z() << endl;

      return res;
    }
  else
    {
      // stave_type = 2

      // Start with the point in sensor local coords 
      TVector3   pos1 = sensor_local;
      //cout << " Start with local coords in sensor:  pos1 = " << pos1.X() << " " << pos1.Y() << " " << pos1.Z() << endl;

      // transform sensor location to location in chip
      TVector3 tr1(loc_sensor_in_chip[0], loc_sensor_in_chip[1], loc_sensor_in_chip[2]);
      TVector3  res = pos1 + tr1;
      //cout << " tr11 = " << tr1.X() << " " << tr1.Y() << " " << tr1.Z() << endl;
      //cout << " Translated  to local coords in chip: res = " << res.X() << " " << res.Y() << " " << res.Z() << endl;
      
      // transform location in chip to location in module
      // for odd numbered chips, the chip is flipped by 180 degrees around the x and z axes
      TRotation Rchip;
      if(chip % 2)
	{
	  Rchip.RotateZ(M_PI);
	  Rchip.RotateX(M_PI);
	  res = Rchip * res;
	  //cout << "Rotaed 180- deg around Z and X axes, now res = "  << res.X() << " " << res.Y() << " " << res.Z() << endl;
	}

      TVector3 tr2(outer_loc_chip_in_module[chip][0], outer_loc_chip_in_module[chip][1], outer_loc_chip_in_module[chip][2]);
      res = res + tr2;
      //cout << " tr2 = " << tr2.X() << " " << tr2.Y() << " " << tr2.Z() << endl;
      //cout << " Translated to local coords in module: res = " << res.X() << " " << res.Y() << " " << res.Z() << endl;

      // module to half stave
      TVector3 tr2a(outer_loc_module_in_halfstave[module][0], outer_loc_module_in_halfstave[module][1], outer_loc_module_in_halfstave[module][2] );
      res = res + tr2a;
      //cout << " tr2a = " << tr2a.X() << " " << tr2a.Y() << " " << tr2a.Z() << endl;
      //cout << " Translated to local coords in half stave: res = " << res.X() << " " << res.Y() << " " << res.Z() << endl;

      // transform location in half stave to location in stave 
      TVector3 tr3(outer_loc_halfstave_in_stave[half_stave][0],  outer_loc_halfstave_in_stave[half_stave][1], outer_loc_halfstave_in_stave[half_stave][2]);
      res = res + tr3;
      //cout << " tr3 = " << tr3.X() << " " << tr3.Y() << " " << tr3.Z() << endl;
      //cout << " Translated to local coords in stave: res = " << res.X() << " " << res.Y() << " " << res.Z() << endl;

      // Rotate stave to its angle in the world
      // This requires rotating it by
      //    90 degrees to make it point to the origin instead of vertically up in y when it is at phi = 0
      //    Rotating it fiurther so that it points at the origin after being translated to the x and y coords of the stave phi location
      //    adding the tilt (for layers 0-2)
      // for a rotation
      TRotation R;
      R.RotateZ(stave_phi + stave_phi_offset + stave_phi_tilt);
      res = R * res;    // rotates res using R
      //cout << "Rotate through phi = " << stave_phi << " + phi_offset = " << stave_phi_offset << " + phitilt = " << stave_phi_tilt << endl;
      //cout << " Rotated stave to point at origin, then tilted it: res = " << res.X() << " " << res.Y() << " " << res.Z() << endl;

      // transform location of stave to its location in the world 
      TVector3 tr4(layer_radius * cos(stave_phi), layer_radius * sin(stave_phi), 0.0);
      res = res + tr4;
      //cout << " tr4 = " << tr4.X() << " " << tr4.Y() << " " << tr4.Z() << endl;
      //cout << " Translated stave to location in world: res = " << res.X() << " " << res.Y() << " " << res.Z() << endl;

      return res;
    }


}

void
PHG4CylinderGeom_MAPS::identify(std::ostream& os) const
{
  os << "PHG4CylinderGeom_MAPS: layer: " << layer 
     << ", layer_radius: " << layer_radius 
     << " , stave_type " << stave_type
     << ", N_staves in layer: " << N_staves
     << ", N_half_staves in layer: " << N_half_staves
     << ", N_modules in layer: " << N_modules
     << ", N_chips in layer: " << N_chips
     << endl;
  return;
}

void PHG4CylinderGeom_MAPS::find_sensor_center(int stave_number, int half_stave_number, int module_number, int chip_number, double location[])
{
  double x_location = 0.0;
  double y_location = 0.0;
  double z_location = 0.0;

  location[0] = x_location;
  location[1] = y_location;
  location[2] = z_location;
}


