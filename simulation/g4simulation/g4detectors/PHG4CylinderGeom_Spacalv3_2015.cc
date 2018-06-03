// $$Id: PHG4CylinderGeom_Spacalv3.cc,v 1.3 2014/08/28 22:18:35 jinhuang Exp $$

/*!
 * \file ${file_name}
 * \brief This file include obsolete design parameter of PHG4CylinderGeom_Spacalv3::load_demo_sector_tower_map_2015_Chris_Cullen_2D_spacal(). Only kept for comparison study purpose.
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $$Revision: 1.3 $$
 * \date $$Date: 2014/08/28 22:18:35 $$
 */

#include "PHG4CylinderGeom_Spacalv3.h"

#include <Geant4/globals.hh>
#include <Geant4/G4PhysicalConstants.hh>

#include <algorithm>
#include <cmath>
#include <cassert>
#include <iostream>
#include <sstream>
#include <limits>       // std::numeric_limits
#include <map>

using namespace std;
using std::make_pair;


void
PHG4CylinderGeom_Spacalv3::load_demo_sector_tower_map_2015_Chris_Cullen_2D_spacal()
{
  cout << "PHG4CylinderGeom_Spacalv3::load_demo_sector_tower_map_2015_Chris_Cullen_2D_spacal - "
      << "load Chris Cullen 2D spacal design July 2015" << endl;

  // Chris Cullen 2D spacal design July 2015
  radius = 90.000000;
  thickness = 26.130000;
  zmax = 149.470000;
  zmin = -zmax;

  azimuthal_n_sec = 32;
  max_phi_bin_in_sec = 8;
  init_default_sector_map();
//  cout << "PHG4CylinderGeom_Spacalv3::load_demo_sector_tower_map_2015_Chris_Cullen_2D_spacal - init_default_sector_map size to "
//      <<sector_map.size()<<"/"<<get_azimuthal_n_sec()<<endl;

  azimuthal_tilt = 0;
  azimuthal_seg_visible = false;
  polar_taper_ratio = 1.;
  assembly_spacing = 0.002500;
  sidewall_thickness = 0.075000;
  sidewall_outer_torr = 0.030000;
  sector_tower_map.clear();
    {
      // tower 1021 based Row/Col = 102/1
      geom_tower geom;
      geom.id = 1021;
      geom.pDz = 6.751948;
      geom.pTheta = 0.086675;
      geom.pPhi = -3.004692;
      geom.pAlp1 = -0.002046;
      geom.pAlp2 = -0.002046;
      geom.pDy1 = 1.121195;
      geom.pDx1 = 1.199084;
      geom.pDx2 = 1.199749;
      geom.pDy2 = 1.281451;
      geom.pDx3 = 1.365897;
      geom.pDx4 = 1.366657;
      geom.centralX = -8.954176;
      geom.centralY = 105.060369;
      geom.centralZ = 3.686651;
      geom.pRotationAngleX = -1.547057;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1022 based Row/Col = 102/2
      geom_tower geom;
      geom.id = 1022;
      geom.pDz = 6.751948;
      geom.pTheta = 0.062464;
      geom.pPhi = -2.950840;
      geom.pAlp1 = -0.001456;
      geom.pAlp2 = -0.001456;
      geom.pDy1 = 1.121195;
      geom.pDx1 = 1.194747;
      geom.pDx2 = 1.195405;
      geom.pDy2 = 1.281451;
      geom.pDx3 = 1.360958;
      geom.pDx4 = 1.361710;
      geom.centralX = -6.388124;
      geom.centralY = 105.060369;
      geom.centralZ = 3.686651;
      geom.pRotationAngleX = -1.547057;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1023 based Row/Col = 102/3
      geom_tower geom;
      geom.id = 1023;
      geom.pDz = 6.751948;
      geom.pTheta = 0.038660;
      geom.pPhi = -2.829992;
      geom.pAlp1 = -0.000872;
      geom.pAlp2 = -0.000872;
      geom.pDy1 = 1.121195;
      geom.pDx1 = 1.191868;
      geom.pDx2 = 1.192521;
      geom.pDy2 = 1.281451;
      geom.pDx3 = 1.357679;
      geom.pDx4 = 1.358425;
      geom.centralX = -3.829796;
      geom.centralY = 105.060369;
      geom.centralZ = 3.686651;
      geom.pRotationAngleX = -1.547057;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1024 based Row/Col = 102/4
      geom_tower geom;
      geom.id = 1024;
      geom.pDz = 6.751948;
      geom.pTheta = 0.017060;
      geom.pPhi = -2.373142;
      geom.pAlp1 = -0.000290;
      geom.pAlp2 = -0.000290;
      geom.pDy1 = 1.121195;
      geom.pDx1 = 1.190432;
      geom.pDx2 = 1.191082;
      geom.pDy2 = 1.281451;
      geom.pDx3 = 1.356043;
      geom.pDx4 = 1.356787;
      geom.centralX = -1.276086;
      geom.centralY = 105.060369;
      geom.centralZ = 3.686651;
      geom.pRotationAngleX = -1.547057;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1025 based Row/Col = 102/5
      geom_tower geom;
      geom.id = 1025;
      geom.pDz = 6.751948;
      geom.pTheta = 0.017060;
      geom.pPhi = -0.768451;
      geom.pAlp1 = 0.000290;
      geom.pAlp2 = 0.000290;
      geom.pDy1 = 1.121195;
      geom.pDx1 = 1.190432;
      geom.pDx2 = 1.191082;
      geom.pDy2 = 1.281451;
      geom.pDx3 = 1.356043;
      geom.pDx4 = 1.356787;
      geom.centralX = 1.276086;
      geom.centralY = 105.060369;
      geom.centralZ = 3.686651;
      geom.pRotationAngleX = -1.547057;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1026 based Row/Col = 102/6
      geom_tower geom;
      geom.id = 1026;
      geom.pDz = 6.751948;
      geom.pTheta = 0.038660;
      geom.pPhi = -0.311600;
      geom.pAlp1 = 0.000872;
      geom.pAlp2 = 0.000872;
      geom.pDy1 = 1.121195;
      geom.pDx1 = 1.191868;
      geom.pDx2 = 1.192521;
      geom.pDy2 = 1.281451;
      geom.pDx3 = 1.357679;
      geom.pDx4 = 1.358425;
      geom.centralX = 3.829796;
      geom.centralY = 105.060369;
      geom.centralZ = 3.686651;
      geom.pRotationAngleX = -1.547057;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1027 based Row/Col = 102/7
      geom_tower geom;
      geom.id = 1027;
      geom.pDz = 6.751948;
      geom.pTheta = 0.062464;
      geom.pPhi = -0.190753;
      geom.pAlp1 = 0.001456;
      geom.pAlp2 = 0.001456;
      geom.pDy1 = 1.121195;
      geom.pDx1 = 1.194747;
      geom.pDx2 = 1.195405;
      geom.pDy2 = 1.281451;
      geom.pDx3 = 1.360958;
      geom.pDx4 = 1.361710;
      geom.centralX = 6.388124;
      geom.centralY = 105.060369;
      geom.centralZ = 3.686651;
      geom.pRotationAngleX = -1.547057;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1028 based Row/Col = 102/8
      geom_tower geom;
      geom.id = 1028;
      geom.pDz = 6.751948;
      geom.pTheta = 0.086675;
      geom.pPhi = -0.136901;
      geom.pAlp1 = 0.002046;
      geom.pAlp2 = 0.002046;
      geom.pDy1 = 1.121195;
      geom.pDx1 = 1.199084;
      geom.pDx2 = 1.199749;
      geom.pDy2 = 1.281451;
      geom.pDx3 = 1.365897;
      geom.pDx4 = 1.366657;
      geom.centralX = 8.954176;
      geom.centralY = 105.060369;
      geom.centralZ = 3.686651;
      geom.pRotationAngleX = -1.547057;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1011 based Row/Col = 101/1
      geom_tower geom;
      geom.id = 1011;
      geom.pDz = 6.751954;
      geom.pTheta = 0.086725;
      geom.pPhi = 3.004566;
      geom.pAlp1 = -0.002046;
      geom.pAlp2 = -0.002046;
      geom.pDy1 = 1.121195;
      geom.pDx1 = 1.199749;
      geom.pDx2 = 1.200414;
      geom.pDy2 = 1.281451;
      geom.pDx3 = 1.366657;
      geom.pDx4 = 1.367417;
      geom.centralX = -8.959092;
      geom.centralY = 105.117519;
      geom.centralZ = 1.279684;
      geom.pRotationAngleX = -1.547057;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1012 based Row/Col = 101/2
      geom_tower geom;
      geom.id = 1012;
      geom.pDz = 6.751954;
      geom.pTheta = 0.062501;
      geom.pPhi = 2.950667;
      geom.pAlp1 = -0.001456;
      geom.pAlp2 = -0.001456;
      geom.pDy1 = 1.121195;
      geom.pDx1 = 1.195405;
      geom.pDx2 = 1.196063;
      geom.pDy2 = 1.281451;
      geom.pDx3 = 1.361710;
      geom.pDx4 = 1.362462;
      geom.centralX = -6.391623;
      geom.centralY = 105.117519;
      geom.centralZ = 1.279684;
      geom.pRotationAngleX = -1.547057;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1013 based Row/Col = 101/3
      geom_tower geom;
      geom.id = 1013;
      geom.pDz = 6.751954;
      geom.pTheta = 0.038686;
      geom.pPhi = 2.829721;
      geom.pAlp1 = -0.000872;
      geom.pAlp2 = -0.000872;
      geom.pDy1 = 1.121195;
      geom.pDx1 = 1.192521;
      geom.pDx2 = 1.193174;
      geom.pDy2 = 1.281451;
      geom.pDx3 = 1.358425;
      geom.pDx4 = 1.359172;
      geom.centralX = -3.831890;
      geom.centralY = 105.117519;
      geom.centralZ = 1.279684;
      geom.pRotationAngleX = -1.547057;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1014 based Row/Col = 101/4
      geom_tower geom;
      geom.id = 1014;
      geom.pDz = 6.751954;
      geom.pTheta = 0.017078;
      geom.pPhi = 2.372677;
      geom.pAlp1 = -0.000290;
      geom.pAlp2 = -0.000290;
      geom.pDy1 = 1.121195;
      geom.pDx1 = 1.191082;
      geom.pDx2 = 1.191733;
      geom.pDy2 = 1.281451;
      geom.pDx3 = 1.356787;
      geom.pDx4 = 1.357531;
      geom.centralX = -1.276783;
      geom.centralY = 105.117519;
      geom.centralZ = 1.279684;
      geom.pRotationAngleX = -1.547057;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1015 based Row/Col = 101/5
      geom_tower geom;
      geom.id = 1015;
      geom.pDz = 6.751954;
      geom.pTheta = 0.017078;
      geom.pPhi = 0.768916;
      geom.pAlp1 = 0.000290;
      geom.pAlp2 = 0.000290;
      geom.pDy1 = 1.121195;
      geom.pDx1 = 1.191082;
      geom.pDx2 = 1.191733;
      geom.pDy2 = 1.281451;
      geom.pDx3 = 1.356787;
      geom.pDx4 = 1.357531;
      geom.centralX = 1.276783;
      geom.centralY = 105.117519;
      geom.centralZ = 1.279684;
      geom.pRotationAngleX = -1.547057;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1016 based Row/Col = 101/6
      geom_tower geom;
      geom.id = 1016;
      geom.pDz = 6.751954;
      geom.pTheta = 0.038686;
      geom.pPhi = 0.311872;
      geom.pAlp1 = 0.000872;
      geom.pAlp2 = 0.000872;
      geom.pDy1 = 1.121195;
      geom.pDx1 = 1.192521;
      geom.pDx2 = 1.193174;
      geom.pDy2 = 1.281451;
      geom.pDx3 = 1.358425;
      geom.pDx4 = 1.359172;
      geom.centralX = 3.831890;
      geom.centralY = 105.117519;
      geom.centralZ = 1.279684;
      geom.pRotationAngleX = -1.547057;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1017 based Row/Col = 101/7
      geom_tower geom;
      geom.id = 1017;
      geom.pDz = 6.751954;
      geom.pTheta = 0.062501;
      geom.pPhi = 0.190926;
      geom.pAlp1 = 0.001456;
      geom.pAlp2 = 0.001456;
      geom.pDy1 = 1.121195;
      geom.pDx1 = 1.195405;
      geom.pDx2 = 1.196063;
      geom.pDy2 = 1.281451;
      geom.pDx3 = 1.361710;
      geom.pDx4 = 1.362462;
      geom.centralX = 6.391623;
      geom.centralY = 105.117519;
      geom.centralZ = 1.279684;
      geom.pRotationAngleX = -1.547057;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1018 based Row/Col = 101/8
      geom_tower geom;
      geom.id = 1018;
      geom.pDz = 6.751954;
      geom.pTheta = 0.086725;
      geom.pPhi = 0.137026;
      geom.pAlp1 = 0.002046;
      geom.pAlp2 = 0.002046;
      geom.pDy1 = 1.121195;
      geom.pDx1 = 1.199749;
      geom.pDx2 = 1.200414;
      geom.pDy2 = 1.281451;
      geom.pDx3 = 1.366657;
      geom.pDx4 = 1.367417;
      geom.centralX = 8.959092;
      geom.centralY = 105.117519;
      geom.centralZ = 1.279684;
      geom.pRotationAngleX = -1.547057;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 991 based Row/Col = 99/1
      geom_tower geom;
      geom.id = 991;
      geom.pDz = 6.751948;
      geom.pTheta = 0.086675;
      geom.pPhi = 3.004692;
      geom.pAlp1 = 0.002046;
      geom.pAlp2 = 0.002046;
      geom.pDy1 = 1.121195;
      geom.pDx1 = 1.199749;
      geom.pDx2 = 1.199084;
      geom.pDy2 = 1.281451;
      geom.pDx3 = 1.366657;
      geom.pDx4 = 1.365897;
      geom.centralX = -8.954176;
      geom.centralY = 105.060369;
      geom.centralZ = -3.686651;
      geom.pRotationAngleX = -1.594535;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 992 based Row/Col = 99/2
      geom_tower geom;
      geom.id = 992;
      geom.pDz = 6.751948;
      geom.pTheta = 0.062464;
      geom.pPhi = 2.950840;
      geom.pAlp1 = 0.001456;
      geom.pAlp2 = 0.001456;
      geom.pDy1 = 1.121195;
      geom.pDx1 = 1.195405;
      geom.pDx2 = 1.194747;
      geom.pDy2 = 1.281451;
      geom.pDx3 = 1.361710;
      geom.pDx4 = 1.360958;
      geom.centralX = -6.388124;
      geom.centralY = 105.060369;
      geom.centralZ = -3.686651;
      geom.pRotationAngleX = -1.594535;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 993 based Row/Col = 99/3
      geom_tower geom;
      geom.id = 993;
      geom.pDz = 6.751948;
      geom.pTheta = 0.038660;
      geom.pPhi = 2.829992;
      geom.pAlp1 = 0.000872;
      geom.pAlp2 = 0.000872;
      geom.pDy1 = 1.121195;
      geom.pDx1 = 1.192521;
      geom.pDx2 = 1.191868;
      geom.pDy2 = 1.281451;
      geom.pDx3 = 1.358425;
      geom.pDx4 = 1.357679;
      geom.centralX = -3.829796;
      geom.centralY = 105.060369;
      geom.centralZ = -3.686651;
      geom.pRotationAngleX = -1.594535;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 994 based Row/Col = 99/4
      geom_tower geom;
      geom.id = 994;
      geom.pDz = 6.751948;
      geom.pTheta = 0.017060;
      geom.pPhi = 2.373142;
      geom.pAlp1 = 0.000290;
      geom.pAlp2 = 0.000290;
      geom.pDy1 = 1.121195;
      geom.pDx1 = 1.191082;
      geom.pDx2 = 1.190432;
      geom.pDy2 = 1.281451;
      geom.pDx3 = 1.356787;
      geom.pDx4 = 1.356043;
      geom.centralX = -1.276086;
      geom.centralY = 105.060369;
      geom.centralZ = -3.686651;
      geom.pRotationAngleX = -1.594535;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 995 based Row/Col = 99/5
      geom_tower geom;
      geom.id = 995;
      geom.pDz = 6.751948;
      geom.pTheta = 0.017060;
      geom.pPhi = 0.768451;
      geom.pAlp1 = -0.000290;
      geom.pAlp2 = -0.000290;
      geom.pDy1 = 1.121195;
      geom.pDx1 = 1.191082;
      geom.pDx2 = 1.190432;
      geom.pDy2 = 1.281451;
      geom.pDx3 = 1.356787;
      geom.pDx4 = 1.356043;
      geom.centralX = 1.276086;
      geom.centralY = 105.060369;
      geom.centralZ = -3.686651;
      geom.pRotationAngleX = -1.594535;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 996 based Row/Col = 99/6
      geom_tower geom;
      geom.id = 996;
      geom.pDz = 6.751948;
      geom.pTheta = 0.038660;
      geom.pPhi = 0.311600;
      geom.pAlp1 = -0.000872;
      geom.pAlp2 = -0.000872;
      geom.pDy1 = 1.121195;
      geom.pDx1 = 1.192521;
      geom.pDx2 = 1.191868;
      geom.pDy2 = 1.281451;
      geom.pDx3 = 1.358425;
      geom.pDx4 = 1.357679;
      geom.centralX = 3.829796;
      geom.centralY = 105.060369;
      geom.centralZ = -3.686651;
      geom.pRotationAngleX = -1.594535;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 997 based Row/Col = 99/7
      geom_tower geom;
      geom.id = 997;
      geom.pDz = 6.751948;
      geom.pTheta = 0.062464;
      geom.pPhi = 0.190753;
      geom.pAlp1 = -0.001456;
      geom.pAlp2 = -0.001456;
      geom.pDy1 = 1.121195;
      geom.pDx1 = 1.195405;
      geom.pDx2 = 1.194747;
      geom.pDy2 = 1.281451;
      geom.pDx3 = 1.361710;
      geom.pDx4 = 1.360958;
      geom.centralX = 6.388124;
      geom.centralY = 105.060369;
      geom.centralZ = -3.686651;
      geom.pRotationAngleX = -1.594535;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 998 based Row/Col = 99/8
      geom_tower geom;
      geom.id = 998;
      geom.pDz = 6.751948;
      geom.pTheta = 0.086675;
      geom.pPhi = 0.136901;
      geom.pAlp1 = -0.002046;
      geom.pAlp2 = -0.002046;
      geom.pDy1 = 1.121195;
      geom.pDx1 = 1.199749;
      geom.pDx2 = 1.199084;
      geom.pDy2 = 1.281451;
      geom.pDx3 = 1.366657;
      geom.pDx4 = 1.365897;
      geom.centralX = 8.954176;
      geom.centralY = 105.060369;
      geom.centralZ = -3.686651;
      geom.pRotationAngleX = -1.594535;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1001 based Row/Col = 100/1
      geom_tower geom;
      geom.id = 1001;
      geom.pDz = 6.751954;
      geom.pTheta = 0.086725;
      geom.pPhi = -3.004566;
      geom.pAlp1 = 0.002046;
      geom.pAlp2 = 0.002046;
      geom.pDy1 = 1.121195;
      geom.pDx1 = 1.200414;
      geom.pDx2 = 1.199749;
      geom.pDy2 = 1.281451;
      geom.pDx3 = 1.367417;
      geom.pDx4 = 1.366657;
      geom.centralX = -8.959092;
      geom.centralY = 105.117519;
      geom.centralZ = -1.279684;
      geom.pRotationAngleX = -1.594535;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1002 based Row/Col = 100/2
      geom_tower geom;
      geom.id = 1002;
      geom.pDz = 6.751954;
      geom.pTheta = 0.062501;
      geom.pPhi = -2.950667;
      geom.pAlp1 = 0.001456;
      geom.pAlp2 = 0.001456;
      geom.pDy1 = 1.121195;
      geom.pDx1 = 1.196063;
      geom.pDx2 = 1.195405;
      geom.pDy2 = 1.281451;
      geom.pDx3 = 1.362462;
      geom.pDx4 = 1.361710;
      geom.centralX = -6.391623;
      geom.centralY = 105.117519;
      geom.centralZ = -1.279684;
      geom.pRotationAngleX = -1.594535;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1003 based Row/Col = 100/3
      geom_tower geom;
      geom.id = 1003;
      geom.pDz = 6.751954;
      geom.pTheta = 0.038686;
      geom.pPhi = -2.829721;
      geom.pAlp1 = 0.000872;
      geom.pAlp2 = 0.000872;
      geom.pDy1 = 1.121195;
      geom.pDx1 = 1.193174;
      geom.pDx2 = 1.192521;
      geom.pDy2 = 1.281451;
      geom.pDx3 = 1.359172;
      geom.pDx4 = 1.358425;
      geom.centralX = -3.831890;
      geom.centralY = 105.117519;
      geom.centralZ = -1.279684;
      geom.pRotationAngleX = -1.594535;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1004 based Row/Col = 100/4
      geom_tower geom;
      geom.id = 1004;
      geom.pDz = 6.751954;
      geom.pTheta = 0.017078;
      geom.pPhi = -2.372677;
      geom.pAlp1 = 0.000290;
      geom.pAlp2 = 0.000290;
      geom.pDy1 = 1.121195;
      geom.pDx1 = 1.191733;
      geom.pDx2 = 1.191082;
      geom.pDy2 = 1.281451;
      geom.pDx3 = 1.357531;
      geom.pDx4 = 1.356787;
      geom.centralX = -1.276783;
      geom.centralY = 105.117519;
      geom.centralZ = -1.279684;
      geom.pRotationAngleX = -1.594535;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1005 based Row/Col = 100/5
      geom_tower geom;
      geom.id = 1005;
      geom.pDz = 6.751954;
      geom.pTheta = 0.017078;
      geom.pPhi = -0.768916;
      geom.pAlp1 = -0.000290;
      geom.pAlp2 = -0.000290;
      geom.pDy1 = 1.121195;
      geom.pDx1 = 1.191733;
      geom.pDx2 = 1.191082;
      geom.pDy2 = 1.281451;
      geom.pDx3 = 1.357531;
      geom.pDx4 = 1.356787;
      geom.centralX = 1.276783;
      geom.centralY = 105.117519;
      geom.centralZ = -1.279684;
      geom.pRotationAngleX = -1.594535;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1006 based Row/Col = 100/6
      geom_tower geom;
      geom.id = 1006;
      geom.pDz = 6.751954;
      geom.pTheta = 0.038686;
      geom.pPhi = -0.311872;
      geom.pAlp1 = -0.000872;
      geom.pAlp2 = -0.000872;
      geom.pDy1 = 1.121195;
      geom.pDx1 = 1.193174;
      geom.pDx2 = 1.192521;
      geom.pDy2 = 1.281451;
      geom.pDx3 = 1.359172;
      geom.pDx4 = 1.358425;
      geom.centralX = 3.831890;
      geom.centralY = 105.117519;
      geom.centralZ = -1.279684;
      geom.pRotationAngleX = -1.594535;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1007 based Row/Col = 100/7
      geom_tower geom;
      geom.id = 1007;
      geom.pDz = 6.751954;
      geom.pTheta = 0.062501;
      geom.pPhi = -0.190926;
      geom.pAlp1 = -0.001456;
      geom.pAlp2 = -0.001456;
      geom.pDy1 = 1.121195;
      geom.pDx1 = 1.196063;
      geom.pDx2 = 1.195405;
      geom.pDy2 = 1.281451;
      geom.pDx3 = 1.362462;
      geom.pDx4 = 1.361710;
      geom.centralX = 6.391623;
      geom.centralY = 105.117519;
      geom.centralZ = -1.279684;
      geom.pRotationAngleX = -1.594535;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1008 based Row/Col = 100/8
      geom_tower geom;
      geom.id = 1008;
      geom.pDz = 6.751954;
      geom.pTheta = 0.086725;
      geom.pPhi = -0.137026;
      geom.pAlp1 = -0.002046;
      geom.pAlp2 = -0.002046;
      geom.pDy1 = 1.121195;
      geom.pDx1 = 1.200414;
      geom.pDx2 = 1.199749;
      geom.pDy2 = 1.281451;
      geom.pDx3 = 1.367417;
      geom.pDx4 = 1.366657;
      geom.centralX = 8.959092;
      geom.centralY = 105.117519;
      geom.centralZ = -1.279684;
      geom.pRotationAngleX = -1.594535;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1061 based Row/Col = 106/1
      geom_tower geom;
      geom.id = 1061;
      geom.pDz = 6.751959;
      geom.pTheta = 0.086041;
      geom.pPhi = -3.004366;
      geom.pAlp1 = -0.009890;
      geom.pAlp2 = -0.009890;
      geom.pDy1 = 1.121760;
      geom.pDx1 = 1.194637;
      geom.pDx2 = 1.197853;
      geom.pDy2 = 1.280866;
      geom.pDx3 = 1.360025;
      geom.pDx4 = 1.363697;
      geom.centralX = -8.928381;
      geom.centralY = 104.759951;
      geom.centralZ = 13.309949;
      geom.pRotationAngleX = -1.455814;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1062 based Row/Col = 106/2
      geom_tower geom;
      geom.id = 1062;
      geom.pDz = 6.751959;
      geom.pTheta = 0.062010;
      geom.pPhi = -2.950394;
      geom.pAlp1 = -0.007039;
      geom.pAlp2 = -0.007039;
      geom.pDy1 = 1.121760;
      geom.pDx1 = 1.190389;
      geom.pDx2 = 1.193571;
      geom.pDy2 = 1.280866;
      geom.pDx3 = 1.355191;
      geom.pDx4 = 1.358823;
      geom.centralX = -6.369834;
      geom.centralY = 104.759951;
      geom.centralZ = 13.309949;
      geom.pRotationAngleX = -1.455814;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1063 based Row/Col = 106/3
      geom_tower geom;
      geom.id = 1063;
      geom.pDz = 6.751959;
      geom.pTheta = 0.038385;
      geom.pPhi = -2.829297;
      geom.pAlp1 = -0.004213;
      geom.pAlp2 = -0.004213;
      geom.pDy1 = 1.121760;
      geom.pDx1 = 1.187569;
      geom.pDx2 = 1.190727;
      geom.pDy2 = 1.280866;
      geom.pDx3 = 1.351980;
      geom.pDx4 = 1.355586;
      geom.centralX = -3.818875;
      geom.centralY = 104.759951;
      geom.centralZ = 13.309949;
      geom.pRotationAngleX = -1.455814;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1064 based Row/Col = 106/4
      geom_tower geom;
      geom.id = 1064;
      geom.pDz = 6.751959;
      geom.pTheta = 0.016954;
      geom.pPhi = -2.371956;
      geom.pAlp1 = -0.001403;
      geom.pAlp2 = -0.001403;
      geom.pDy1 = 1.121760;
      geom.pDx1 = 1.186161;
      geom.pDx2 = 1.189308;
      geom.pDy2 = 1.280866;
      geom.pDx3 = 1.350379;
      geom.pDx4 = 1.353972;
      geom.centralX = -1.272455;
      geom.centralY = 104.759951;
      geom.centralZ = 13.309949;
      geom.pRotationAngleX = -1.455814;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1065 based Row/Col = 106/5
      geom_tower geom;
      geom.id = 1065;
      geom.pDz = 6.751959;
      geom.pTheta = 0.016954;
      geom.pPhi = -0.769637;
      geom.pAlp1 = 0.001403;
      geom.pAlp2 = 0.001403;
      geom.pDy1 = 1.121760;
      geom.pDx1 = 1.186161;
      geom.pDx2 = 1.189308;
      geom.pDy2 = 1.280866;
      geom.pDx3 = 1.350379;
      geom.pDx4 = 1.353972;
      geom.centralX = 1.272455;
      geom.centralY = 104.759951;
      geom.centralZ = 13.309949;
      geom.pRotationAngleX = -1.455814;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1066 based Row/Col = 106/6
      geom_tower geom;
      geom.id = 1066;
      geom.pDz = 6.751959;
      geom.pTheta = 0.038385;
      geom.pPhi = -0.312295;
      geom.pAlp1 = 0.004213;
      geom.pAlp2 = 0.004213;
      geom.pDy1 = 1.121760;
      geom.pDx1 = 1.187569;
      geom.pDx2 = 1.190727;
      geom.pDy2 = 1.280866;
      geom.pDx3 = 1.351980;
      geom.pDx4 = 1.355586;
      geom.centralX = 3.818875;
      geom.centralY = 104.759951;
      geom.centralZ = 13.309949;
      geom.pRotationAngleX = -1.455814;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1067 based Row/Col = 106/7
      geom_tower geom;
      geom.id = 1067;
      geom.pDz = 6.751959;
      geom.pTheta = 0.062010;
      geom.pPhi = -0.191198;
      geom.pAlp1 = 0.007039;
      geom.pAlp2 = 0.007039;
      geom.pDy1 = 1.121760;
      geom.pDx1 = 1.190389;
      geom.pDx2 = 1.193571;
      geom.pDy2 = 1.280866;
      geom.pDx3 = 1.355191;
      geom.pDx4 = 1.358823;
      geom.centralX = 6.369834;
      geom.centralY = 104.759951;
      geom.centralZ = 13.309949;
      geom.pRotationAngleX = -1.455814;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1068 based Row/Col = 106/8
      geom_tower geom;
      geom.id = 1068;
      geom.pDz = 6.751959;
      geom.pTheta = 0.086041;
      geom.pPhi = -0.137227;
      geom.pAlp1 = 0.009890;
      geom.pAlp2 = 0.009890;
      geom.pDy1 = 1.121760;
      geom.pDx1 = 1.194637;
      geom.pDx2 = 1.197853;
      geom.pDy2 = 1.280866;
      geom.pDx3 = 1.360025;
      geom.pDx4 = 1.363697;
      geom.centralX = 8.928381;
      geom.centralY = 104.759951;
      geom.centralZ = 13.309949;
      geom.pRotationAngleX = -1.455814;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1051 based Row/Col = 105/1
      geom_tower geom;
      geom.id = 1051;
      geom.pDz = 6.751993;
      geom.pTheta = 0.086265;
      geom.pPhi = 3.005126;
      geom.pAlp1 = -0.009890;
      geom.pAlp2 = -0.009889;
      geom.pDy1 = 1.121760;
      geom.pDx1 = 1.197853;
      geom.pDx2 = 1.201069;
      geom.pDy2 = 1.280866;
      geom.pDx3 = 1.363697;
      geom.pDx4 = 1.367369;
      geom.centralX = -8.952143;
      geom.centralY = 105.036176;
      geom.centralZ = 10.918222;
      geom.pRotationAngleX = -1.455814;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1052 based Row/Col = 105/2
      geom_tower geom;
      geom.id = 1052;
      geom.pDz = 6.751993;
      geom.pTheta = 0.062166;
      geom.pPhi = 2.951439;
      geom.pAlp1 = -0.007038;
      geom.pAlp2 = -0.007038;
      geom.pDy1 = 1.121760;
      geom.pDx1 = 1.193571;
      geom.pDx2 = 1.196752;
      geom.pDy2 = 1.280866;
      geom.pDx3 = 1.358823;
      geom.pDx4 = 1.362455;
      geom.centralX = -6.386746;
      geom.centralY = 105.036176;
      geom.centralZ = 10.918222;
      geom.pRotationAngleX = -1.455814;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1053 based Row/Col = 105/3
      geom_tower geom;
      geom.id = 1053;
      geom.pDz = 6.751993;
      geom.pTheta = 0.038469;
      geom.pPhi = 2.830934;
      geom.pAlp1 = -0.004213;
      geom.pAlp2 = -0.004213;
      geom.pDy1 = 1.121760;
      geom.pDx1 = 1.190727;
      geom.pDx2 = 1.193885;
      geom.pDy2 = 1.280866;
      geom.pDx3 = 1.355586;
      geom.pDx4 = 1.359192;
      geom.centralX = -3.828998;
      geom.centralY = 105.036176;
      geom.centralZ = 10.918222;
      geom.pRotationAngleX = -1.455814;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1054 based Row/Col = 105/4
      geom_tower geom;
      geom.id = 1054;
      geom.pDz = 6.751993;
      geom.pTheta = 0.016954;
      geom.pPhi = 2.374759;
      geom.pAlp1 = -0.001403;
      geom.pAlp2 = -0.001403;
      geom.pDy1 = 1.121760;
      geom.pDx1 = 1.189308;
      geom.pDx2 = 1.192455;
      geom.pDy2 = 1.280866;
      geom.pDx3 = 1.353972;
      geom.pDx4 = 1.357565;
      geom.centralX = -1.275825;
      geom.centralY = 105.036176;
      geom.centralZ = 10.918222;
      geom.pRotationAngleX = -1.455814;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1055 based Row/Col = 105/5
      geom_tower geom;
      geom.id = 1055;
      geom.pDz = 6.751993;
      geom.pTheta = 0.016954;
      geom.pPhi = 0.766834;
      geom.pAlp1 = 0.001403;
      geom.pAlp2 = 0.001403;
      geom.pDy1 = 1.121760;
      geom.pDx1 = 1.189308;
      geom.pDx2 = 1.192455;
      geom.pDy2 = 1.280866;
      geom.pDx3 = 1.353972;
      geom.pDx4 = 1.357565;
      geom.centralX = 1.275825;
      geom.centralY = 105.036176;
      geom.centralZ = 10.918222;
      geom.pRotationAngleX = -1.455814;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1056 based Row/Col = 105/6
      geom_tower geom;
      geom.id = 1056;
      geom.pDz = 6.751993;
      geom.pTheta = 0.038469;
      geom.pPhi = 0.310658;
      geom.pAlp1 = 0.004213;
      geom.pAlp2 = 0.004213;
      geom.pDy1 = 1.121760;
      geom.pDx1 = 1.190727;
      geom.pDx2 = 1.193885;
      geom.pDy2 = 1.280866;
      geom.pDx3 = 1.355586;
      geom.pDx4 = 1.359192;
      geom.centralX = 3.828998;
      geom.centralY = 105.036176;
      geom.centralZ = 10.918222;
      geom.pRotationAngleX = -1.455814;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1057 based Row/Col = 105/7
      geom_tower geom;
      geom.id = 1057;
      geom.pDz = 6.751993;
      geom.pTheta = 0.062166;
      geom.pPhi = 0.190153;
      geom.pAlp1 = 0.007038;
      geom.pAlp2 = 0.007038;
      geom.pDy1 = 1.121760;
      geom.pDx1 = 1.193571;
      geom.pDx2 = 1.196752;
      geom.pDy2 = 1.280866;
      geom.pDx3 = 1.358823;
      geom.pDx4 = 1.362455;
      geom.centralX = 6.386746;
      geom.centralY = 105.036176;
      geom.centralZ = 10.918222;
      geom.pRotationAngleX = -1.455814;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1058 based Row/Col = 105/8
      geom_tower geom;
      geom.id = 1058;
      geom.pDz = 6.751993;
      geom.pTheta = 0.086265;
      geom.pPhi = 0.136467;
      geom.pAlp1 = 0.009890;
      geom.pAlp2 = 0.009889;
      geom.pDy1 = 1.121760;
      geom.pDx1 = 1.197853;
      geom.pDx2 = 1.201069;
      geom.pDy2 = 1.280866;
      geom.pDx3 = 1.363697;
      geom.pDx4 = 1.367369;
      geom.centralX = 8.952143;
      geom.centralY = 105.036176;
      geom.centralZ = 10.918222;
      geom.pRotationAngleX = -1.455814;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 951 based Row/Col = 95/1
      geom_tower geom;
      geom.id = 951;
      geom.pDz = 6.751959;
      geom.pTheta = 0.086041;
      geom.pPhi = 3.004366;
      geom.pAlp1 = 0.009890;
      geom.pAlp2 = 0.009890;
      geom.pDy1 = 1.121760;
      geom.pDx1 = 1.197853;
      geom.pDx2 = 1.194637;
      geom.pDy2 = 1.280866;
      geom.pDx3 = 1.363697;
      geom.pDx4 = 1.360025;
      geom.centralX = -8.928381;
      geom.centralY = 104.759951;
      geom.centralZ = -13.309949;
      geom.pRotationAngleX = -1.685779;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 952 based Row/Col = 95/2
      geom_tower geom;
      geom.id = 952;
      geom.pDz = 6.751959;
      geom.pTheta = 0.062010;
      geom.pPhi = 2.950394;
      geom.pAlp1 = 0.007039;
      geom.pAlp2 = 0.007039;
      geom.pDy1 = 1.121760;
      geom.pDx1 = 1.193571;
      geom.pDx2 = 1.190389;
      geom.pDy2 = 1.280866;
      geom.pDx3 = 1.358823;
      geom.pDx4 = 1.355191;
      geom.centralX = -6.369834;
      geom.centralY = 104.759951;
      geom.centralZ = -13.309949;
      geom.pRotationAngleX = -1.685779;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 953 based Row/Col = 95/3
      geom_tower geom;
      geom.id = 953;
      geom.pDz = 6.751959;
      geom.pTheta = 0.038385;
      geom.pPhi = 2.829297;
      geom.pAlp1 = 0.004213;
      geom.pAlp2 = 0.004213;
      geom.pDy1 = 1.121760;
      geom.pDx1 = 1.190727;
      geom.pDx2 = 1.187569;
      geom.pDy2 = 1.280866;
      geom.pDx3 = 1.355586;
      geom.pDx4 = 1.351980;
      geom.centralX = -3.818875;
      geom.centralY = 104.759951;
      geom.centralZ = -13.309949;
      geom.pRotationAngleX = -1.685779;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 954 based Row/Col = 95/4
      geom_tower geom;
      geom.id = 954;
      geom.pDz = 6.751959;
      geom.pTheta = 0.016954;
      geom.pPhi = 2.371956;
      geom.pAlp1 = 0.001403;
      geom.pAlp2 = 0.001403;
      geom.pDy1 = 1.121760;
      geom.pDx1 = 1.189308;
      geom.pDx2 = 1.186161;
      geom.pDy2 = 1.280866;
      geom.pDx3 = 1.353972;
      geom.pDx4 = 1.350379;
      geom.centralX = -1.272455;
      geom.centralY = 104.759951;
      geom.centralZ = -13.309949;
      geom.pRotationAngleX = -1.685779;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 955 based Row/Col = 95/5
      geom_tower geom;
      geom.id = 955;
      geom.pDz = 6.751959;
      geom.pTheta = 0.016954;
      geom.pPhi = 0.769637;
      geom.pAlp1 = -0.001403;
      geom.pAlp2 = -0.001403;
      geom.pDy1 = 1.121760;
      geom.pDx1 = 1.189308;
      geom.pDx2 = 1.186161;
      geom.pDy2 = 1.280866;
      geom.pDx3 = 1.353972;
      geom.pDx4 = 1.350379;
      geom.centralX = 1.272455;
      geom.centralY = 104.759951;
      geom.centralZ = -13.309949;
      geom.pRotationAngleX = -1.685779;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 956 based Row/Col = 95/6
      geom_tower geom;
      geom.id = 956;
      geom.pDz = 6.751959;
      geom.pTheta = 0.038385;
      geom.pPhi = 0.312295;
      geom.pAlp1 = -0.004213;
      geom.pAlp2 = -0.004213;
      geom.pDy1 = 1.121760;
      geom.pDx1 = 1.190727;
      geom.pDx2 = 1.187569;
      geom.pDy2 = 1.280866;
      geom.pDx3 = 1.355586;
      geom.pDx4 = 1.351980;
      geom.centralX = 3.818875;
      geom.centralY = 104.759951;
      geom.centralZ = -13.309949;
      geom.pRotationAngleX = -1.685779;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 957 based Row/Col = 95/7
      geom_tower geom;
      geom.id = 957;
      geom.pDz = 6.751959;
      geom.pTheta = 0.062010;
      geom.pPhi = 0.191198;
      geom.pAlp1 = -0.007039;
      geom.pAlp2 = -0.007039;
      geom.pDy1 = 1.121760;
      geom.pDx1 = 1.193571;
      geom.pDx2 = 1.190389;
      geom.pDy2 = 1.280866;
      geom.pDx3 = 1.358823;
      geom.pDx4 = 1.355191;
      geom.centralX = 6.369834;
      geom.centralY = 104.759951;
      geom.centralZ = -13.309949;
      geom.pRotationAngleX = -1.685779;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 958 based Row/Col = 95/8
      geom_tower geom;
      geom.id = 958;
      geom.pDz = 6.751959;
      geom.pTheta = 0.086041;
      geom.pPhi = 0.137227;
      geom.pAlp1 = -0.009890;
      geom.pAlp2 = -0.009890;
      geom.pDy1 = 1.121760;
      geom.pDx1 = 1.197853;
      geom.pDx2 = 1.194637;
      geom.pDy2 = 1.280866;
      geom.pDx3 = 1.363697;
      geom.pDx4 = 1.360025;
      geom.centralX = 8.928381;
      geom.centralY = 104.759951;
      geom.centralZ = -13.309949;
      geom.pRotationAngleX = -1.685779;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 961 based Row/Col = 96/1
      geom_tower geom;
      geom.id = 961;
      geom.pDz = 6.751993;
      geom.pTheta = 0.086265;
      geom.pPhi = -3.005126;
      geom.pAlp1 = 0.009890;
      geom.pAlp2 = 0.009889;
      geom.pDy1 = 1.121760;
      geom.pDx1 = 1.201069;
      geom.pDx2 = 1.197853;
      geom.pDy2 = 1.280866;
      geom.pDx3 = 1.367369;
      geom.pDx4 = 1.363697;
      geom.centralX = -8.952143;
      geom.centralY = 105.036176;
      geom.centralZ = -10.918222;
      geom.pRotationAngleX = -1.685779;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 962 based Row/Col = 96/2
      geom_tower geom;
      geom.id = 962;
      geom.pDz = 6.751993;
      geom.pTheta = 0.062166;
      geom.pPhi = -2.951439;
      geom.pAlp1 = 0.007038;
      geom.pAlp2 = 0.007038;
      geom.pDy1 = 1.121760;
      geom.pDx1 = 1.196752;
      geom.pDx2 = 1.193571;
      geom.pDy2 = 1.280866;
      geom.pDx3 = 1.362455;
      geom.pDx4 = 1.358823;
      geom.centralX = -6.386746;
      geom.centralY = 105.036176;
      geom.centralZ = -10.918222;
      geom.pRotationAngleX = -1.685779;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 963 based Row/Col = 96/3
      geom_tower geom;
      geom.id = 963;
      geom.pDz = 6.751993;
      geom.pTheta = 0.038469;
      geom.pPhi = -2.830934;
      geom.pAlp1 = 0.004213;
      geom.pAlp2 = 0.004213;
      geom.pDy1 = 1.121760;
      geom.pDx1 = 1.193885;
      geom.pDx2 = 1.190727;
      geom.pDy2 = 1.280866;
      geom.pDx3 = 1.359192;
      geom.pDx4 = 1.355586;
      geom.centralX = -3.828998;
      geom.centralY = 105.036176;
      geom.centralZ = -10.918222;
      geom.pRotationAngleX = -1.685779;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 964 based Row/Col = 96/4
      geom_tower geom;
      geom.id = 964;
      geom.pDz = 6.751993;
      geom.pTheta = 0.016954;
      geom.pPhi = -2.374759;
      geom.pAlp1 = 0.001403;
      geom.pAlp2 = 0.001403;
      geom.pDy1 = 1.121760;
      geom.pDx1 = 1.192455;
      geom.pDx2 = 1.189308;
      geom.pDy2 = 1.280866;
      geom.pDx3 = 1.357565;
      geom.pDx4 = 1.353972;
      geom.centralX = -1.275825;
      geom.centralY = 105.036176;
      geom.centralZ = -10.918222;
      geom.pRotationAngleX = -1.685779;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 965 based Row/Col = 96/5
      geom_tower geom;
      geom.id = 965;
      geom.pDz = 6.751993;
      geom.pTheta = 0.016954;
      geom.pPhi = -0.766834;
      geom.pAlp1 = -0.001403;
      geom.pAlp2 = -0.001403;
      geom.pDy1 = 1.121760;
      geom.pDx1 = 1.192455;
      geom.pDx2 = 1.189308;
      geom.pDy2 = 1.280866;
      geom.pDx3 = 1.357565;
      geom.pDx4 = 1.353972;
      geom.centralX = 1.275825;
      geom.centralY = 105.036176;
      geom.centralZ = -10.918222;
      geom.pRotationAngleX = -1.685779;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 966 based Row/Col = 96/6
      geom_tower geom;
      geom.id = 966;
      geom.pDz = 6.751993;
      geom.pTheta = 0.038469;
      geom.pPhi = -0.310658;
      geom.pAlp1 = -0.004213;
      geom.pAlp2 = -0.004213;
      geom.pDy1 = 1.121760;
      geom.pDx1 = 1.193885;
      geom.pDx2 = 1.190727;
      geom.pDy2 = 1.280866;
      geom.pDx3 = 1.359192;
      geom.pDx4 = 1.355586;
      geom.centralX = 3.828998;
      geom.centralY = 105.036176;
      geom.centralZ = -10.918222;
      geom.pRotationAngleX = -1.685779;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 967 based Row/Col = 96/7
      geom_tower geom;
      geom.id = 967;
      geom.pDz = 6.751993;
      geom.pTheta = 0.062166;
      geom.pPhi = -0.190153;
      geom.pAlp1 = -0.007038;
      geom.pAlp2 = -0.007038;
      geom.pDy1 = 1.121760;
      geom.pDx1 = 1.196752;
      geom.pDx2 = 1.193571;
      geom.pDy2 = 1.280866;
      geom.pDx3 = 1.362455;
      geom.pDx4 = 1.358823;
      geom.centralX = 6.386746;
      geom.centralY = 105.036176;
      geom.centralZ = -10.918222;
      geom.pRotationAngleX = -1.685779;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 968 based Row/Col = 96/8
      geom_tower geom;
      geom.id = 968;
      geom.pDz = 6.751993;
      geom.pTheta = 0.086265;
      geom.pPhi = -0.136467;
      geom.pAlp1 = -0.009890;
      geom.pAlp2 = -0.009889;
      geom.pDy1 = 1.121760;
      geom.pDx1 = 1.201069;
      geom.pDx2 = 1.197853;
      geom.pDy2 = 1.280866;
      geom.pDx3 = 1.367369;
      geom.pDx4 = 1.363697;
      geom.centralX = 8.952143;
      geom.centralY = 105.036176;
      geom.centralZ = -10.918222;
      geom.pRotationAngleX = -1.685779;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1101 based Row/Col = 110/1
      geom_tower geom;
      geom.id = 1101;
      geom.pDz = 6.751953;
      geom.pTheta = 0.084718;
      geom.pPhi = -3.003391;
      geom.pAlp1 = -0.017594;
      geom.pAlp2 = -0.017591;
      geom.pDy1 = 1.122326;
      geom.pDx1 = 1.191567;
      geom.pDx2 = 1.197289;
      geom.pDy2 = 1.280215;
      geom.pDx3 = 1.354176;
      geom.pDx4 = 1.360702;
      geom.centralX = -8.907948;
      geom.centralY = 104.520810;
      geom.centralZ = 23.016845;
      geom.pRotationAngleX = -1.365259;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1102 based Row/Col = 110/2
      geom_tower geom;
      geom.id = 1102;
      geom.pDz = 6.751953;
      geom.pTheta = 0.061064;
      geom.pPhi = -2.949059;
      geom.pAlp1 = -0.012523;
      geom.pAlp2 = -0.012522;
      geom.pDy1 = 1.122326;
      geom.pDx1 = 1.187469;
      geom.pDx2 = 1.193131;
      geom.pDy2 = 1.280215;
      geom.pDx3 = 1.349520;
      geom.pDx4 = 1.355978;
      geom.centralX = -6.355490;
      geom.centralY = 104.520810;
      geom.centralZ = 23.016845;
      geom.pRotationAngleX = -1.365259;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1103 based Row/Col = 110/3
      geom_tower geom;
      geom.id = 1103;
      geom.pDz = 6.751953;
      geom.pTheta = 0.037815;
      geom.pPhi = -2.827213;
      geom.pAlp1 = -0.007497;
      geom.pAlp2 = -0.007496;
      geom.pDy1 = 1.122326;
      geom.pDx1 = 1.184748;
      geom.pDx2 = 1.190370;
      geom.pDy2 = 1.280215;
      geom.pDx3 = 1.346428;
      geom.pDx4 = 1.352841;
      geom.centralX = -3.810368;
      geom.centralY = 104.520810;
      geom.centralZ = 23.016845;
      geom.pRotationAngleX = -1.365259;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1104 based Row/Col = 110/4
      geom_tower geom;
      geom.id = 1104;
      geom.pDz = 6.751953;
      geom.pTheta = 0.016749;
      geom.pPhi = -2.368410;
      geom.pAlp1 = -0.002496;
      geom.pAlp2 = -0.002496;
      geom.pDy1 = 1.122326;
      geom.pDx1 = 1.183390;
      geom.pDx2 = 1.188993;
      geom.pDy2 = 1.280215;
      geom.pDx3 = 1.344885;
      geom.pDx4 = 1.351276;
      geom.centralX = -1.269636;
      geom.centralY = 104.520810;
      geom.centralZ = 23.016845;
      geom.pRotationAngleX = -1.365259;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1105 based Row/Col = 110/5
      geom_tower geom;
      geom.id = 1105;
      geom.pDz = 6.751953;
      geom.pTheta = 0.016749;
      geom.pPhi = -0.773183;
      geom.pAlp1 = 0.002496;
      geom.pAlp2 = 0.002496;
      geom.pDy1 = 1.122326;
      geom.pDx1 = 1.183390;
      geom.pDx2 = 1.188993;
      geom.pDy2 = 1.280215;
      geom.pDx3 = 1.344885;
      geom.pDx4 = 1.351276;
      geom.centralX = 1.269636;
      geom.centralY = 104.520810;
      geom.centralZ = 23.016845;
      geom.pRotationAngleX = -1.365259;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1106 based Row/Col = 110/6
      geom_tower geom;
      geom.id = 1106;
      geom.pDz = 6.751953;
      geom.pTheta = 0.037815;
      geom.pPhi = -0.314379;
      geom.pAlp1 = 0.007497;
      geom.pAlp2 = 0.007496;
      geom.pDy1 = 1.122326;
      geom.pDx1 = 1.184748;
      geom.pDx2 = 1.190370;
      geom.pDy2 = 1.280215;
      geom.pDx3 = 1.346428;
      geom.pDx4 = 1.352841;
      geom.centralX = 3.810368;
      geom.centralY = 104.520810;
      geom.centralZ = 23.016845;
      geom.pRotationAngleX = -1.365259;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1107 based Row/Col = 110/7
      geom_tower geom;
      geom.id = 1107;
      geom.pDz = 6.751953;
      geom.pTheta = 0.061064;
      geom.pPhi = -0.192533;
      geom.pAlp1 = 0.012523;
      geom.pAlp2 = 0.012522;
      geom.pDy1 = 1.122326;
      geom.pDx1 = 1.187469;
      geom.pDx2 = 1.193131;
      geom.pDy2 = 1.280215;
      geom.pDx3 = 1.349520;
      geom.pDx4 = 1.355978;
      geom.centralX = 6.355490;
      geom.centralY = 104.520810;
      geom.centralZ = 23.016845;
      geom.pRotationAngleX = -1.365259;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1108 based Row/Col = 110/8
      geom_tower geom;
      geom.id = 1108;
      geom.pDz = 6.751953;
      geom.pTheta = 0.084718;
      geom.pPhi = -0.138202;
      geom.pAlp1 = 0.017594;
      geom.pAlp2 = 0.017591;
      geom.pDy1 = 1.122326;
      geom.pDx1 = 1.191567;
      geom.pDx2 = 1.197289;
      geom.pDy2 = 1.280215;
      geom.pDx3 = 1.354176;
      geom.pDx4 = 1.360702;
      geom.centralX = 8.907948;
      geom.centralY = 104.520810;
      geom.centralZ = 23.016845;
      geom.pRotationAngleX = -1.365259;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1091 based Row/Col = 109/1
      geom_tower geom;
      geom.id = 1091;
      geom.pDz = 6.751988;
      geom.pTheta = 0.085120;
      geom.pPhi = 3.004213;
      geom.pAlp1 = -0.017593;
      geom.pAlp2 = -0.017591;
      geom.pDy1 = 1.122326;
      geom.pDx1 = 1.197289;
      geom.pDx2 = 1.203012;
      geom.pDy2 = 1.280215;
      geom.pDx3 = 1.360702;
      geom.pDx4 = 1.367228;
      geom.centralX = -8.950218;
      geom.centralY = 105.012173;
      geom.centralZ = 20.659979;
      geom.pRotationAngleX = -1.365259;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1092 based Row/Col = 109/2
      geom_tower geom;
      geom.id = 1092;
      geom.pDz = 6.751988;
      geom.pTheta = 0.061348;
      geom.pPhi = 2.950189;
      geom.pAlp1 = -0.012523;
      geom.pAlp2 = -0.012521;
      geom.pDy1 = 1.122326;
      geom.pDx1 = 1.193131;
      geom.pDx2 = 1.198793;
      geom.pDy2 = 1.280215;
      geom.pDx3 = 1.355978;
      geom.pDx4 = 1.362436;
      geom.centralX = -6.385576;
      geom.centralY = 105.012173;
      geom.centralZ = 20.659979;
      geom.pRotationAngleX = -1.365259;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1093 based Row/Col = 109/3
      geom_tower geom;
      geom.id = 1093;
      geom.pDz = 6.751988;
      geom.pTheta = 0.037977;
      geom.pPhi = 2.828981;
      geom.pAlp1 = -0.007496;
      geom.pAlp2 = -0.007495;
      geom.pDy1 = 1.122326;
      geom.pDx1 = 1.190370;
      geom.pDx2 = 1.195992;
      geom.pDy2 = 1.280215;
      geom.pDx3 = 1.352841;
      geom.pDx4 = 1.359253;
      geom.centralX = -3.828378;
      geom.centralY = 105.012173;
      geom.centralZ = 20.659979;
      geom.pRotationAngleX = -1.365259;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1094 based Row/Col = 109/4
      geom_tower geom;
      geom.id = 1094;
      geom.pDz = 6.751988;
      geom.pTheta = 0.016781;
      geom.pPhi = 2.371419;
      geom.pAlp1 = -0.002496;
      geom.pAlp2 = -0.002495;
      geom.pDy1 = 1.122326;
      geom.pDx1 = 1.188993;
      geom.pDx2 = 1.194595;
      geom.pDy2 = 1.280215;
      geom.pDx3 = 1.351276;
      geom.pDx4 = 1.357665;
      geom.centralX = -1.275632;
      geom.centralY = 105.012173;
      geom.centralZ = 20.659979;
      geom.pRotationAngleX = -1.365259;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1095 based Row/Col = 109/5
      geom_tower geom;
      geom.id = 1095;
      geom.pDz = 6.751988;
      geom.pTheta = 0.016781;
      geom.pPhi = 0.770174;
      geom.pAlp1 = 0.002496;
      geom.pAlp2 = 0.002495;
      geom.pDy1 = 1.122326;
      geom.pDx1 = 1.188993;
      geom.pDx2 = 1.194595;
      geom.pDy2 = 1.280215;
      geom.pDx3 = 1.351276;
      geom.pDx4 = 1.357665;
      geom.centralX = 1.275632;
      geom.centralY = 105.012173;
      geom.centralZ = 20.659979;
      geom.pRotationAngleX = -1.365259;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1096 based Row/Col = 109/6
      geom_tower geom;
      geom.id = 1096;
      geom.pDz = 6.751988;
      geom.pTheta = 0.037977;
      geom.pPhi = 0.312612;
      geom.pAlp1 = 0.007496;
      geom.pAlp2 = 0.007495;
      geom.pDy1 = 1.122326;
      geom.pDx1 = 1.190370;
      geom.pDx2 = 1.195992;
      geom.pDy2 = 1.280215;
      geom.pDx3 = 1.352841;
      geom.pDx4 = 1.359253;
      geom.centralX = 3.828378;
      geom.centralY = 105.012173;
      geom.centralZ = 20.659979;
      geom.pRotationAngleX = -1.365259;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1097 based Row/Col = 109/7
      geom_tower geom;
      geom.id = 1097;
      geom.pDz = 6.751988;
      geom.pTheta = 0.061348;
      geom.pPhi = 0.191404;
      geom.pAlp1 = 0.012523;
      geom.pAlp2 = 0.012521;
      geom.pDy1 = 1.122326;
      geom.pDx1 = 1.193131;
      geom.pDx2 = 1.198793;
      geom.pDy2 = 1.280215;
      geom.pDx3 = 1.355978;
      geom.pDx4 = 1.362436;
      geom.centralX = 6.385576;
      geom.centralY = 105.012173;
      geom.centralZ = 20.659979;
      geom.pRotationAngleX = -1.365259;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1098 based Row/Col = 109/8
      geom_tower geom;
      geom.id = 1098;
      geom.pDz = 6.751988;
      geom.pTheta = 0.085120;
      geom.pPhi = 0.137379;
      geom.pAlp1 = 0.017593;
      geom.pAlp2 = 0.017591;
      geom.pDy1 = 1.122326;
      geom.pDx1 = 1.197289;
      geom.pDx2 = 1.203012;
      geom.pDy2 = 1.280215;
      geom.pDx3 = 1.360702;
      geom.pDx4 = 1.367228;
      geom.centralX = 8.950218;
      geom.centralY = 105.012173;
      geom.centralZ = 20.659979;
      geom.pRotationAngleX = -1.365259;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 911 based Row/Col = 91/1
      geom_tower geom;
      geom.id = 911;
      geom.pDz = 6.751953;
      geom.pTheta = 0.084718;
      geom.pPhi = 3.003391;
      geom.pAlp1 = 0.017594;
      geom.pAlp2 = 0.017591;
      geom.pDy1 = 1.122326;
      geom.pDx1 = 1.197289;
      geom.pDx2 = 1.191567;
      geom.pDy2 = 1.280215;
      geom.pDx3 = 1.360702;
      geom.pDx4 = 1.354176;
      geom.centralX = -8.907948;
      geom.centralY = 104.520810;
      geom.centralZ = -23.016845;
      geom.pRotationAngleX = -1.776334;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 912 based Row/Col = 91/2
      geom_tower geom;
      geom.id = 912;
      geom.pDz = 6.751953;
      geom.pTheta = 0.061064;
      geom.pPhi = 2.949059;
      geom.pAlp1 = 0.012523;
      geom.pAlp2 = 0.012522;
      geom.pDy1 = 1.122326;
      geom.pDx1 = 1.193131;
      geom.pDx2 = 1.187469;
      geom.pDy2 = 1.280215;
      geom.pDx3 = 1.355978;
      geom.pDx4 = 1.349520;
      geom.centralX = -6.355490;
      geom.centralY = 104.520810;
      geom.centralZ = -23.016845;
      geom.pRotationAngleX = -1.776334;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 913 based Row/Col = 91/3
      geom_tower geom;
      geom.id = 913;
      geom.pDz = 6.751953;
      geom.pTheta = 0.037815;
      geom.pPhi = 2.827213;
      geom.pAlp1 = 0.007497;
      geom.pAlp2 = 0.007496;
      geom.pDy1 = 1.122326;
      geom.pDx1 = 1.190370;
      geom.pDx2 = 1.184748;
      geom.pDy2 = 1.280215;
      geom.pDx3 = 1.352841;
      geom.pDx4 = 1.346428;
      geom.centralX = -3.810368;
      geom.centralY = 104.520810;
      geom.centralZ = -23.016845;
      geom.pRotationAngleX = -1.776334;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 914 based Row/Col = 91/4
      geom_tower geom;
      geom.id = 914;
      geom.pDz = 6.751953;
      geom.pTheta = 0.016749;
      geom.pPhi = 2.368410;
      geom.pAlp1 = 0.002496;
      geom.pAlp2 = 0.002496;
      geom.pDy1 = 1.122326;
      geom.pDx1 = 1.188993;
      geom.pDx2 = 1.183390;
      geom.pDy2 = 1.280215;
      geom.pDx3 = 1.351276;
      geom.pDx4 = 1.344885;
      geom.centralX = -1.269636;
      geom.centralY = 104.520810;
      geom.centralZ = -23.016845;
      geom.pRotationAngleX = -1.776334;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 915 based Row/Col = 91/5
      geom_tower geom;
      geom.id = 915;
      geom.pDz = 6.751953;
      geom.pTheta = 0.016749;
      geom.pPhi = 0.773183;
      geom.pAlp1 = -0.002496;
      geom.pAlp2 = -0.002496;
      geom.pDy1 = 1.122326;
      geom.pDx1 = 1.188993;
      geom.pDx2 = 1.183390;
      geom.pDy2 = 1.280215;
      geom.pDx3 = 1.351276;
      geom.pDx4 = 1.344885;
      geom.centralX = 1.269636;
      geom.centralY = 104.520810;
      geom.centralZ = -23.016845;
      geom.pRotationAngleX = -1.776334;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 916 based Row/Col = 91/6
      geom_tower geom;
      geom.id = 916;
      geom.pDz = 6.751953;
      geom.pTheta = 0.037815;
      geom.pPhi = 0.314379;
      geom.pAlp1 = -0.007497;
      geom.pAlp2 = -0.007496;
      geom.pDy1 = 1.122326;
      geom.pDx1 = 1.190370;
      geom.pDx2 = 1.184748;
      geom.pDy2 = 1.280215;
      geom.pDx3 = 1.352841;
      geom.pDx4 = 1.346428;
      geom.centralX = 3.810368;
      geom.centralY = 104.520810;
      geom.centralZ = -23.016845;
      geom.pRotationAngleX = -1.776334;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 917 based Row/Col = 91/7
      geom_tower geom;
      geom.id = 917;
      geom.pDz = 6.751953;
      geom.pTheta = 0.061064;
      geom.pPhi = 0.192533;
      geom.pAlp1 = -0.012523;
      geom.pAlp2 = -0.012522;
      geom.pDy1 = 1.122326;
      geom.pDx1 = 1.193131;
      geom.pDx2 = 1.187469;
      geom.pDy2 = 1.280215;
      geom.pDx3 = 1.355978;
      geom.pDx4 = 1.349520;
      geom.centralX = 6.355490;
      geom.centralY = 104.520810;
      geom.centralZ = -23.016845;
      geom.pRotationAngleX = -1.776334;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 918 based Row/Col = 91/8
      geom_tower geom;
      geom.id = 918;
      geom.pDz = 6.751953;
      geom.pTheta = 0.084718;
      geom.pPhi = 0.138202;
      geom.pAlp1 = -0.017594;
      geom.pAlp2 = -0.017591;
      geom.pDy1 = 1.122326;
      geom.pDx1 = 1.197289;
      geom.pDx2 = 1.191567;
      geom.pDy2 = 1.280215;
      geom.pDx3 = 1.360702;
      geom.pDx4 = 1.354176;
      geom.centralX = 8.907948;
      geom.centralY = 104.520810;
      geom.centralZ = -23.016845;
      geom.pRotationAngleX = -1.776334;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 921 based Row/Col = 92/1
      geom_tower geom;
      geom.id = 921;
      geom.pDz = 6.751988;
      geom.pTheta = 0.085120;
      geom.pPhi = -3.004213;
      geom.pAlp1 = 0.017593;
      geom.pAlp2 = 0.017591;
      geom.pDy1 = 1.122326;
      geom.pDx1 = 1.203012;
      geom.pDx2 = 1.197289;
      geom.pDy2 = 1.280215;
      geom.pDx3 = 1.367228;
      geom.pDx4 = 1.360702;
      geom.centralX = -8.950218;
      geom.centralY = 105.012173;
      geom.centralZ = -20.659979;
      geom.pRotationAngleX = -1.776334;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 922 based Row/Col = 92/2
      geom_tower geom;
      geom.id = 922;
      geom.pDz = 6.751988;
      geom.pTheta = 0.061348;
      geom.pPhi = -2.950189;
      geom.pAlp1 = 0.012523;
      geom.pAlp2 = 0.012521;
      geom.pDy1 = 1.122326;
      geom.pDx1 = 1.198793;
      geom.pDx2 = 1.193131;
      geom.pDy2 = 1.280215;
      geom.pDx3 = 1.362436;
      geom.pDx4 = 1.355978;
      geom.centralX = -6.385576;
      geom.centralY = 105.012173;
      geom.centralZ = -20.659979;
      geom.pRotationAngleX = -1.776334;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 923 based Row/Col = 92/3
      geom_tower geom;
      geom.id = 923;
      geom.pDz = 6.751988;
      geom.pTheta = 0.037977;
      geom.pPhi = -2.828981;
      geom.pAlp1 = 0.007496;
      geom.pAlp2 = 0.007495;
      geom.pDy1 = 1.122326;
      geom.pDx1 = 1.195992;
      geom.pDx2 = 1.190370;
      geom.pDy2 = 1.280215;
      geom.pDx3 = 1.359253;
      geom.pDx4 = 1.352841;
      geom.centralX = -3.828378;
      geom.centralY = 105.012173;
      geom.centralZ = -20.659979;
      geom.pRotationAngleX = -1.776334;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 924 based Row/Col = 92/4
      geom_tower geom;
      geom.id = 924;
      geom.pDz = 6.751988;
      geom.pTheta = 0.016781;
      geom.pPhi = -2.371419;
      geom.pAlp1 = 0.002496;
      geom.pAlp2 = 0.002495;
      geom.pDy1 = 1.122326;
      geom.pDx1 = 1.194595;
      geom.pDx2 = 1.188993;
      geom.pDy2 = 1.280215;
      geom.pDx3 = 1.357665;
      geom.pDx4 = 1.351276;
      geom.centralX = -1.275632;
      geom.centralY = 105.012173;
      geom.centralZ = -20.659979;
      geom.pRotationAngleX = -1.776334;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 925 based Row/Col = 92/5
      geom_tower geom;
      geom.id = 925;
      geom.pDz = 6.751988;
      geom.pTheta = 0.016781;
      geom.pPhi = -0.770174;
      geom.pAlp1 = -0.002496;
      geom.pAlp2 = -0.002495;
      geom.pDy1 = 1.122326;
      geom.pDx1 = 1.194595;
      geom.pDx2 = 1.188993;
      geom.pDy2 = 1.280215;
      geom.pDx3 = 1.357665;
      geom.pDx4 = 1.351276;
      geom.centralX = 1.275632;
      geom.centralY = 105.012173;
      geom.centralZ = -20.659979;
      geom.pRotationAngleX = -1.776334;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 926 based Row/Col = 92/6
      geom_tower geom;
      geom.id = 926;
      geom.pDz = 6.751988;
      geom.pTheta = 0.037977;
      geom.pPhi = -0.312612;
      geom.pAlp1 = -0.007496;
      geom.pAlp2 = -0.007495;
      geom.pDy1 = 1.122326;
      geom.pDx1 = 1.195992;
      geom.pDx2 = 1.190370;
      geom.pDy2 = 1.280215;
      geom.pDx3 = 1.359253;
      geom.pDx4 = 1.352841;
      geom.centralX = 3.828378;
      geom.centralY = 105.012173;
      geom.centralZ = -20.659979;
      geom.pRotationAngleX = -1.776334;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 927 based Row/Col = 92/7
      geom_tower geom;
      geom.id = 927;
      geom.pDz = 6.751988;
      geom.pTheta = 0.061348;
      geom.pPhi = -0.191404;
      geom.pAlp1 = -0.012523;
      geom.pAlp2 = -0.012521;
      geom.pDy1 = 1.122326;
      geom.pDx1 = 1.198793;
      geom.pDx2 = 1.193131;
      geom.pDy2 = 1.280215;
      geom.pDx3 = 1.362436;
      geom.pDx4 = 1.355978;
      geom.centralX = 6.385576;
      geom.centralY = 105.012173;
      geom.centralZ = -20.659979;
      geom.pRotationAngleX = -1.776334;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 928 based Row/Col = 92/8
      geom_tower geom;
      geom.id = 928;
      geom.pDz = 6.751988;
      geom.pTheta = 0.085120;
      geom.pPhi = -0.137379;
      geom.pAlp1 = -0.017593;
      geom.pAlp2 = -0.017591;
      geom.pDy1 = 1.122326;
      geom.pDx1 = 1.203012;
      geom.pDx2 = 1.197289;
      geom.pDy2 = 1.280215;
      geom.pDx3 = 1.367228;
      geom.pDx4 = 1.360702;
      geom.centralX = 8.950218;
      geom.centralY = 105.012173;
      geom.centralZ = -20.659979;
      geom.pRotationAngleX = -1.776334;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1141 based Row/Col = 114/1
      geom_tower geom;
      geom.id = 1141;
      geom.pDz = 6.751882;
      geom.pTheta = 0.082746;
      geom.pPhi = -3.003226;
      geom.pAlp1 = -0.025008;
      geom.pAlp2 = -0.025006;
      geom.pDy1 = 1.123313;
      geom.pDx1 = 1.189896;
      geom.pDx2 = 1.198032;
      geom.pDy2 = 1.277607;
      geom.pDx3 = 1.348497;
      geom.pDx4 = 1.357749;
      geom.centralX = -8.893182;
      geom.centralY = 104.346629;
      geom.centralZ = 32.896588;
      geom.pRotationAngleX = -1.276450;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1142 based Row/Col = 114/2
      geom_tower geom;
      geom.id = 1142;
      geom.pDz = 6.751882;
      geom.pTheta = 0.059644;
      geom.pPhi = -2.948842;
      geom.pAlp1 = -0.017804;
      geom.pAlp2 = -0.017803;
      geom.pDy1 = 1.123313;
      geom.pDx1 = 1.186001;
      geom.pDx2 = 1.194055;
      geom.pDy2 = 1.277607;
      geom.pDx3 = 1.344083;
      geom.pDx4 = 1.353243;
      geom.centralX = -6.345293;
      geom.centralY = 104.346629;
      geom.centralZ = 32.896588;
      geom.pRotationAngleX = -1.276450;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1143 based Row/Col = 114/3
      geom_tower geom;
      geom.id = 1143;
      geom.pDz = 6.751882;
      geom.pTheta = 0.036938;
      geom.pPhi = -2.826883;
      geom.pAlp1 = -0.010659;
      geom.pAlp2 = -0.010659;
      geom.pDy1 = 1.123313;
      geom.pDx1 = 1.183413;
      geom.pDx2 = 1.191414;
      geom.pDy2 = 1.277607;
      geom.pDx3 = 1.341151;
      geom.pDx4 = 1.350251;
      geom.centralX = -3.804390;
      geom.centralY = 104.346629;
      geom.centralZ = 32.896588;
      geom.pRotationAngleX = -1.276450;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1144 based Row/Col = 114/4
      geom_tower geom;
      geom.id = 1144;
      geom.pDz = 6.751882;
      geom.pTheta = 0.016368;
      geom.pPhi = -2.367857;
      geom.pAlp1 = -0.003549;
      geom.pAlp2 = -0.003549;
      geom.pDy1 = 1.123313;
      geom.pDx1 = 1.182122;
      geom.pDx2 = 1.190096;
      geom.pDy2 = 1.277607;
      geom.pDx3 = 1.339689;
      geom.pDx4 = 1.348757;
      geom.centralX = -1.267666;
      geom.centralY = 104.346629;
      geom.centralZ = 32.896588;
      geom.pRotationAngleX = -1.276450;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1145 based Row/Col = 114/5
      geom_tower geom;
      geom.id = 1145;
      geom.pDz = 6.751882;
      geom.pTheta = 0.016368;
      geom.pPhi = -0.773735;
      geom.pAlp1 = 0.003549;
      geom.pAlp2 = 0.003549;
      geom.pDy1 = 1.123313;
      geom.pDx1 = 1.182122;
      geom.pDx2 = 1.190096;
      geom.pDy2 = 1.277607;
      geom.pDx3 = 1.339689;
      geom.pDx4 = 1.348757;
      geom.centralX = 1.267666;
      geom.centralY = 104.346629;
      geom.centralZ = 32.896588;
      geom.pRotationAngleX = -1.276450;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1146 based Row/Col = 114/6
      geom_tower geom;
      geom.id = 1146;
      geom.pDz = 6.751882;
      geom.pTheta = 0.036938;
      geom.pPhi = -0.314710;
      geom.pAlp1 = 0.010659;
      geom.pAlp2 = 0.010659;
      geom.pDy1 = 1.123313;
      geom.pDx1 = 1.183413;
      geom.pDx2 = 1.191414;
      geom.pDy2 = 1.277607;
      geom.pDx3 = 1.341151;
      geom.pDx4 = 1.350251;
      geom.centralX = 3.804390;
      geom.centralY = 104.346629;
      geom.centralZ = 32.896588;
      geom.pRotationAngleX = -1.276450;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1147 based Row/Col = 114/7
      geom_tower geom;
      geom.id = 1147;
      geom.pDz = 6.751882;
      geom.pTheta = 0.059644;
      geom.pPhi = -0.192751;
      geom.pAlp1 = 0.017804;
      geom.pAlp2 = 0.017803;
      geom.pDy1 = 1.123313;
      geom.pDx1 = 1.186001;
      geom.pDx2 = 1.194055;
      geom.pDy2 = 1.277607;
      geom.pDx3 = 1.344083;
      geom.pDx4 = 1.353243;
      geom.centralX = 6.345293;
      geom.centralY = 104.346629;
      geom.centralZ = 32.896588;
      geom.pRotationAngleX = -1.276450;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1148 based Row/Col = 114/8
      geom_tower geom;
      geom.id = 1148;
      geom.pDz = 6.751882;
      geom.pTheta = 0.082746;
      geom.pPhi = -0.138367;
      geom.pAlp1 = 0.025008;
      geom.pAlp2 = 0.025006;
      geom.pDy1 = 1.123313;
      geom.pDx1 = 1.189896;
      geom.pDx2 = 1.198032;
      geom.pDy2 = 1.277607;
      geom.pDx3 = 1.348497;
      geom.pDx4 = 1.357749;
      geom.centralX = 8.893182;
      geom.centralY = 104.346629;
      geom.centralZ = 32.896588;
      geom.pRotationAngleX = -1.276450;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1131 based Row/Col = 113/1
      geom_tower geom;
      geom.id = 1131;
      geom.pDz = 6.751953;
      geom.pTheta = 0.083303;
      geom.pPhi = 3.004475;
      geom.pAlp1 = -0.025007;
      geom.pAlp2 = -0.025005;
      geom.pDy1 = 1.123313;
      geom.pDx1 = 1.198032;
      geom.pDx2 = 1.206169;
      geom.pDy2 = 1.277607;
      geom.pDx3 = 1.357749;
      geom.pDx4 = 1.367003;
      geom.centralX = -8.953233;
      geom.centralY = 105.044621;
      geom.centralZ = 30.594141;
      geom.pRotationAngleX = -1.276450;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1132 based Row/Col = 113/2
      geom_tower geom;
      geom.id = 1132;
      geom.pDz = 6.751953;
      geom.pTheta = 0.060036;
      geom.pPhi = 2.950559;
      geom.pAlp1 = -0.017803;
      geom.pAlp2 = -0.017802;
      geom.pDy1 = 1.123313;
      geom.pDx1 = 1.194055;
      geom.pDx2 = 1.202110;
      geom.pDy2 = 1.277607;
      geom.pDx3 = 1.353243;
      geom.pDx4 = 1.362404;
      geom.centralX = -6.388042;
      geom.centralY = 105.044621;
      geom.centralZ = 30.594141;
      geom.pRotationAngleX = -1.276450;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1133 based Row/Col = 113/3
      geom_tower geom;
      geom.id = 1133;
      geom.pDz = 6.751953;
      geom.pTheta = 0.037161;
      geom.pPhi = 2.829570;
      geom.pAlp1 = -0.010658;
      geom.pAlp2 = -0.010658;
      geom.pDy1 = 1.123313;
      geom.pDx1 = 1.191414;
      geom.pDx2 = 1.199414;
      geom.pDy2 = 1.277607;
      geom.pDx3 = 1.350251;
      geom.pDx4 = 1.359349;
      geom.centralX = -3.829981;
      geom.centralY = 105.044621;
      geom.centralZ = 30.594141;
      geom.pRotationAngleX = -1.276450;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1134 based Row/Col = 113/4
      geom_tower geom;
      geom.id = 1134;
      geom.pDz = 6.751953;
      geom.pTheta = 0.016408;
      geom.pPhi = 2.372434;
      geom.pAlp1 = -0.003549;
      geom.pAlp2 = -0.003549;
      geom.pDy1 = 1.123313;
      geom.pDx1 = 1.190096;
      geom.pDx2 = 1.198069;
      geom.pDy2 = 1.277607;
      geom.pDx3 = 1.348758;
      geom.pDx4 = 1.357825;
      geom.centralX = -1.276187;
      geom.centralY = 105.044621;
      geom.centralZ = 30.594141;
      geom.pRotationAngleX = -1.276450;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1135 based Row/Col = 113/5
      geom_tower geom;
      geom.id = 1135;
      geom.pDz = 6.751953;
      geom.pTheta = 0.016408;
      geom.pPhi = 0.769159;
      geom.pAlp1 = 0.003549;
      geom.pAlp2 = 0.003549;
      geom.pDy1 = 1.123313;
      geom.pDx1 = 1.190096;
      geom.pDx2 = 1.198069;
      geom.pDy2 = 1.277607;
      geom.pDx3 = 1.348758;
      geom.pDx4 = 1.357825;
      geom.centralX = 1.276187;
      geom.centralY = 105.044621;
      geom.centralZ = 30.594141;
      geom.pRotationAngleX = -1.276450;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1136 based Row/Col = 113/6
      geom_tower geom;
      geom.id = 1136;
      geom.pDz = 6.751953;
      geom.pTheta = 0.037161;
      geom.pPhi = 0.312023;
      geom.pAlp1 = 0.010658;
      geom.pAlp2 = 0.010658;
      geom.pDy1 = 1.123313;
      geom.pDx1 = 1.191414;
      geom.pDx2 = 1.199414;
      geom.pDy2 = 1.277607;
      geom.pDx3 = 1.350251;
      geom.pDx4 = 1.359349;
      geom.centralX = 3.829981;
      geom.centralY = 105.044621;
      geom.centralZ = 30.594141;
      geom.pRotationAngleX = -1.276450;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1137 based Row/Col = 113/7
      geom_tower geom;
      geom.id = 1137;
      geom.pDz = 6.751953;
      geom.pTheta = 0.060036;
      geom.pPhi = 0.191034;
      geom.pAlp1 = 0.017803;
      geom.pAlp2 = 0.017802;
      geom.pDy1 = 1.123313;
      geom.pDx1 = 1.194055;
      geom.pDx2 = 1.202110;
      geom.pDy2 = 1.277607;
      geom.pDx3 = 1.353243;
      geom.pDx4 = 1.362404;
      geom.centralX = 6.388042;
      geom.centralY = 105.044621;
      geom.centralZ = 30.594141;
      geom.pRotationAngleX = -1.276450;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1138 based Row/Col = 113/8
      geom_tower geom;
      geom.id = 1138;
      geom.pDz = 6.751953;
      geom.pTheta = 0.083303;
      geom.pPhi = 0.137118;
      geom.pAlp1 = 0.025007;
      geom.pAlp2 = 0.025005;
      geom.pDy1 = 1.123313;
      geom.pDx1 = 1.198032;
      geom.pDx2 = 1.206169;
      geom.pDy2 = 1.277607;
      geom.pDx3 = 1.357749;
      geom.pDx4 = 1.367003;
      geom.centralX = 8.953233;
      geom.centralY = 105.044621;
      geom.centralZ = 30.594141;
      geom.pRotationAngleX = -1.276450;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 871 based Row/Col = 87/1
      geom_tower geom;
      geom.id = 871;
      geom.pDz = 6.751882;
      geom.pTheta = 0.082746;
      geom.pPhi = 3.003226;
      geom.pAlp1 = 0.025008;
      geom.pAlp2 = 0.025006;
      geom.pDy1 = 1.123313;
      geom.pDx1 = 1.198032;
      geom.pDx2 = 1.189896;
      geom.pDy2 = 1.277607;
      geom.pDx3 = 1.357749;
      geom.pDx4 = 1.348497;
      geom.centralX = -8.893182;
      geom.centralY = 104.346629;
      geom.centralZ = -32.896588;
      geom.pRotationAngleX = -1.865143;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 872 based Row/Col = 87/2
      geom_tower geom;
      geom.id = 872;
      geom.pDz = 6.751882;
      geom.pTheta = 0.059644;
      geom.pPhi = 2.948842;
      geom.pAlp1 = 0.017804;
      geom.pAlp2 = 0.017803;
      geom.pDy1 = 1.123313;
      geom.pDx1 = 1.194055;
      geom.pDx2 = 1.186001;
      geom.pDy2 = 1.277607;
      geom.pDx3 = 1.353243;
      geom.pDx4 = 1.344083;
      geom.centralX = -6.345293;
      geom.centralY = 104.346629;
      geom.centralZ = -32.896588;
      geom.pRotationAngleX = -1.865143;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 873 based Row/Col = 87/3
      geom_tower geom;
      geom.id = 873;
      geom.pDz = 6.751882;
      geom.pTheta = 0.036938;
      geom.pPhi = 2.826883;
      geom.pAlp1 = 0.010659;
      geom.pAlp2 = 0.010659;
      geom.pDy1 = 1.123313;
      geom.pDx1 = 1.191414;
      geom.pDx2 = 1.183413;
      geom.pDy2 = 1.277607;
      geom.pDx3 = 1.350251;
      geom.pDx4 = 1.341151;
      geom.centralX = -3.804390;
      geom.centralY = 104.346629;
      geom.centralZ = -32.896588;
      geom.pRotationAngleX = -1.865143;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 874 based Row/Col = 87/4
      geom_tower geom;
      geom.id = 874;
      geom.pDz = 6.751882;
      geom.pTheta = 0.016368;
      geom.pPhi = 2.367857;
      geom.pAlp1 = 0.003549;
      geom.pAlp2 = 0.003549;
      geom.pDy1 = 1.123313;
      geom.pDx1 = 1.190096;
      geom.pDx2 = 1.182122;
      geom.pDy2 = 1.277607;
      geom.pDx3 = 1.348757;
      geom.pDx4 = 1.339689;
      geom.centralX = -1.267666;
      geom.centralY = 104.346629;
      geom.centralZ = -32.896588;
      geom.pRotationAngleX = -1.865143;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 875 based Row/Col = 87/5
      geom_tower geom;
      geom.id = 875;
      geom.pDz = 6.751882;
      geom.pTheta = 0.016368;
      geom.pPhi = 0.773735;
      geom.pAlp1 = -0.003549;
      geom.pAlp2 = -0.003549;
      geom.pDy1 = 1.123313;
      geom.pDx1 = 1.190096;
      geom.pDx2 = 1.182122;
      geom.pDy2 = 1.277607;
      geom.pDx3 = 1.348757;
      geom.pDx4 = 1.339689;
      geom.centralX = 1.267666;
      geom.centralY = 104.346629;
      geom.centralZ = -32.896588;
      geom.pRotationAngleX = -1.865143;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 876 based Row/Col = 87/6
      geom_tower geom;
      geom.id = 876;
      geom.pDz = 6.751882;
      geom.pTheta = 0.036938;
      geom.pPhi = 0.314710;
      geom.pAlp1 = -0.010659;
      geom.pAlp2 = -0.010659;
      geom.pDy1 = 1.123313;
      geom.pDx1 = 1.191414;
      geom.pDx2 = 1.183413;
      geom.pDy2 = 1.277607;
      geom.pDx3 = 1.350251;
      geom.pDx4 = 1.341151;
      geom.centralX = 3.804390;
      geom.centralY = 104.346629;
      geom.centralZ = -32.896588;
      geom.pRotationAngleX = -1.865143;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 877 based Row/Col = 87/7
      geom_tower geom;
      geom.id = 877;
      geom.pDz = 6.751882;
      geom.pTheta = 0.059644;
      geom.pPhi = 0.192751;
      geom.pAlp1 = -0.017804;
      geom.pAlp2 = -0.017803;
      geom.pDy1 = 1.123313;
      geom.pDx1 = 1.194055;
      geom.pDx2 = 1.186001;
      geom.pDy2 = 1.277607;
      geom.pDx3 = 1.353243;
      geom.pDx4 = 1.344083;
      geom.centralX = 6.345293;
      geom.centralY = 104.346629;
      geom.centralZ = -32.896588;
      geom.pRotationAngleX = -1.865143;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 878 based Row/Col = 87/8
      geom_tower geom;
      geom.id = 878;
      geom.pDz = 6.751882;
      geom.pTheta = 0.082746;
      geom.pPhi = 0.138367;
      geom.pAlp1 = -0.025008;
      geom.pAlp2 = -0.025006;
      geom.pDy1 = 1.123313;
      geom.pDx1 = 1.198032;
      geom.pDx2 = 1.189896;
      geom.pDy2 = 1.277607;
      geom.pDx3 = 1.357749;
      geom.pDx4 = 1.348497;
      geom.centralX = 8.893182;
      geom.centralY = 104.346629;
      geom.centralZ = -32.896588;
      geom.pRotationAngleX = -1.865143;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 881 based Row/Col = 88/1
      geom_tower geom;
      geom.id = 881;
      geom.pDz = 6.751953;
      geom.pTheta = 0.083303;
      geom.pPhi = -3.004475;
      geom.pAlp1 = 0.025007;
      geom.pAlp2 = 0.025005;
      geom.pDy1 = 1.123313;
      geom.pDx1 = 1.206169;
      geom.pDx2 = 1.198032;
      geom.pDy2 = 1.277607;
      geom.pDx3 = 1.367003;
      geom.pDx4 = 1.357749;
      geom.centralX = -8.953233;
      geom.centralY = 105.044621;
      geom.centralZ = -30.594141;
      geom.pRotationAngleX = -1.865143;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 882 based Row/Col = 88/2
      geom_tower geom;
      geom.id = 882;
      geom.pDz = 6.751953;
      geom.pTheta = 0.060036;
      geom.pPhi = -2.950559;
      geom.pAlp1 = 0.017803;
      geom.pAlp2 = 0.017802;
      geom.pDy1 = 1.123313;
      geom.pDx1 = 1.202110;
      geom.pDx2 = 1.194055;
      geom.pDy2 = 1.277607;
      geom.pDx3 = 1.362404;
      geom.pDx4 = 1.353243;
      geom.centralX = -6.388042;
      geom.centralY = 105.044621;
      geom.centralZ = -30.594141;
      geom.pRotationAngleX = -1.865143;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 883 based Row/Col = 88/3
      geom_tower geom;
      geom.id = 883;
      geom.pDz = 6.751953;
      geom.pTheta = 0.037161;
      geom.pPhi = -2.829570;
      geom.pAlp1 = 0.010658;
      geom.pAlp2 = 0.010658;
      geom.pDy1 = 1.123313;
      geom.pDx1 = 1.199414;
      geom.pDx2 = 1.191414;
      geom.pDy2 = 1.277607;
      geom.pDx3 = 1.359349;
      geom.pDx4 = 1.350251;
      geom.centralX = -3.829981;
      geom.centralY = 105.044621;
      geom.centralZ = -30.594141;
      geom.pRotationAngleX = -1.865143;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 884 based Row/Col = 88/4
      geom_tower geom;
      geom.id = 884;
      geom.pDz = 6.751953;
      geom.pTheta = 0.016408;
      geom.pPhi = -2.372434;
      geom.pAlp1 = 0.003549;
      geom.pAlp2 = 0.003549;
      geom.pDy1 = 1.123313;
      geom.pDx1 = 1.198069;
      geom.pDx2 = 1.190096;
      geom.pDy2 = 1.277607;
      geom.pDx3 = 1.357825;
      geom.pDx4 = 1.348758;
      geom.centralX = -1.276187;
      geom.centralY = 105.044621;
      geom.centralZ = -30.594141;
      geom.pRotationAngleX = -1.865143;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 885 based Row/Col = 88/5
      geom_tower geom;
      geom.id = 885;
      geom.pDz = 6.751953;
      geom.pTheta = 0.016408;
      geom.pPhi = -0.769159;
      geom.pAlp1 = -0.003549;
      geom.pAlp2 = -0.003549;
      geom.pDy1 = 1.123313;
      geom.pDx1 = 1.198069;
      geom.pDx2 = 1.190096;
      geom.pDy2 = 1.277607;
      geom.pDx3 = 1.357825;
      geom.pDx4 = 1.348758;
      geom.centralX = 1.276187;
      geom.centralY = 105.044621;
      geom.centralZ = -30.594141;
      geom.pRotationAngleX = -1.865143;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 886 based Row/Col = 88/6
      geom_tower geom;
      geom.id = 886;
      geom.pDz = 6.751953;
      geom.pTheta = 0.037161;
      geom.pPhi = -0.312023;
      geom.pAlp1 = -0.010658;
      geom.pAlp2 = -0.010658;
      geom.pDy1 = 1.123313;
      geom.pDx1 = 1.199414;
      geom.pDx2 = 1.191414;
      geom.pDy2 = 1.277607;
      geom.pDx3 = 1.359349;
      geom.pDx4 = 1.350251;
      geom.centralX = 3.829981;
      geom.centralY = 105.044621;
      geom.centralZ = -30.594141;
      geom.pRotationAngleX = -1.865143;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 887 based Row/Col = 88/7
      geom_tower geom;
      geom.id = 887;
      geom.pDz = 6.751953;
      geom.pTheta = 0.060036;
      geom.pPhi = -0.191034;
      geom.pAlp1 = -0.017803;
      geom.pAlp2 = -0.017802;
      geom.pDy1 = 1.123313;
      geom.pDx1 = 1.202110;
      geom.pDx2 = 1.194055;
      geom.pDy2 = 1.277607;
      geom.pDx3 = 1.362404;
      geom.pDx4 = 1.353243;
      geom.centralX = 6.388042;
      geom.centralY = 105.044621;
      geom.centralZ = -30.594141;
      geom.pRotationAngleX = -1.865143;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 888 based Row/Col = 88/8
      geom_tower geom;
      geom.id = 888;
      geom.pDz = 6.751953;
      geom.pTheta = 0.083303;
      geom.pPhi = -0.137118;
      geom.pAlp1 = -0.025007;
      geom.pAlp2 = -0.025005;
      geom.pDy1 = 1.123313;
      geom.pDx1 = 1.206169;
      geom.pDx2 = 1.198032;
      geom.pDy2 = 1.277607;
      geom.pDx3 = 1.367003;
      geom.pDx4 = 1.357749;
      geom.centralX = 8.953233;
      geom.centralY = 105.044621;
      geom.centralZ = -30.594141;
      geom.pRotationAngleX = -1.865143;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1181 based Row/Col = 118/1
      geom_tower geom;
      geom.id = 1181;
      geom.pDz = 6.751871;
      geom.pTheta = 0.080204;
      geom.pPhi = -3.003105;
      geom.pAlp1 = -0.032046;
      geom.pAlp2 = -0.032039;
      geom.pDy1 = 1.124216;
      geom.pDx1 = 1.189569;
      geom.pDx2 = 1.199996;
      geom.pDy2 = 1.273817;
      geom.pDx3 = 1.343082;
      geom.pDx4 = 1.354894;
      geom.centralX = -8.883983;
      geom.centralY = 104.236393;
      geom.centralZ = 43.040967;
      geom.pRotationAngleX = -1.189907;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1182 based Row/Col = 118/2
      geom_tower geom;
      geom.id = 1182;
      geom.pDz = 6.751871;
      geom.pTheta = 0.057813;
      geom.pPhi = -2.948688;
      geom.pAlp1 = -0.022822;
      geom.pAlp2 = -0.022817;
      geom.pDy1 = 1.124216;
      geom.pDx1 = 1.185917;
      geom.pDx2 = 1.196246;
      geom.pDy2 = 1.273817;
      geom.pDx3 = 1.338960;
      geom.pDx4 = 1.350661;
      geom.centralX = -6.339152;
      geom.centralY = 104.236393;
      geom.centralZ = 43.040967;
      geom.pRotationAngleX = -1.189907;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1183 based Row/Col = 118/3
      geom_tower geom;
      geom.id = 1183;
      geom.pDz = 6.751871;
      geom.pTheta = 0.035805;
      geom.pPhi = -2.826656;
      geom.pAlp1 = -0.013666;
      geom.pAlp2 = -0.013663;
      geom.pDy1 = 1.124216;
      geom.pDx1 = 1.183491;
      geom.pDx2 = 1.193755;
      geom.pDy2 = 1.273817;
      geom.pDx3 = 1.336222;
      geom.pDx4 = 1.347849;
      geom.centralX = -3.800876;
      geom.centralY = 104.236393;
      geom.centralZ = 43.040967;
      geom.pRotationAngleX = -1.189907;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1184 based Row/Col = 118/4
      geom_tower geom;
      geom.id = 1184;
      geom.pDz = 6.751871;
      geom.pTheta = 0.015870;
      geom.pPhi = -2.367482;
      geom.pAlp1 = -0.004551;
      geom.pAlp2 = -0.004550;
      geom.pDy1 = 1.124216;
      geom.pDx1 = 1.182280;
      geom.pDx2 = 1.192512;
      geom.pDy2 = 1.273817;
      geom.pDx3 = 1.334855;
      geom.pDx4 = 1.346446;
      geom.centralX = -1.266524;
      geom.centralY = 104.236393;
      geom.centralZ = 43.040967;
      geom.pRotationAngleX = -1.189907;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1185 based Row/Col = 118/5
      geom_tower geom;
      geom.id = 1185;
      geom.pDz = 6.751871;
      geom.pTheta = 0.015870;
      geom.pPhi = -0.774110;
      geom.pAlp1 = 0.004551;
      geom.pAlp2 = 0.004550;
      geom.pDy1 = 1.124216;
      geom.pDx1 = 1.182280;
      geom.pDx2 = 1.192512;
      geom.pDy2 = 1.273817;
      geom.pDx3 = 1.334855;
      geom.pDx4 = 1.346446;
      geom.centralX = 1.266524;
      geom.centralY = 104.236393;
      geom.centralZ = 43.040967;
      geom.pRotationAngleX = -1.189907;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1186 based Row/Col = 118/6
      geom_tower geom;
      geom.id = 1186;
      geom.pDz = 6.751871;
      geom.pTheta = 0.035805;
      geom.pPhi = -0.314937;
      geom.pAlp1 = 0.013666;
      geom.pAlp2 = 0.013663;
      geom.pDy1 = 1.124216;
      geom.pDx1 = 1.183491;
      geom.pDx2 = 1.193755;
      geom.pDy2 = 1.273817;
      geom.pDx3 = 1.336222;
      geom.pDx4 = 1.347849;
      geom.centralX = 3.800876;
      geom.centralY = 104.236393;
      geom.centralZ = 43.040967;
      geom.pRotationAngleX = -1.189907;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1187 based Row/Col = 118/7
      geom_tower geom;
      geom.id = 1187;
      geom.pDz = 6.751871;
      geom.pTheta = 0.057813;
      geom.pPhi = -0.192905;
      geom.pAlp1 = 0.022822;
      geom.pAlp2 = 0.022817;
      geom.pDy1 = 1.124216;
      geom.pDx1 = 1.185917;
      geom.pDx2 = 1.196246;
      geom.pDy2 = 1.273817;
      geom.pDx3 = 1.338960;
      geom.pDx4 = 1.350661;
      geom.centralX = 6.339152;
      geom.centralY = 104.236393;
      geom.centralZ = 43.040967;
      geom.pRotationAngleX = -1.189907;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1188 based Row/Col = 118/8
      geom_tower geom;
      geom.id = 1188;
      geom.pDz = 6.751871;
      geom.pTheta = 0.080204;
      geom.pPhi = -0.138488;
      geom.pAlp1 = 0.032046;
      geom.pAlp2 = 0.032039;
      geom.pDy1 = 1.124216;
      geom.pDx1 = 1.189569;
      geom.pDx2 = 1.199996;
      geom.pDy2 = 1.273817;
      geom.pDx3 = 1.343082;
      geom.pDx4 = 1.354894;
      geom.centralX = 8.883983;
      geom.centralY = 104.236393;
      geom.centralZ = 43.040967;
      geom.pRotationAngleX = -1.189907;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1171 based Row/Col = 117/1
      geom_tower geom;
      geom.id = 1171;
      geom.pDz = 6.751885;
      geom.pTheta = 0.080897;
      geom.pPhi = 3.004732;
      geom.pAlp1 = -0.032045;
      geom.pAlp2 = -0.032038;
      geom.pDy1 = 1.124216;
      geom.pDx1 = 1.199996;
      geom.pDx2 = 1.210425;
      geom.pDy2 = 1.273817;
      geom.pDx3 = 1.354894;
      geom.pDx4 = 1.366708;
      geom.centralX = -8.960847;
      geom.centralY = 105.129711;
      geom.centralZ = 40.810148;
      geom.pRotationAngleX = -1.189907;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1172 based Row/Col = 117/2
      geom_tower geom;
      geom.id = 1172;
      geom.pDz = 6.751885;
      geom.pTheta = 0.058300;
      geom.pPhi = 2.950924;
      geom.pAlp1 = -0.022820;
      geom.pAlp2 = -0.022815;
      geom.pDy1 = 1.124216;
      geom.pDx1 = 1.196246;
      geom.pDx2 = 1.206576;
      geom.pDy2 = 1.273817;
      geom.pDx3 = 1.350661;
      geom.pDx4 = 1.362362;
      geom.centralX = -6.393880;
      geom.centralY = 105.129711;
      geom.centralZ = 40.810148;
      geom.pRotationAngleX = -1.189907;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1173 based Row/Col = 117/3
      geom_tower geom;
      geom.id = 1173;
      geom.pDz = 6.751885;
      geom.pTheta = 0.036082;
      geom.pPhi = 2.830155;
      geom.pAlp1 = -0.013664;
      geom.pAlp2 = -0.013661;
      geom.pDy1 = 1.124216;
      geom.pDx1 = 1.193755;
      geom.pDx2 = 1.204019;
      geom.pDy2 = 1.273817;
      geom.pDx3 = 1.347849;
      geom.pDx4 = 1.359476;
      geom.centralX = -3.833644;
      geom.centralY = 105.129711;
      geom.centralZ = 40.810148;
      geom.pRotationAngleX = -1.189907;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1174 based Row/Col = 117/4
      geom_tower geom;
      geom.id = 1174;
      geom.pDz = 6.751885;
      geom.pTheta = 0.015919;
      geom.pPhi = 2.373446;
      geom.pAlp1 = -0.004550;
      geom.pAlp2 = -0.004549;
      geom.pDy1 = 1.124216;
      geom.pDx1 = 1.192512;
      geom.pDx2 = 1.202743;
      geom.pDy2 = 1.273817;
      geom.pDx3 = 1.346446;
      geom.pDx4 = 1.358036;
      geom.centralX = -1.277434;
      geom.centralY = 105.129711;
      geom.centralZ = 40.810148;
      geom.pRotationAngleX = -1.189907;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1175 based Row/Col = 117/5
      geom_tower geom;
      geom.id = 1175;
      geom.pDz = 6.751885;
      geom.pTheta = 0.015919;
      geom.pPhi = 0.768147;
      geom.pAlp1 = 0.004550;
      geom.pAlp2 = 0.004549;
      geom.pDy1 = 1.124216;
      geom.pDx1 = 1.192512;
      geom.pDx2 = 1.202743;
      geom.pDy2 = 1.273817;
      geom.pDx3 = 1.346446;
      geom.pDx4 = 1.358036;
      geom.centralX = 1.277434;
      geom.centralY = 105.129711;
      geom.centralZ = 40.810148;
      geom.pRotationAngleX = -1.189907;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1176 based Row/Col = 117/6
      geom_tower geom;
      geom.id = 1176;
      geom.pDz = 6.751885;
      geom.pTheta = 0.036082;
      geom.pPhi = 0.311438;
      geom.pAlp1 = 0.013664;
      geom.pAlp2 = 0.013661;
      geom.pDy1 = 1.124216;
      geom.pDx1 = 1.193755;
      geom.pDx2 = 1.204019;
      geom.pDy2 = 1.273817;
      geom.pDx3 = 1.347849;
      geom.pDx4 = 1.359476;
      geom.centralX = 3.833644;
      geom.centralY = 105.129711;
      geom.centralZ = 40.810148;
      geom.pRotationAngleX = -1.189907;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1177 based Row/Col = 117/7
      geom_tower geom;
      geom.id = 1177;
      geom.pDz = 6.751885;
      geom.pTheta = 0.058300;
      geom.pPhi = 0.190669;
      geom.pAlp1 = 0.022820;
      geom.pAlp2 = 0.022815;
      geom.pDy1 = 1.124216;
      geom.pDx1 = 1.196246;
      geom.pDx2 = 1.206576;
      geom.pDy2 = 1.273817;
      geom.pDx3 = 1.350661;
      geom.pDx4 = 1.362362;
      geom.centralX = 6.393880;
      geom.centralY = 105.129711;
      geom.centralZ = 40.810148;
      geom.pRotationAngleX = -1.189907;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1178 based Row/Col = 117/8
      geom_tower geom;
      geom.id = 1178;
      geom.pDz = 6.751885;
      geom.pTheta = 0.080897;
      geom.pPhi = 0.136861;
      geom.pAlp1 = 0.032045;
      geom.pAlp2 = 0.032038;
      geom.pDy1 = 1.124216;
      geom.pDx1 = 1.199996;
      geom.pDx2 = 1.210425;
      geom.pDy2 = 1.273817;
      geom.pDx3 = 1.354894;
      geom.pDx4 = 1.366708;
      geom.centralX = 8.960847;
      geom.centralY = 105.129711;
      geom.centralZ = 40.810148;
      geom.pRotationAngleX = -1.189907;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 831 based Row/Col = 83/1
      geom_tower geom;
      geom.id = 831;
      geom.pDz = 6.751871;
      geom.pTheta = 0.080204;
      geom.pPhi = 3.003105;
      geom.pAlp1 = 0.032046;
      geom.pAlp2 = 0.032039;
      geom.pDy1 = 1.124216;
      geom.pDx1 = 1.199996;
      geom.pDx2 = 1.189569;
      geom.pDy2 = 1.273817;
      geom.pDx3 = 1.354894;
      geom.pDx4 = 1.343082;
      geom.centralX = -8.883983;
      geom.centralY = 104.236393;
      geom.centralZ = -43.040967;
      geom.pRotationAngleX = -1.951685;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 832 based Row/Col = 83/2
      geom_tower geom;
      geom.id = 832;
      geom.pDz = 6.751871;
      geom.pTheta = 0.057813;
      geom.pPhi = 2.948688;
      geom.pAlp1 = 0.022822;
      geom.pAlp2 = 0.022817;
      geom.pDy1 = 1.124216;
      geom.pDx1 = 1.196246;
      geom.pDx2 = 1.185917;
      geom.pDy2 = 1.273817;
      geom.pDx3 = 1.350661;
      geom.pDx4 = 1.338960;
      geom.centralX = -6.339152;
      geom.centralY = 104.236393;
      geom.centralZ = -43.040967;
      geom.pRotationAngleX = -1.951685;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 833 based Row/Col = 83/3
      geom_tower geom;
      geom.id = 833;
      geom.pDz = 6.751871;
      geom.pTheta = 0.035805;
      geom.pPhi = 2.826656;
      geom.pAlp1 = 0.013666;
      geom.pAlp2 = 0.013663;
      geom.pDy1 = 1.124216;
      geom.pDx1 = 1.193755;
      geom.pDx2 = 1.183491;
      geom.pDy2 = 1.273817;
      geom.pDx3 = 1.347849;
      geom.pDx4 = 1.336222;
      geom.centralX = -3.800876;
      geom.centralY = 104.236393;
      geom.centralZ = -43.040967;
      geom.pRotationAngleX = -1.951685;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 834 based Row/Col = 83/4
      geom_tower geom;
      geom.id = 834;
      geom.pDz = 6.751871;
      geom.pTheta = 0.015870;
      geom.pPhi = 2.367482;
      geom.pAlp1 = 0.004551;
      geom.pAlp2 = 0.004550;
      geom.pDy1 = 1.124216;
      geom.pDx1 = 1.192512;
      geom.pDx2 = 1.182280;
      geom.pDy2 = 1.273817;
      geom.pDx3 = 1.346446;
      geom.pDx4 = 1.334855;
      geom.centralX = -1.266524;
      geom.centralY = 104.236393;
      geom.centralZ = -43.040967;
      geom.pRotationAngleX = -1.951685;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 835 based Row/Col = 83/5
      geom_tower geom;
      geom.id = 835;
      geom.pDz = 6.751871;
      geom.pTheta = 0.015870;
      geom.pPhi = 0.774110;
      geom.pAlp1 = -0.004551;
      geom.pAlp2 = -0.004550;
      geom.pDy1 = 1.124216;
      geom.pDx1 = 1.192512;
      geom.pDx2 = 1.182280;
      geom.pDy2 = 1.273817;
      geom.pDx3 = 1.346446;
      geom.pDx4 = 1.334855;
      geom.centralX = 1.266524;
      geom.centralY = 104.236393;
      geom.centralZ = -43.040967;
      geom.pRotationAngleX = -1.951685;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 836 based Row/Col = 83/6
      geom_tower geom;
      geom.id = 836;
      geom.pDz = 6.751871;
      geom.pTheta = 0.035805;
      geom.pPhi = 0.314937;
      geom.pAlp1 = -0.013666;
      geom.pAlp2 = -0.013663;
      geom.pDy1 = 1.124216;
      geom.pDx1 = 1.193755;
      geom.pDx2 = 1.183491;
      geom.pDy2 = 1.273817;
      geom.pDx3 = 1.347849;
      geom.pDx4 = 1.336222;
      geom.centralX = 3.800876;
      geom.centralY = 104.236393;
      geom.centralZ = -43.040967;
      geom.pRotationAngleX = -1.951685;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 837 based Row/Col = 83/7
      geom_tower geom;
      geom.id = 837;
      geom.pDz = 6.751871;
      geom.pTheta = 0.057813;
      geom.pPhi = 0.192905;
      geom.pAlp1 = -0.022822;
      geom.pAlp2 = -0.022817;
      geom.pDy1 = 1.124216;
      geom.pDx1 = 1.196246;
      geom.pDx2 = 1.185917;
      geom.pDy2 = 1.273817;
      geom.pDx3 = 1.350661;
      geom.pDx4 = 1.338960;
      geom.centralX = 6.339152;
      geom.centralY = 104.236393;
      geom.centralZ = -43.040967;
      geom.pRotationAngleX = -1.951685;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 838 based Row/Col = 83/8
      geom_tower geom;
      geom.id = 838;
      geom.pDz = 6.751871;
      geom.pTheta = 0.080204;
      geom.pPhi = 0.138488;
      geom.pAlp1 = -0.032046;
      geom.pAlp2 = -0.032039;
      geom.pDy1 = 1.124216;
      geom.pDx1 = 1.199996;
      geom.pDx2 = 1.189569;
      geom.pDy2 = 1.273817;
      geom.pDx3 = 1.354894;
      geom.pDx4 = 1.343082;
      geom.centralX = 8.883983;
      geom.centralY = 104.236393;
      geom.centralZ = -43.040967;
      geom.pRotationAngleX = -1.951685;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 841 based Row/Col = 84/1
      geom_tower geom;
      geom.id = 841;
      geom.pDz = 6.751885;
      geom.pTheta = 0.080897;
      geom.pPhi = -3.004732;
      geom.pAlp1 = 0.032045;
      geom.pAlp2 = 0.032038;
      geom.pDy1 = 1.124216;
      geom.pDx1 = 1.210425;
      geom.pDx2 = 1.199996;
      geom.pDy2 = 1.273817;
      geom.pDx3 = 1.366708;
      geom.pDx4 = 1.354894;
      geom.centralX = -8.960847;
      geom.centralY = 105.129711;
      geom.centralZ = -40.810148;
      geom.pRotationAngleX = -1.951685;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 842 based Row/Col = 84/2
      geom_tower geom;
      geom.id = 842;
      geom.pDz = 6.751885;
      geom.pTheta = 0.058300;
      geom.pPhi = -2.950924;
      geom.pAlp1 = 0.022820;
      geom.pAlp2 = 0.022815;
      geom.pDy1 = 1.124216;
      geom.pDx1 = 1.206576;
      geom.pDx2 = 1.196246;
      geom.pDy2 = 1.273817;
      geom.pDx3 = 1.362362;
      geom.pDx4 = 1.350661;
      geom.centralX = -6.393880;
      geom.centralY = 105.129711;
      geom.centralZ = -40.810148;
      geom.pRotationAngleX = -1.951685;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 843 based Row/Col = 84/3
      geom_tower geom;
      geom.id = 843;
      geom.pDz = 6.751885;
      geom.pTheta = 0.036082;
      geom.pPhi = -2.830155;
      geom.pAlp1 = 0.013664;
      geom.pAlp2 = 0.013661;
      geom.pDy1 = 1.124216;
      geom.pDx1 = 1.204019;
      geom.pDx2 = 1.193755;
      geom.pDy2 = 1.273817;
      geom.pDx3 = 1.359476;
      geom.pDx4 = 1.347849;
      geom.centralX = -3.833644;
      geom.centralY = 105.129711;
      geom.centralZ = -40.810148;
      geom.pRotationAngleX = -1.951685;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 844 based Row/Col = 84/4
      geom_tower geom;
      geom.id = 844;
      geom.pDz = 6.751885;
      geom.pTheta = 0.015919;
      geom.pPhi = -2.373446;
      geom.pAlp1 = 0.004550;
      geom.pAlp2 = 0.004549;
      geom.pDy1 = 1.124216;
      geom.pDx1 = 1.202743;
      geom.pDx2 = 1.192512;
      geom.pDy2 = 1.273817;
      geom.pDx3 = 1.358036;
      geom.pDx4 = 1.346446;
      geom.centralX = -1.277434;
      geom.centralY = 105.129711;
      geom.centralZ = -40.810148;
      geom.pRotationAngleX = -1.951685;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 845 based Row/Col = 84/5
      geom_tower geom;
      geom.id = 845;
      geom.pDz = 6.751885;
      geom.pTheta = 0.015919;
      geom.pPhi = -0.768147;
      geom.pAlp1 = -0.004550;
      geom.pAlp2 = -0.004549;
      geom.pDy1 = 1.124216;
      geom.pDx1 = 1.202743;
      geom.pDx2 = 1.192512;
      geom.pDy2 = 1.273817;
      geom.pDx3 = 1.358036;
      geom.pDx4 = 1.346446;
      geom.centralX = 1.277434;
      geom.centralY = 105.129711;
      geom.centralZ = -40.810148;
      geom.pRotationAngleX = -1.951685;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 846 based Row/Col = 84/6
      geom_tower geom;
      geom.id = 846;
      geom.pDz = 6.751885;
      geom.pTheta = 0.036082;
      geom.pPhi = -0.311438;
      geom.pAlp1 = -0.013664;
      geom.pAlp2 = -0.013661;
      geom.pDy1 = 1.124216;
      geom.pDx1 = 1.204019;
      geom.pDx2 = 1.193755;
      geom.pDy2 = 1.273817;
      geom.pDx3 = 1.359476;
      geom.pDx4 = 1.347849;
      geom.centralX = 3.833644;
      geom.centralY = 105.129711;
      geom.centralZ = -40.810148;
      geom.pRotationAngleX = -1.951685;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 847 based Row/Col = 84/7
      geom_tower geom;
      geom.id = 847;
      geom.pDz = 6.751885;
      geom.pTheta = 0.058300;
      geom.pPhi = -0.190669;
      geom.pAlp1 = -0.022820;
      geom.pAlp2 = -0.022815;
      geom.pDy1 = 1.124216;
      geom.pDx1 = 1.206576;
      geom.pDx2 = 1.196246;
      geom.pDy2 = 1.273817;
      geom.pDx3 = 1.362362;
      geom.pDx4 = 1.350661;
      geom.centralX = 6.393880;
      geom.centralY = 105.129711;
      geom.centralZ = -40.810148;
      geom.pRotationAngleX = -1.951685;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 848 based Row/Col = 84/8
      geom_tower geom;
      geom.id = 848;
      geom.pDz = 6.751885;
      geom.pTheta = 0.080897;
      geom.pPhi = -0.136861;
      geom.pAlp1 = -0.032045;
      geom.pAlp2 = -0.032038;
      geom.pDy1 = 1.124216;
      geom.pDx1 = 1.210425;
      geom.pDx2 = 1.199996;
      geom.pDy2 = 1.273817;
      geom.pDx3 = 1.366708;
      geom.pDx4 = 1.354894;
      geom.centralX = 8.960847;
      geom.centralY = 105.129711;
      geom.centralZ = -40.810148;
      geom.pRotationAngleX = -1.951685;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1221 based Row/Col = 122/1
      geom_tower geom;
      geom.id = 1221;
      geom.pDz = 6.751869;
      geom.pTheta = 0.077209;
      geom.pPhi = -3.001570;
      geom.pAlp1 = -0.038599;
      geom.pAlp2 = -0.038586;
      geom.pDy1 = 1.125755;
      geom.pDx1 = 1.190458;
      geom.pDx2 = 1.203026;
      geom.pDy2 = 1.270671;
      geom.pDx3 = 1.337997;
      geom.pDx4 = 1.352177;
      geom.centralX = -8.879916;
      geom.centralY = 104.185244;
      geom.centralZ = 53.539930;
      geom.pRotationAngleX = -1.106546;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1222 based Row/Col = 122/2
      geom_tower geom;
      geom.id = 1222;
      geom.pDz = 6.751869;
      geom.pTheta = 0.055664;
      geom.pPhi = -2.946590;
      geom.pAlp1 = -0.027497;
      geom.pAlp2 = -0.027487;
      geom.pDy1 = 1.125755;
      geom.pDx1 = 1.187080;
      geom.pDx2 = 1.199537;
      geom.pDy2 = 1.270671;
      geom.pDx3 = 1.334201;
      geom.pDx4 = 1.348257;
      geom.centralX = -6.336733;
      geom.centralY = 104.185244;
      geom.centralZ = 53.539930;
      geom.pRotationAngleX = -1.106546;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1223 based Row/Col = 122/3
      geom_tower geom;
      geom.id = 1223;
      geom.pDz = 6.751869;
      geom.pTheta = 0.034497;
      geom.pPhi = -2.823389;
      geom.pAlp1 = -0.016468;
      geom.pAlp2 = -0.016463;
      geom.pDy1 = 1.125755;
      geom.pDx1 = 1.184834;
      geom.pDx2 = 1.197219;
      geom.pDy2 = 1.270671;
      geom.pDx3 = 1.331678;
      geom.pDx4 = 1.345652;
      geom.centralX = -3.799619;
      geom.centralY = 104.185244;
      geom.centralZ = 53.539930;
      geom.pRotationAngleX = -1.106546;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1224 based Row/Col = 122/4
      geom_tower geom;
      geom.id = 1224;
      geom.pDz = 6.751869;
      geom.pTheta = 0.015357;
      geom.pPhi = -2.361975;
      geom.pAlp1 = -0.005485;
      geom.pAlp2 = -0.005483;
      geom.pDy1 = 1.125755;
      geom.pDx1 = 1.183714;
      geom.pDx2 = 1.196062;
      geom.pDy2 = 1.270671;
      geom.pDx3 = 1.330418;
      geom.pDx4 = 1.344352;
      geom.centralX = -1.266137;
      geom.centralY = 104.185244;
      geom.centralZ = 53.539930;
      geom.pRotationAngleX = -1.106546;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1225 based Row/Col = 122/5
      geom_tower geom;
      geom.id = 1225;
      geom.pDz = 6.751869;
      geom.pTheta = 0.015357;
      geom.pPhi = -0.779618;
      geom.pAlp1 = 0.005485;
      geom.pAlp2 = 0.005483;
      geom.pDy1 = 1.125755;
      geom.pDx1 = 1.183714;
      geom.pDx2 = 1.196062;
      geom.pDy2 = 1.270671;
      geom.pDx3 = 1.330418;
      geom.pDx4 = 1.344352;
      geom.centralX = 1.266137;
      geom.centralY = 104.185244;
      geom.centralZ = 53.539930;
      geom.pRotationAngleX = -1.106546;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1226 based Row/Col = 122/6
      geom_tower geom;
      geom.id = 1226;
      geom.pDz = 6.751869;
      geom.pTheta = 0.034497;
      geom.pPhi = -0.318204;
      geom.pAlp1 = 0.016468;
      geom.pAlp2 = 0.016463;
      geom.pDy1 = 1.125755;
      geom.pDx1 = 1.184834;
      geom.pDx2 = 1.197219;
      geom.pDy2 = 1.270671;
      geom.pDx3 = 1.331678;
      geom.pDx4 = 1.345652;
      geom.centralX = 3.799619;
      geom.centralY = 104.185244;
      geom.centralZ = 53.539930;
      geom.pRotationAngleX = -1.106546;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1227 based Row/Col = 122/7
      geom_tower geom;
      geom.id = 1227;
      geom.pDz = 6.751869;
      geom.pTheta = 0.055664;
      geom.pPhi = -0.195002;
      geom.pAlp1 = 0.027497;
      geom.pAlp2 = 0.027487;
      geom.pDy1 = 1.125755;
      geom.pDx1 = 1.187080;
      geom.pDx2 = 1.199537;
      geom.pDy2 = 1.270671;
      geom.pDx3 = 1.334201;
      geom.pDx4 = 1.348257;
      geom.centralX = 6.336733;
      geom.centralY = 104.185244;
      geom.centralZ = 53.539930;
      geom.pRotationAngleX = -1.106546;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1228 based Row/Col = 122/8
      geom_tower geom;
      geom.id = 1228;
      geom.pDz = 6.751869;
      geom.pTheta = 0.077209;
      geom.pPhi = -0.140023;
      geom.pAlp1 = 0.038599;
      geom.pAlp2 = 0.038586;
      geom.pDy1 = 1.125755;
      geom.pDx1 = 1.190458;
      geom.pDx2 = 1.203026;
      geom.pDy2 = 1.270671;
      geom.pDx3 = 1.337997;
      geom.pDx4 = 1.352177;
      geom.centralX = 8.879916;
      geom.centralY = 104.185244;
      geom.centralZ = 53.539930;
      geom.pRotationAngleX = -1.106546;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1211 based Row/Col = 121/1
      geom_tower geom;
      geom.id = 1211;
      geom.pDz = 6.751810;
      geom.pTheta = 0.078005;
      geom.pPhi = 3.004708;
      geom.pAlp1 = -0.038598;
      geom.pAlp2 = -0.038585;
      geom.pDy1 = 1.125755;
      geom.pDx1 = 1.203026;
      geom.pDx2 = 1.215596;
      geom.pDy2 = 1.270671;
      geom.pDx3 = 1.352177;
      geom.pDx4 = 1.366360;
      geom.centralX = -8.972444;
      geom.centralY = 105.260489;
      geom.centralZ = 51.392677;
      geom.pRotationAngleX = -1.106546;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1212 based Row/Col = 121/2
      geom_tower geom;
      geom.id = 1212;
      geom.pDz = 6.751810;
      geom.pTheta = 0.056215;
      geom.pPhi = 2.950906;
      geom.pAlp1 = -0.027494;
      geom.pAlp2 = -0.027485;
      geom.pDy1 = 1.125755;
      geom.pDx1 = 1.199537;
      geom.pDx2 = 1.211995;
      geom.pDy2 = 1.270671;
      geom.pDx3 = 1.348257;
      geom.pDx4 = 1.362314;
      geom.centralX = -6.402628;
      geom.centralY = 105.260489;
      geom.centralZ = 51.392677;
      geom.pRotationAngleX = -1.106546;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1213 based Row/Col = 121/3
      geom_tower geom;
      geom.id = 1213;
      geom.pDz = 6.751810;
      geom.pTheta = 0.034792;
      geom.pPhi = 2.830140;
      geom.pAlp1 = -0.016466;
      geom.pAlp2 = -0.016461;
      geom.pDy1 = 1.125755;
      geom.pDx1 = 1.197219;
      geom.pDx2 = 1.209603;
      geom.pDy2 = 1.270671;
      geom.pDx3 = 1.345652;
      geom.pDx4 = 1.359625;
      geom.centralX = -3.839078;
      geom.centralY = 105.260489;
      geom.centralZ = 51.392677;
      geom.pRotationAngleX = -1.106546;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1214 based Row/Col = 121/4
      geom_tower geom;
      geom.id = 1214;
      geom.pDz = 6.751810;
      geom.pTheta = 0.015350;
      geom.pPhi = 2.373433;
      geom.pAlp1 = -0.005484;
      geom.pAlp2 = -0.005482;
      geom.pDy1 = 1.125755;
      geom.pDx1 = 1.196062;
      geom.pDx2 = 1.208409;
      geom.pDy2 = 1.270671;
      geom.pDx3 = 1.344352;
      geom.pDx4 = 1.358283;
      geom.centralX = -1.279277;
      geom.centralY = 105.260489;
      geom.centralZ = 51.392677;
      geom.pRotationAngleX = -1.106546;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1215 based Row/Col = 121/5
      geom_tower geom;
      geom.id = 1215;
      geom.pDz = 6.751810;
      geom.pTheta = 0.015350;
      geom.pPhi = 0.768160;
      geom.pAlp1 = 0.005484;
      geom.pAlp2 = 0.005482;
      geom.pDy1 = 1.125755;
      geom.pDx1 = 1.196062;
      geom.pDx2 = 1.208409;
      geom.pDy2 = 1.270671;
      geom.pDx3 = 1.344352;
      geom.pDx4 = 1.358283;
      geom.centralX = 1.279277;
      geom.centralY = 105.260489;
      geom.centralZ = 51.392677;
      geom.pRotationAngleX = -1.106546;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1216 based Row/Col = 121/6
      geom_tower geom;
      geom.id = 1216;
      geom.pDz = 6.751810;
      geom.pTheta = 0.034792;
      geom.pPhi = 0.311453;
      geom.pAlp1 = 0.016466;
      geom.pAlp2 = 0.016461;
      geom.pDy1 = 1.125755;
      geom.pDx1 = 1.197219;
      geom.pDx2 = 1.209603;
      geom.pDy2 = 1.270671;
      geom.pDx3 = 1.345652;
      geom.pDx4 = 1.359625;
      geom.centralX = 3.839078;
      geom.centralY = 105.260489;
      geom.centralZ = 51.392677;
      geom.pRotationAngleX = -1.106546;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1217 based Row/Col = 121/7
      geom_tower geom;
      geom.id = 1217;
      geom.pDz = 6.751810;
      geom.pTheta = 0.056215;
      geom.pPhi = 0.190687;
      geom.pAlp1 = 0.027494;
      geom.pAlp2 = 0.027485;
      geom.pDy1 = 1.125755;
      geom.pDx1 = 1.199537;
      geom.pDx2 = 1.211995;
      geom.pDy2 = 1.270671;
      geom.pDx3 = 1.348257;
      geom.pDx4 = 1.362314;
      geom.centralX = 6.402628;
      geom.centralY = 105.260489;
      geom.centralZ = 51.392677;
      geom.pRotationAngleX = -1.106546;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1218 based Row/Col = 121/8
      geom_tower geom;
      geom.id = 1218;
      geom.pDz = 6.751810;
      geom.pTheta = 0.078005;
      geom.pPhi = 0.136884;
      geom.pAlp1 = 0.038598;
      geom.pAlp2 = 0.038585;
      geom.pDy1 = 1.125755;
      geom.pDx1 = 1.203026;
      geom.pDx2 = 1.215596;
      geom.pDy2 = 1.270671;
      geom.pDx3 = 1.352177;
      geom.pDx4 = 1.366360;
      geom.centralX = 8.972444;
      geom.centralY = 105.260489;
      geom.centralZ = 51.392677;
      geom.pRotationAngleX = -1.106546;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 791 based Row/Col = 79/1
      geom_tower geom;
      geom.id = 791;
      geom.pDz = 6.751869;
      geom.pTheta = 0.077209;
      geom.pPhi = 3.001570;
      geom.pAlp1 = 0.038599;
      geom.pAlp2 = 0.038586;
      geom.pDy1 = 1.125755;
      geom.pDx1 = 1.203026;
      geom.pDx2 = 1.190458;
      geom.pDy2 = 1.270671;
      geom.pDx3 = 1.352177;
      geom.pDx4 = 1.337997;
      geom.centralX = -8.879916;
      geom.centralY = 104.185244;
      geom.centralZ = -53.539930;
      geom.pRotationAngleX = -2.035047;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 792 based Row/Col = 79/2
      geom_tower geom;
      geom.id = 792;
      geom.pDz = 6.751869;
      geom.pTheta = 0.055664;
      geom.pPhi = 2.946590;
      geom.pAlp1 = 0.027497;
      geom.pAlp2 = 0.027487;
      geom.pDy1 = 1.125755;
      geom.pDx1 = 1.199537;
      geom.pDx2 = 1.187080;
      geom.pDy2 = 1.270671;
      geom.pDx3 = 1.348257;
      geom.pDx4 = 1.334201;
      geom.centralX = -6.336733;
      geom.centralY = 104.185244;
      geom.centralZ = -53.539930;
      geom.pRotationAngleX = -2.035047;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 793 based Row/Col = 79/3
      geom_tower geom;
      geom.id = 793;
      geom.pDz = 6.751869;
      geom.pTheta = 0.034497;
      geom.pPhi = 2.823389;
      geom.pAlp1 = 0.016468;
      geom.pAlp2 = 0.016463;
      geom.pDy1 = 1.125755;
      geom.pDx1 = 1.197219;
      geom.pDx2 = 1.184834;
      geom.pDy2 = 1.270671;
      geom.pDx3 = 1.345652;
      geom.pDx4 = 1.331678;
      geom.centralX = -3.799619;
      geom.centralY = 104.185244;
      geom.centralZ = -53.539930;
      geom.pRotationAngleX = -2.035047;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 794 based Row/Col = 79/4
      geom_tower geom;
      geom.id = 794;
      geom.pDz = 6.751869;
      geom.pTheta = 0.015357;
      geom.pPhi = 2.361975;
      geom.pAlp1 = 0.005485;
      geom.pAlp2 = 0.005483;
      geom.pDy1 = 1.125755;
      geom.pDx1 = 1.196062;
      geom.pDx2 = 1.183714;
      geom.pDy2 = 1.270671;
      geom.pDx3 = 1.344352;
      geom.pDx4 = 1.330418;
      geom.centralX = -1.266137;
      geom.centralY = 104.185244;
      geom.centralZ = -53.539930;
      geom.pRotationAngleX = -2.035047;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 795 based Row/Col = 79/5
      geom_tower geom;
      geom.id = 795;
      geom.pDz = 6.751869;
      geom.pTheta = 0.015357;
      geom.pPhi = 0.779618;
      geom.pAlp1 = -0.005485;
      geom.pAlp2 = -0.005483;
      geom.pDy1 = 1.125755;
      geom.pDx1 = 1.196062;
      geom.pDx2 = 1.183714;
      geom.pDy2 = 1.270671;
      geom.pDx3 = 1.344352;
      geom.pDx4 = 1.330418;
      geom.centralX = 1.266137;
      geom.centralY = 104.185244;
      geom.centralZ = -53.539930;
      geom.pRotationAngleX = -2.035047;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 796 based Row/Col = 79/6
      geom_tower geom;
      geom.id = 796;
      geom.pDz = 6.751869;
      geom.pTheta = 0.034497;
      geom.pPhi = 0.318204;
      geom.pAlp1 = -0.016468;
      geom.pAlp2 = -0.016463;
      geom.pDy1 = 1.125755;
      geom.pDx1 = 1.197219;
      geom.pDx2 = 1.184834;
      geom.pDy2 = 1.270671;
      geom.pDx3 = 1.345652;
      geom.pDx4 = 1.331678;
      geom.centralX = 3.799619;
      geom.centralY = 104.185244;
      geom.centralZ = -53.539930;
      geom.pRotationAngleX = -2.035047;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 797 based Row/Col = 79/7
      geom_tower geom;
      geom.id = 797;
      geom.pDz = 6.751869;
      geom.pTheta = 0.055664;
      geom.pPhi = 0.195002;
      geom.pAlp1 = -0.027497;
      geom.pAlp2 = -0.027487;
      geom.pDy1 = 1.125755;
      geom.pDx1 = 1.199537;
      geom.pDx2 = 1.187080;
      geom.pDy2 = 1.270671;
      geom.pDx3 = 1.348257;
      geom.pDx4 = 1.334201;
      geom.centralX = 6.336733;
      geom.centralY = 104.185244;
      geom.centralZ = -53.539930;
      geom.pRotationAngleX = -2.035047;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 798 based Row/Col = 79/8
      geom_tower geom;
      geom.id = 798;
      geom.pDz = 6.751869;
      geom.pTheta = 0.077209;
      geom.pPhi = 0.140023;
      geom.pAlp1 = -0.038599;
      geom.pAlp2 = -0.038586;
      geom.pDy1 = 1.125755;
      geom.pDx1 = 1.203026;
      geom.pDx2 = 1.190458;
      geom.pDy2 = 1.270671;
      geom.pDx3 = 1.352177;
      geom.pDx4 = 1.337997;
      geom.centralX = 8.879916;
      geom.centralY = 104.185244;
      geom.centralZ = -53.539930;
      geom.pRotationAngleX = -2.035047;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 801 based Row/Col = 80/1
      geom_tower geom;
      geom.id = 801;
      geom.pDz = 6.751810;
      geom.pTheta = 0.078005;
      geom.pPhi = -3.004708;
      geom.pAlp1 = 0.038598;
      geom.pAlp2 = 0.038585;
      geom.pDy1 = 1.125755;
      geom.pDx1 = 1.215596;
      geom.pDx2 = 1.203026;
      geom.pDy2 = 1.270671;
      geom.pDx3 = 1.366360;
      geom.pDx4 = 1.352177;
      geom.centralX = -8.972444;
      geom.centralY = 105.260489;
      geom.centralZ = -51.392677;
      geom.pRotationAngleX = -2.035047;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 802 based Row/Col = 80/2
      geom_tower geom;
      geom.id = 802;
      geom.pDz = 6.751810;
      geom.pTheta = 0.056215;
      geom.pPhi = -2.950906;
      geom.pAlp1 = 0.027494;
      geom.pAlp2 = 0.027485;
      geom.pDy1 = 1.125755;
      geom.pDx1 = 1.211995;
      geom.pDx2 = 1.199537;
      geom.pDy2 = 1.270671;
      geom.pDx3 = 1.362314;
      geom.pDx4 = 1.348257;
      geom.centralX = -6.402628;
      geom.centralY = 105.260489;
      geom.centralZ = -51.392677;
      geom.pRotationAngleX = -2.035047;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 803 based Row/Col = 80/3
      geom_tower geom;
      geom.id = 803;
      geom.pDz = 6.751810;
      geom.pTheta = 0.034792;
      geom.pPhi = -2.830140;
      geom.pAlp1 = 0.016466;
      geom.pAlp2 = 0.016461;
      geom.pDy1 = 1.125755;
      geom.pDx1 = 1.209603;
      geom.pDx2 = 1.197219;
      geom.pDy2 = 1.270671;
      geom.pDx3 = 1.359625;
      geom.pDx4 = 1.345652;
      geom.centralX = -3.839078;
      geom.centralY = 105.260489;
      geom.centralZ = -51.392677;
      geom.pRotationAngleX = -2.035047;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 804 based Row/Col = 80/4
      geom_tower geom;
      geom.id = 804;
      geom.pDz = 6.751810;
      geom.pTheta = 0.015350;
      geom.pPhi = -2.373433;
      geom.pAlp1 = 0.005484;
      geom.pAlp2 = 0.005482;
      geom.pDy1 = 1.125755;
      geom.pDx1 = 1.208409;
      geom.pDx2 = 1.196062;
      geom.pDy2 = 1.270671;
      geom.pDx3 = 1.358283;
      geom.pDx4 = 1.344352;
      geom.centralX = -1.279277;
      geom.centralY = 105.260489;
      geom.centralZ = -51.392677;
      geom.pRotationAngleX = -2.035047;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 805 based Row/Col = 80/5
      geom_tower geom;
      geom.id = 805;
      geom.pDz = 6.751810;
      geom.pTheta = 0.015350;
      geom.pPhi = -0.768160;
      geom.pAlp1 = -0.005484;
      geom.pAlp2 = -0.005482;
      geom.pDy1 = 1.125755;
      geom.pDx1 = 1.208409;
      geom.pDx2 = 1.196062;
      geom.pDy2 = 1.270671;
      geom.pDx3 = 1.358283;
      geom.pDx4 = 1.344352;
      geom.centralX = 1.279277;
      geom.centralY = 105.260489;
      geom.centralZ = -51.392677;
      geom.pRotationAngleX = -2.035047;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 806 based Row/Col = 80/6
      geom_tower geom;
      geom.id = 806;
      geom.pDz = 6.751810;
      geom.pTheta = 0.034792;
      geom.pPhi = -0.311453;
      geom.pAlp1 = -0.016466;
      geom.pAlp2 = -0.016461;
      geom.pDy1 = 1.125755;
      geom.pDx1 = 1.209603;
      geom.pDx2 = 1.197219;
      geom.pDy2 = 1.270671;
      geom.pDx3 = 1.359625;
      geom.pDx4 = 1.345652;
      geom.centralX = 3.839078;
      geom.centralY = 105.260489;
      geom.centralZ = -51.392677;
      geom.pRotationAngleX = -2.035047;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 807 based Row/Col = 80/7
      geom_tower geom;
      geom.id = 807;
      geom.pDz = 6.751810;
      geom.pTheta = 0.056215;
      geom.pPhi = -0.190687;
      geom.pAlp1 = -0.027494;
      geom.pAlp2 = -0.027485;
      geom.pDy1 = 1.125755;
      geom.pDx1 = 1.211995;
      geom.pDx2 = 1.199537;
      geom.pDy2 = 1.270671;
      geom.pDx3 = 1.362314;
      geom.pDx4 = 1.348257;
      geom.centralX = 6.402628;
      geom.centralY = 105.260489;
      geom.centralZ = -51.392677;
      geom.pRotationAngleX = -2.035047;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 808 based Row/Col = 80/8
      geom_tower geom;
      geom.id = 808;
      geom.pDz = 6.751810;
      geom.pTheta = 0.078005;
      geom.pPhi = -0.136884;
      geom.pAlp1 = -0.038598;
      geom.pAlp2 = -0.038585;
      geom.pDy1 = 1.125755;
      geom.pDx1 = 1.215596;
      geom.pDx2 = 1.203026;
      geom.pDy2 = 1.270671;
      geom.pDx3 = 1.366360;
      geom.pDx4 = 1.352177;
      geom.centralX = 8.972444;
      geom.centralY = 105.260489;
      geom.centralZ = -51.392677;
      geom.pRotationAngleX = -2.035047;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1261 based Row/Col = 126/1
      geom_tower geom;
      geom.id = 1261;
      geom.pDz = 6.751789;
      geom.pTheta = 0.073828;
      geom.pPhi = -3.000754;
      geom.pAlp1 = -0.044637;
      geom.pAlp2 = -0.044625;
      geom.pDy1 = 1.127101;
      geom.pDx1 = 1.192450;
      geom.pDx2 = 1.206989;
      geom.pDy2 = 1.266192;
      geom.pDx3 = 1.333316;
      geom.pDx4 = 1.349644;
      geom.centralX = -8.880634;
      geom.centralY = 104.189340;
      geom.centralZ = 64.494347;
      geom.pRotationAngleX = -1.026472;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1262 based Row/Col = 126/2
      geom_tower geom;
      geom.id = 1262;
      geom.pDz = 6.751789;
      geom.pTheta = 0.053232;
      geom.pPhi = -2.945484;
      geom.pAlp1 = -0.031808;
      geom.pAlp2 = -0.031800;
      geom.pDy1 = 1.127101;
      geom.pDx1 = 1.189362;
      geom.pDx2 = 1.203784;
      geom.pDy2 = 1.266192;
      geom.pDx3 = 1.329863;
      geom.pDx4 = 1.346061;
      geom.centralX = -6.337767;
      geom.centralY = 104.189340;
      geom.centralZ = 64.494347;
      geom.pRotationAngleX = -1.026472;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1263 based Row/Col = 126/3
      geom_tower geom;
      geom.id = 1263;
      geom.pDz = 6.751789;
      geom.pTheta = 0.033002;
      geom.pPhi = -2.821676;
      geom.pAlp1 = -0.019055;
      geom.pAlp2 = -0.019050;
      geom.pDy1 = 1.127101;
      geom.pDx1 = 1.187309;
      geom.pDx2 = 1.201654;
      geom.pDy2 = 1.266192;
      geom.pDx3 = 1.327568;
      geom.pDx4 = 1.343680;
      geom.centralX = -3.800447;
      geom.centralY = 104.189340;
      geom.centralZ = 64.494347;
      geom.pRotationAngleX = -1.026472;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1264 based Row/Col = 126/4
      geom_tower geom;
      geom.id = 1264;
      geom.pDz = 6.751789;
      geom.pTheta = 0.014725;
      geom.pPhi = -2.359114;
      geom.pAlp1 = -0.006347;
      geom.pAlp2 = -0.006345;
      geom.pDy1 = 1.127101;
      geom.pDx1 = 1.186284;
      geom.pDx2 = 1.200591;
      geom.pDy2 = 1.266192;
      geom.pDx3 = 1.326423;
      geom.pDx4 = 1.342491;
      geom.centralX = -1.266447;
      geom.centralY = 104.189340;
      geom.centralZ = 64.494347;
      geom.pRotationAngleX = -1.026472;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1265 based Row/Col = 126/5
      geom_tower geom;
      geom.id = 1265;
      geom.pDz = 6.751789;
      geom.pTheta = 0.014725;
      geom.pPhi = -0.782479;
      geom.pAlp1 = 0.006347;
      geom.pAlp2 = 0.006345;
      geom.pDy1 = 1.127101;
      geom.pDx1 = 1.186284;
      geom.pDx2 = 1.200591;
      geom.pDy2 = 1.266192;
      geom.pDx3 = 1.326423;
      geom.pDx4 = 1.342491;
      geom.centralX = 1.266447;
      geom.centralY = 104.189340;
      geom.centralZ = 64.494347;
      geom.pRotationAngleX = -1.026472;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1266 based Row/Col = 126/6
      geom_tower geom;
      geom.id = 1266;
      geom.pDz = 6.751789;
      geom.pTheta = 0.033002;
      geom.pPhi = -0.319916;
      geom.pAlp1 = 0.019055;
      geom.pAlp2 = 0.019050;
      geom.pDy1 = 1.127101;
      geom.pDx1 = 1.187309;
      geom.pDx2 = 1.201654;
      geom.pDy2 = 1.266192;
      geom.pDx3 = 1.327568;
      geom.pDx4 = 1.343680;
      geom.centralX = 3.800447;
      geom.centralY = 104.189340;
      geom.centralZ = 64.494347;
      geom.pRotationAngleX = -1.026472;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1267 based Row/Col = 126/7
      geom_tower geom;
      geom.id = 1267;
      geom.pDz = 6.751789;
      geom.pTheta = 0.053232;
      geom.pPhi = -0.196109;
      geom.pAlp1 = 0.031808;
      geom.pAlp2 = 0.031800;
      geom.pDy1 = 1.127101;
      geom.pDx1 = 1.189362;
      geom.pDx2 = 1.203784;
      geom.pDy2 = 1.266192;
      geom.pDx3 = 1.329863;
      geom.pDx4 = 1.346061;
      geom.centralX = 6.337767;
      geom.centralY = 104.189340;
      geom.centralZ = 64.494347;
      geom.pRotationAngleX = -1.026472;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1268 based Row/Col = 126/8
      geom_tower geom;
      geom.id = 1268;
      geom.pDz = 6.751789;
      geom.pTheta = 0.073828;
      geom.pPhi = -0.140839;
      geom.pAlp1 = 0.044637;
      geom.pAlp2 = 0.044625;
      geom.pDy1 = 1.127101;
      geom.pDx1 = 1.192450;
      geom.pDx2 = 1.206989;
      geom.pDy2 = 1.266192;
      geom.pDx3 = 1.333316;
      geom.pDx4 = 1.349644;
      geom.centralX = 8.880634;
      geom.centralY = 104.189340;
      geom.centralZ = 64.494347;
      geom.pRotationAngleX = -1.026472;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1251 based Row/Col = 125/1
      geom_tower geom;
      geom.id = 1251;
      geom.pDz = 6.751768;
      geom.pTheta = 0.074709;
      geom.pPhi = 3.004647;
      geom.pAlp1 = -0.044636;
      geom.pAlp2 = -0.044624;
      geom.pDy1 = 1.127101;
      geom.pDx1 = 1.206989;
      geom.pDx2 = 1.221532;
      geom.pDy2 = 1.266192;
      geom.pDx3 = 1.349644;
      geom.pDx4 = 1.365976;
      geom.centralX = -8.987518;
      geom.centralY = 105.431273;
      geom.centralZ = 62.442662;
      geom.pRotationAngleX = -1.026472;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1252 based Row/Col = 125/2
      geom_tower geom;
      geom.id = 1252;
      geom.pDz = 6.751768;
      geom.pTheta = 0.053841;
      geom.pPhi = 2.950836;
      geom.pAlp1 = -0.031805;
      geom.pAlp2 = -0.031797;
      geom.pDy1 = 1.127101;
      geom.pDx1 = 1.203784;
      geom.pDx2 = 1.218207;
      geom.pDy2 = 1.266192;
      geom.pDx3 = 1.346061;
      geom.pDx4 = 1.362260;
      geom.centralX = -6.413905;
      geom.centralY = 105.431273;
      geom.centralZ = 62.442662;
      geom.pRotationAngleX = -1.026472;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1253 based Row/Col = 125/3
      geom_tower geom;
      geom.id = 1253;
      geom.pDz = 6.751768;
      geom.pTheta = 0.033323;
      geom.pPhi = 2.830046;
      geom.pAlp1 = -0.019052;
      geom.pAlp2 = -0.019047;
      geom.pDy1 = 1.127101;
      geom.pDx1 = 1.201654;
      geom.pDx2 = 1.215998;
      geom.pDy2 = 1.266192;
      geom.pDx3 = 1.343680;
      geom.pDx4 = 1.359789;
      geom.centralX = -3.846047;
      geom.centralY = 105.431273;
      geom.centralZ = 62.442662;
      geom.pRotationAngleX = -1.026472;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1254 based Row/Col = 125/4
      geom_tower geom;
      geom.id = 1254;
      geom.pDz = 6.751768;
      geom.pTheta = 0.014703;
      geom.pPhi = 2.373286;
      geom.pAlp1 = -0.006346;
      geom.pAlp2 = -0.006344;
      geom.pDy1 = 1.127101;
      geom.pDx1 = 1.200591;
      geom.pDx2 = 1.214895;
      geom.pDy2 = 1.266192;
      geom.pDx3 = 1.342491;
      geom.pDx4 = 1.358556;
      geom.centralX = -1.281633;
      geom.centralY = 105.431273;
      geom.centralZ = 62.442662;
      geom.pRotationAngleX = -1.026472;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1255 based Row/Col = 125/5
      geom_tower geom;
      geom.id = 1255;
      geom.pDz = 6.751768;
      geom.pTheta = 0.014703;
      geom.pPhi = 0.768307;
      geom.pAlp1 = 0.006346;
      geom.pAlp2 = 0.006344;
      geom.pDy1 = 1.127101;
      geom.pDx1 = 1.200591;
      geom.pDx2 = 1.214895;
      geom.pDy2 = 1.266192;
      geom.pDx3 = 1.342491;
      geom.pDx4 = 1.358556;
      geom.centralX = 1.281633;
      geom.centralY = 105.431273;
      geom.centralZ = 62.442662;
      geom.pRotationAngleX = -1.026472;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1256 based Row/Col = 125/6
      geom_tower geom;
      geom.id = 1256;
      geom.pDz = 6.751768;
      geom.pTheta = 0.033323;
      geom.pPhi = 0.311546;
      geom.pAlp1 = 0.019052;
      geom.pAlp2 = 0.019047;
      geom.pDy1 = 1.127101;
      geom.pDx1 = 1.201654;
      geom.pDx2 = 1.215998;
      geom.pDy2 = 1.266192;
      geom.pDx3 = 1.343680;
      geom.pDx4 = 1.359789;
      geom.centralX = 3.846047;
      geom.centralY = 105.431273;
      geom.centralZ = 62.442662;
      geom.pRotationAngleX = -1.026472;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1257 based Row/Col = 125/7
      geom_tower geom;
      geom.id = 1257;
      geom.pDz = 6.751768;
      geom.pTheta = 0.053841;
      geom.pPhi = 0.190757;
      geom.pAlp1 = 0.031805;
      geom.pAlp2 = 0.031797;
      geom.pDy1 = 1.127101;
      geom.pDx1 = 1.203784;
      geom.pDx2 = 1.218207;
      geom.pDy2 = 1.266192;
      geom.pDx3 = 1.346061;
      geom.pDx4 = 1.362260;
      geom.centralX = 6.413905;
      geom.centralY = 105.431273;
      geom.centralZ = 62.442662;
      geom.pRotationAngleX = -1.026472;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1258 based Row/Col = 125/8
      geom_tower geom;
      geom.id = 1258;
      geom.pDz = 6.751768;
      geom.pTheta = 0.074709;
      geom.pPhi = 0.136946;
      geom.pAlp1 = 0.044636;
      geom.pAlp2 = 0.044624;
      geom.pDy1 = 1.127101;
      geom.pDx1 = 1.206989;
      geom.pDx2 = 1.221532;
      geom.pDy2 = 1.266192;
      geom.pDx3 = 1.349644;
      geom.pDx4 = 1.365976;
      geom.centralX = 8.987518;
      geom.centralY = 105.431273;
      geom.centralZ = 62.442662;
      geom.pRotationAngleX = -1.026472;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 751 based Row/Col = 75/1
      geom_tower geom;
      geom.id = 751;
      geom.pDz = 6.751789;
      geom.pTheta = 0.073828;
      geom.pPhi = 3.000754;
      geom.pAlp1 = 0.044637;
      geom.pAlp2 = 0.044625;
      geom.pDy1 = 1.127101;
      geom.pDx1 = 1.206989;
      geom.pDx2 = 1.192450;
      geom.pDy2 = 1.266192;
      geom.pDx3 = 1.349644;
      geom.pDx4 = 1.333316;
      geom.centralX = -8.880634;
      geom.centralY = 104.189340;
      geom.centralZ = -64.494347;
      geom.pRotationAngleX = -2.115121;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 752 based Row/Col = 75/2
      geom_tower geom;
      geom.id = 752;
      geom.pDz = 6.751789;
      geom.pTheta = 0.053232;
      geom.pPhi = 2.945484;
      geom.pAlp1 = 0.031808;
      geom.pAlp2 = 0.031800;
      geom.pDy1 = 1.127101;
      geom.pDx1 = 1.203784;
      geom.pDx2 = 1.189362;
      geom.pDy2 = 1.266192;
      geom.pDx3 = 1.346061;
      geom.pDx4 = 1.329863;
      geom.centralX = -6.337767;
      geom.centralY = 104.189340;
      geom.centralZ = -64.494347;
      geom.pRotationAngleX = -2.115121;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 753 based Row/Col = 75/3
      geom_tower geom;
      geom.id = 753;
      geom.pDz = 6.751789;
      geom.pTheta = 0.033002;
      geom.pPhi = 2.821676;
      geom.pAlp1 = 0.019055;
      geom.pAlp2 = 0.019050;
      geom.pDy1 = 1.127101;
      geom.pDx1 = 1.201654;
      geom.pDx2 = 1.187309;
      geom.pDy2 = 1.266192;
      geom.pDx3 = 1.343680;
      geom.pDx4 = 1.327568;
      geom.centralX = -3.800447;
      geom.centralY = 104.189340;
      geom.centralZ = -64.494347;
      geom.pRotationAngleX = -2.115121;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 754 based Row/Col = 75/4
      geom_tower geom;
      geom.id = 754;
      geom.pDz = 6.751789;
      geom.pTheta = 0.014725;
      geom.pPhi = 2.359114;
      geom.pAlp1 = 0.006347;
      geom.pAlp2 = 0.006345;
      geom.pDy1 = 1.127101;
      geom.pDx1 = 1.200591;
      geom.pDx2 = 1.186284;
      geom.pDy2 = 1.266192;
      geom.pDx3 = 1.342491;
      geom.pDx4 = 1.326423;
      geom.centralX = -1.266447;
      geom.centralY = 104.189340;
      geom.centralZ = -64.494347;
      geom.pRotationAngleX = -2.115121;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 755 based Row/Col = 75/5
      geom_tower geom;
      geom.id = 755;
      geom.pDz = 6.751789;
      geom.pTheta = 0.014725;
      geom.pPhi = 0.782479;
      geom.pAlp1 = -0.006347;
      geom.pAlp2 = -0.006345;
      geom.pDy1 = 1.127101;
      geom.pDx1 = 1.200591;
      geom.pDx2 = 1.186284;
      geom.pDy2 = 1.266192;
      geom.pDx3 = 1.342491;
      geom.pDx4 = 1.326423;
      geom.centralX = 1.266447;
      geom.centralY = 104.189340;
      geom.centralZ = -64.494347;
      geom.pRotationAngleX = -2.115121;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 756 based Row/Col = 75/6
      geom_tower geom;
      geom.id = 756;
      geom.pDz = 6.751789;
      geom.pTheta = 0.033002;
      geom.pPhi = 0.319916;
      geom.pAlp1 = -0.019055;
      geom.pAlp2 = -0.019050;
      geom.pDy1 = 1.127101;
      geom.pDx1 = 1.201654;
      geom.pDx2 = 1.187309;
      geom.pDy2 = 1.266192;
      geom.pDx3 = 1.343680;
      geom.pDx4 = 1.327568;
      geom.centralX = 3.800447;
      geom.centralY = 104.189340;
      geom.centralZ = -64.494347;
      geom.pRotationAngleX = -2.115121;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 757 based Row/Col = 75/7
      geom_tower geom;
      geom.id = 757;
      geom.pDz = 6.751789;
      geom.pTheta = 0.053232;
      geom.pPhi = 0.196109;
      geom.pAlp1 = -0.031808;
      geom.pAlp2 = -0.031800;
      geom.pDy1 = 1.127101;
      geom.pDx1 = 1.203784;
      geom.pDx2 = 1.189362;
      geom.pDy2 = 1.266192;
      geom.pDx3 = 1.346061;
      geom.pDx4 = 1.329863;
      geom.centralX = 6.337767;
      geom.centralY = 104.189340;
      geom.centralZ = -64.494347;
      geom.pRotationAngleX = -2.115121;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 758 based Row/Col = 75/8
      geom_tower geom;
      geom.id = 758;
      geom.pDz = 6.751789;
      geom.pTheta = 0.073828;
      geom.pPhi = 0.140839;
      geom.pAlp1 = -0.044637;
      geom.pAlp2 = -0.044625;
      geom.pDy1 = 1.127101;
      geom.pDx1 = 1.206989;
      geom.pDx2 = 1.192450;
      geom.pDy2 = 1.266192;
      geom.pDx3 = 1.349644;
      geom.pDx4 = 1.333316;
      geom.centralX = 8.880634;
      geom.centralY = 104.189340;
      geom.centralZ = -64.494347;
      geom.pRotationAngleX = -2.115121;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 761 based Row/Col = 76/1
      geom_tower geom;
      geom.id = 761;
      geom.pDz = 6.751768;
      geom.pTheta = 0.074709;
      geom.pPhi = -3.004647;
      geom.pAlp1 = 0.044636;
      geom.pAlp2 = 0.044624;
      geom.pDy1 = 1.127101;
      geom.pDx1 = 1.221532;
      geom.pDx2 = 1.206989;
      geom.pDy2 = 1.266192;
      geom.pDx3 = 1.365976;
      geom.pDx4 = 1.349644;
      geom.centralX = -8.987518;
      geom.centralY = 105.431273;
      geom.centralZ = -62.442662;
      geom.pRotationAngleX = -2.115121;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 762 based Row/Col = 76/2
      geom_tower geom;
      geom.id = 762;
      geom.pDz = 6.751768;
      geom.pTheta = 0.053841;
      geom.pPhi = -2.950836;
      geom.pAlp1 = 0.031805;
      geom.pAlp2 = 0.031797;
      geom.pDy1 = 1.127101;
      geom.pDx1 = 1.218207;
      geom.pDx2 = 1.203784;
      geom.pDy2 = 1.266192;
      geom.pDx3 = 1.362260;
      geom.pDx4 = 1.346061;
      geom.centralX = -6.413905;
      geom.centralY = 105.431273;
      geom.centralZ = -62.442662;
      geom.pRotationAngleX = -2.115121;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 763 based Row/Col = 76/3
      geom_tower geom;
      geom.id = 763;
      geom.pDz = 6.751768;
      geom.pTheta = 0.033323;
      geom.pPhi = -2.830046;
      geom.pAlp1 = 0.019052;
      geom.pAlp2 = 0.019047;
      geom.pDy1 = 1.127101;
      geom.pDx1 = 1.215998;
      geom.pDx2 = 1.201654;
      geom.pDy2 = 1.266192;
      geom.pDx3 = 1.359789;
      geom.pDx4 = 1.343680;
      geom.centralX = -3.846047;
      geom.centralY = 105.431273;
      geom.centralZ = -62.442662;
      geom.pRotationAngleX = -2.115121;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 764 based Row/Col = 76/4
      geom_tower geom;
      geom.id = 764;
      geom.pDz = 6.751768;
      geom.pTheta = 0.014703;
      geom.pPhi = -2.373286;
      geom.pAlp1 = 0.006346;
      geom.pAlp2 = 0.006344;
      geom.pDy1 = 1.127101;
      geom.pDx1 = 1.214895;
      geom.pDx2 = 1.200591;
      geom.pDy2 = 1.266192;
      geom.pDx3 = 1.358556;
      geom.pDx4 = 1.342491;
      geom.centralX = -1.281633;
      geom.centralY = 105.431273;
      geom.centralZ = -62.442662;
      geom.pRotationAngleX = -2.115121;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 765 based Row/Col = 76/5
      geom_tower geom;
      geom.id = 765;
      geom.pDz = 6.751768;
      geom.pTheta = 0.014703;
      geom.pPhi = -0.768307;
      geom.pAlp1 = -0.006346;
      geom.pAlp2 = -0.006344;
      geom.pDy1 = 1.127101;
      geom.pDx1 = 1.214895;
      geom.pDx2 = 1.200591;
      geom.pDy2 = 1.266192;
      geom.pDx3 = 1.358556;
      geom.pDx4 = 1.342491;
      geom.centralX = 1.281633;
      geom.centralY = 105.431273;
      geom.centralZ = -62.442662;
      geom.pRotationAngleX = -2.115121;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 766 based Row/Col = 76/6
      geom_tower geom;
      geom.id = 766;
      geom.pDz = 6.751768;
      geom.pTheta = 0.033323;
      geom.pPhi = -0.311546;
      geom.pAlp1 = -0.019052;
      geom.pAlp2 = -0.019047;
      geom.pDy1 = 1.127101;
      geom.pDx1 = 1.215998;
      geom.pDx2 = 1.201654;
      geom.pDy2 = 1.266192;
      geom.pDx3 = 1.359789;
      geom.pDx4 = 1.343680;
      geom.centralX = 3.846047;
      geom.centralY = 105.431273;
      geom.centralZ = -62.442662;
      geom.pRotationAngleX = -2.115121;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 767 based Row/Col = 76/7
      geom_tower geom;
      geom.id = 767;
      geom.pDz = 6.751768;
      geom.pTheta = 0.053841;
      geom.pPhi = -0.190757;
      geom.pAlp1 = -0.031805;
      geom.pAlp2 = -0.031797;
      geom.pDy1 = 1.127101;
      geom.pDx1 = 1.218207;
      geom.pDx2 = 1.203784;
      geom.pDy2 = 1.266192;
      geom.pDx3 = 1.362260;
      geom.pDx4 = 1.346061;
      geom.centralX = 6.413905;
      geom.centralY = 105.431273;
      geom.centralZ = -62.442662;
      geom.pRotationAngleX = -2.115121;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 768 based Row/Col = 76/8
      geom_tower geom;
      geom.id = 768;
      geom.pDz = 6.751768;
      geom.pTheta = 0.074709;
      geom.pPhi = -0.136946;
      geom.pAlp1 = -0.044636;
      geom.pAlp2 = -0.044624;
      geom.pDy1 = 1.127101;
      geom.pDx1 = 1.221532;
      geom.pDx2 = 1.206989;
      geom.pDy2 = 1.266192;
      geom.pDx3 = 1.365976;
      geom.pDx4 = 1.349644;
      geom.centralX = 8.987518;
      geom.centralY = 105.431273;
      geom.centralZ = -62.442662;
      geom.pRotationAngleX = -2.115121;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1301 based Row/Col = 130/1
      geom_tower geom;
      geom.id = 1301;
      geom.pDz = 6.827976;
      geom.pTheta = 0.070180;
      geom.pPhi = -3.000649;
      geom.pAlp1 = -0.050110;
      geom.pAlp2 = -0.050096;
      geom.pDy1 = 1.127657;
      geom.pDx1 = 1.193847;
      geom.pDx2 = 1.210163;
      geom.pDy2 = 1.261082;
      geom.pDx3 = 1.329090;
      geom.pDx4 = 1.347331;
      geom.centralX = -8.880191;
      geom.centralY = 104.179846;
      geom.centralZ = 75.951969;
      geom.pRotationAngleX = -0.950408;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1302 based Row/Col = 130/2
      geom_tower geom;
      geom.id = 1302;
      geom.pDz = 6.827976;
      geom.pTheta = 0.050603;
      geom.pPhi = -2.945356;
      geom.pAlp1 = -0.035721;
      geom.pAlp2 = -0.035711;
      geom.pDy1 = 1.127657;
      geom.pDx1 = 1.191058;
      geom.pDx2 = 1.207255;
      geom.pDy2 = 1.261082;
      geom.pDx3 = 1.325985;
      geom.pDx4 = 1.344094;
      geom.centralX = -6.337986;
      geom.centralY = 104.179846;
      geom.centralZ = 75.951969;
      geom.pRotationAngleX = -0.950408;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1303 based Row/Col = 130/3
      geom_tower geom;
      geom.id = 1303;
      geom.pDz = 6.827976;
      geom.pTheta = 0.031372;
      geom.pPhi = -2.821493;
      geom.pAlp1 = -0.021403;
      geom.pAlp2 = -0.021398;
      geom.pDy1 = 1.127657;
      geom.pDx1 = 1.189203;
      geom.pDx2 = 1.205322;
      geom.pDy2 = 1.261082;
      geom.pDx3 = 1.323920;
      geom.pDx4 = 1.341942;
      geom.centralX = -3.800792;
      geom.centralY = 104.179846;
      geom.centralZ = 75.951969;
      geom.pRotationAngleX = -0.950408;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1304 based Row/Col = 130/4
      geom_tower geom;
      geom.id = 1304;
      geom.pDz = 6.827976;
      geom.pTheta = 0.014001;
      geom.pPhi = -2.358821;
      geom.pAlp1 = -0.007130;
      geom.pAlp2 = -0.007128;
      geom.pDy1 = 1.127657;
      geom.pDx1 = 1.188277;
      geom.pDx2 = 1.204357;
      geom.pDy2 = 1.261082;
      geom.pDx3 = 1.322890;
      geom.pDx4 = 1.340867;
      geom.centralX = -1.266597;
      geom.centralY = 104.179846;
      geom.centralZ = 75.951969;
      geom.pRotationAngleX = -0.950408;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1305 based Row/Col = 130/5
      geom_tower geom;
      geom.id = 1305;
      geom.pDz = 6.827976;
      geom.pTheta = 0.014001;
      geom.pPhi = -0.782771;
      geom.pAlp1 = 0.007130;
      geom.pAlp2 = 0.007128;
      geom.pDy1 = 1.127657;
      geom.pDx1 = 1.188277;
      geom.pDx2 = 1.204357;
      geom.pDy2 = 1.261082;
      geom.pDx3 = 1.322890;
      geom.pDx4 = 1.340867;
      geom.centralX = 1.266597;
      geom.centralY = 104.179846;
      geom.centralZ = 75.951969;
      geom.pRotationAngleX = -0.950408;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1306 based Row/Col = 130/6
      geom_tower geom;
      geom.id = 1306;
      geom.pDz = 6.827976;
      geom.pTheta = 0.031372;
      geom.pPhi = -0.320099;
      geom.pAlp1 = 0.021403;
      geom.pAlp2 = 0.021398;
      geom.pDy1 = 1.127657;
      geom.pDx1 = 1.189203;
      geom.pDx2 = 1.205322;
      geom.pDy2 = 1.261082;
      geom.pDx3 = 1.323920;
      geom.pDx4 = 1.341942;
      geom.centralX = 3.800792;
      geom.centralY = 104.179846;
      geom.centralZ = 75.951969;
      geom.pRotationAngleX = -0.950408;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1307 based Row/Col = 130/7
      geom_tower geom;
      geom.id = 1307;
      geom.pDz = 6.827976;
      geom.pTheta = 0.050603;
      geom.pPhi = -0.196237;
      geom.pAlp1 = 0.035721;
      geom.pAlp2 = 0.035711;
      geom.pDy1 = 1.127657;
      geom.pDx1 = 1.191058;
      geom.pDx2 = 1.207255;
      geom.pDy2 = 1.261082;
      geom.pDx3 = 1.325985;
      geom.pDx4 = 1.344094;
      geom.centralX = 6.337986;
      geom.centralY = 104.179846;
      geom.centralZ = 75.951969;
      geom.pRotationAngleX = -0.950408;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1308 based Row/Col = 130/8
      geom_tower geom;
      geom.id = 1308;
      geom.pDz = 6.827976;
      geom.pTheta = 0.070180;
      geom.pPhi = -0.140944;
      geom.pAlp1 = 0.050110;
      geom.pAlp2 = 0.050096;
      geom.pDy1 = 1.127657;
      geom.pDx1 = 1.193847;
      geom.pDx2 = 1.210163;
      geom.pDy2 = 1.261082;
      geom.pDx3 = 1.329090;
      geom.pDx4 = 1.347331;
      geom.centralX = 8.880191;
      geom.centralY = 104.179846;
      geom.centralZ = 75.951969;
      geom.pRotationAngleX = -0.950408;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1291 based Row/Col = 129/1
      geom_tower geom;
      geom.id = 1291;
      geom.pDz = 6.827935;
      geom.pTheta = 0.071115;
      geom.pPhi = 3.005481;
      geom.pAlp1 = -0.050108;
      geom.pAlp2 = -0.050094;
      geom.pDy1 = 1.127657;
      geom.pDx1 = 1.210163;
      geom.pDx2 = 1.226483;
      geom.pDy2 = 1.261082;
      geom.pDx3 = 1.347331;
      geom.pDx4 = 1.365576;
      geom.centralX = -8.999972;
      geom.centralY = 105.571449;
      geom.centralZ = 74.004297;
      geom.pRotationAngleX = -0.950408;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1292 based Row/Col = 129/2
      geom_tower geom;
      geom.id = 1292;
      geom.pDz = 6.827935;
      geom.pTheta = 0.051245;
      geom.pPhi = 2.952000;
      geom.pAlp1 = -0.035717;
      geom.pAlp2 = -0.035707;
      geom.pDy1 = 1.127657;
      geom.pDx1 = 1.207255;
      geom.pDx2 = 1.223453;
      geom.pDy2 = 1.261082;
      geom.pDx3 = 1.344094;
      geom.pDx4 = 1.362203;
      geom.centralX = -6.423333;
      geom.centralY = 105.571449;
      geom.centralZ = 74.004297;
      geom.pRotationAngleX = -0.950408;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1293 based Row/Col = 129/3
      geom_tower geom;
      geom.id = 1293;
      geom.pDz = 6.827935;
      geom.pTheta = 0.031704;
      geom.pPhi = 2.831887;
      geom.pAlp1 = -0.021400;
      geom.pAlp2 = -0.021395;
      geom.pDy1 = 1.127657;
      geom.pDx1 = 1.205322;
      geom.pDx2 = 1.221439;
      geom.pDy2 = 1.261082;
      geom.pDx3 = 1.341942;
      geom.pDx4 = 1.359961;
      geom.centralX = -3.851915;
      geom.centralY = 105.571449;
      geom.centralZ = 74.004297;
      geom.pRotationAngleX = -0.950408;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1294 based Row/Col = 129/4
      geom_tower geom;
      geom.id = 1294;
      geom.pDz = 6.827935;
      geom.pTheta = 0.013955;
      geom.pPhi = 2.376460;
      geom.pAlp1 = -0.007128;
      geom.pAlp2 = -0.007127;
      geom.pDy1 = 1.127657;
      geom.pDx1 = 1.204357;
      geom.pDx2 = 1.220434;
      geom.pDy2 = 1.261082;
      geom.pDx3 = 1.340867;
      geom.pDx4 = 1.358842;
      geom.centralX = -1.283625;
      geom.centralY = 105.571449;
      geom.centralZ = 74.004297;
      geom.pRotationAngleX = -0.950408;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1295 based Row/Col = 129/5
      geom_tower geom;
      geom.id = 1295;
      geom.pDz = 6.827935;
      geom.pTheta = 0.013955;
      geom.pPhi = 0.765132;
      geom.pAlp1 = 0.007128;
      geom.pAlp2 = 0.007127;
      geom.pDy1 = 1.127657;
      geom.pDx1 = 1.204357;
      geom.pDx2 = 1.220434;
      geom.pDy2 = 1.261082;
      geom.pDx3 = 1.340867;
      geom.pDx4 = 1.358842;
      geom.centralX = 1.283625;
      geom.centralY = 105.571449;
      geom.centralZ = 74.004297;
      geom.pRotationAngleX = -0.950408;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1296 based Row/Col = 129/6
      geom_tower geom;
      geom.id = 1296;
      geom.pDz = 6.827935;
      geom.pTheta = 0.031704;
      geom.pPhi = 0.309706;
      geom.pAlp1 = 0.021400;
      geom.pAlp2 = 0.021395;
      geom.pDy1 = 1.127657;
      geom.pDx1 = 1.205322;
      geom.pDx2 = 1.221439;
      geom.pDy2 = 1.261082;
      geom.pDx3 = 1.341942;
      geom.pDx4 = 1.359961;
      geom.centralX = 3.851915;
      geom.centralY = 105.571449;
      geom.centralZ = 74.004297;
      geom.pRotationAngleX = -0.950408;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1297 based Row/Col = 129/7
      geom_tower geom;
      geom.id = 1297;
      geom.pDz = 6.827935;
      geom.pTheta = 0.051245;
      geom.pPhi = 0.189593;
      geom.pAlp1 = 0.035717;
      geom.pAlp2 = 0.035707;
      geom.pDy1 = 1.127657;
      geom.pDx1 = 1.207255;
      geom.pDx2 = 1.223453;
      geom.pDy2 = 1.261082;
      geom.pDx3 = 1.344094;
      geom.pDx4 = 1.362203;
      geom.centralX = 6.423333;
      geom.centralY = 105.571449;
      geom.centralZ = 74.004297;
      geom.pRotationAngleX = -0.950408;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1298 based Row/Col = 129/8
      geom_tower geom;
      geom.id = 1298;
      geom.pDz = 6.827935;
      geom.pTheta = 0.071115;
      geom.pPhi = 0.136112;
      geom.pAlp1 = 0.050108;
      geom.pAlp2 = 0.050094;
      geom.pDy1 = 1.127657;
      geom.pDx1 = 1.210163;
      geom.pDx2 = 1.226483;
      geom.pDy2 = 1.261082;
      geom.pDx3 = 1.347331;
      geom.pDx4 = 1.365576;
      geom.centralX = 8.999972;
      geom.centralY = 105.571449;
      geom.centralZ = 74.004297;
      geom.pRotationAngleX = -0.950408;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 711 based Row/Col = 71/1
      geom_tower geom;
      geom.id = 711;
      geom.pDz = 6.827976;
      geom.pTheta = 0.070180;
      geom.pPhi = 3.000649;
      geom.pAlp1 = 0.050110;
      geom.pAlp2 = 0.050096;
      geom.pDy1 = 1.127657;
      geom.pDx1 = 1.210163;
      geom.pDx2 = 1.193847;
      geom.pDy2 = 1.261082;
      geom.pDx3 = 1.347331;
      geom.pDx4 = 1.329090;
      geom.centralX = -8.880191;
      geom.centralY = 104.179846;
      geom.centralZ = -75.951969;
      geom.pRotationAngleX = -2.191185;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 712 based Row/Col = 71/2
      geom_tower geom;
      geom.id = 712;
      geom.pDz = 6.827976;
      geom.pTheta = 0.050603;
      geom.pPhi = 2.945356;
      geom.pAlp1 = 0.035721;
      geom.pAlp2 = 0.035711;
      geom.pDy1 = 1.127657;
      geom.pDx1 = 1.207255;
      geom.pDx2 = 1.191058;
      geom.pDy2 = 1.261082;
      geom.pDx3 = 1.344094;
      geom.pDx4 = 1.325985;
      geom.centralX = -6.337986;
      geom.centralY = 104.179846;
      geom.centralZ = -75.951969;
      geom.pRotationAngleX = -2.191185;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 713 based Row/Col = 71/3
      geom_tower geom;
      geom.id = 713;
      geom.pDz = 6.827976;
      geom.pTheta = 0.031372;
      geom.pPhi = 2.821493;
      geom.pAlp1 = 0.021403;
      geom.pAlp2 = 0.021398;
      geom.pDy1 = 1.127657;
      geom.pDx1 = 1.205322;
      geom.pDx2 = 1.189203;
      geom.pDy2 = 1.261082;
      geom.pDx3 = 1.341942;
      geom.pDx4 = 1.323920;
      geom.centralX = -3.800792;
      geom.centralY = 104.179846;
      geom.centralZ = -75.951969;
      geom.pRotationAngleX = -2.191185;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 714 based Row/Col = 71/4
      geom_tower geom;
      geom.id = 714;
      geom.pDz = 6.827976;
      geom.pTheta = 0.014001;
      geom.pPhi = 2.358821;
      geom.pAlp1 = 0.007130;
      geom.pAlp2 = 0.007128;
      geom.pDy1 = 1.127657;
      geom.pDx1 = 1.204357;
      geom.pDx2 = 1.188277;
      geom.pDy2 = 1.261082;
      geom.pDx3 = 1.340867;
      geom.pDx4 = 1.322890;
      geom.centralX = -1.266597;
      geom.centralY = 104.179846;
      geom.centralZ = -75.951969;
      geom.pRotationAngleX = -2.191185;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 715 based Row/Col = 71/5
      geom_tower geom;
      geom.id = 715;
      geom.pDz = 6.827976;
      geom.pTheta = 0.014001;
      geom.pPhi = 0.782771;
      geom.pAlp1 = -0.007130;
      geom.pAlp2 = -0.007128;
      geom.pDy1 = 1.127657;
      geom.pDx1 = 1.204357;
      geom.pDx2 = 1.188277;
      geom.pDy2 = 1.261082;
      geom.pDx3 = 1.340867;
      geom.pDx4 = 1.322890;
      geom.centralX = 1.266597;
      geom.centralY = 104.179846;
      geom.centralZ = -75.951969;
      geom.pRotationAngleX = -2.191185;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 716 based Row/Col = 71/6
      geom_tower geom;
      geom.id = 716;
      geom.pDz = 6.827976;
      geom.pTheta = 0.031372;
      geom.pPhi = 0.320099;
      geom.pAlp1 = -0.021403;
      geom.pAlp2 = -0.021398;
      geom.pDy1 = 1.127657;
      geom.pDx1 = 1.205322;
      geom.pDx2 = 1.189203;
      geom.pDy2 = 1.261082;
      geom.pDx3 = 1.341942;
      geom.pDx4 = 1.323920;
      geom.centralX = 3.800792;
      geom.centralY = 104.179846;
      geom.centralZ = -75.951969;
      geom.pRotationAngleX = -2.191185;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 717 based Row/Col = 71/7
      geom_tower geom;
      geom.id = 717;
      geom.pDz = 6.827976;
      geom.pTheta = 0.050603;
      geom.pPhi = 0.196237;
      geom.pAlp1 = -0.035721;
      geom.pAlp2 = -0.035711;
      geom.pDy1 = 1.127657;
      geom.pDx1 = 1.207255;
      geom.pDx2 = 1.191058;
      geom.pDy2 = 1.261082;
      geom.pDx3 = 1.344094;
      geom.pDx4 = 1.325985;
      geom.centralX = 6.337986;
      geom.centralY = 104.179846;
      geom.centralZ = -75.951969;
      geom.pRotationAngleX = -2.191185;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 718 based Row/Col = 71/8
      geom_tower geom;
      geom.id = 718;
      geom.pDz = 6.827976;
      geom.pTheta = 0.070180;
      geom.pPhi = 0.140944;
      geom.pAlp1 = -0.050110;
      geom.pAlp2 = -0.050096;
      geom.pDy1 = 1.127657;
      geom.pDx1 = 1.210163;
      geom.pDx2 = 1.193847;
      geom.pDy2 = 1.261082;
      geom.pDx3 = 1.347331;
      geom.pDx4 = 1.329090;
      geom.centralX = 8.880191;
      geom.centralY = 104.179846;
      geom.centralZ = -75.951969;
      geom.pRotationAngleX = -2.191185;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 721 based Row/Col = 72/1
      geom_tower geom;
      geom.id = 721;
      geom.pDz = 6.827935;
      geom.pTheta = 0.071115;
      geom.pPhi = -3.005481;
      geom.pAlp1 = 0.050108;
      geom.pAlp2 = 0.050094;
      geom.pDy1 = 1.127657;
      geom.pDx1 = 1.226483;
      geom.pDx2 = 1.210163;
      geom.pDy2 = 1.261082;
      geom.pDx3 = 1.365576;
      geom.pDx4 = 1.347331;
      geom.centralX = -8.999972;
      geom.centralY = 105.571449;
      geom.centralZ = -74.004297;
      geom.pRotationAngleX = -2.191185;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 722 based Row/Col = 72/2
      geom_tower geom;
      geom.id = 722;
      geom.pDz = 6.827935;
      geom.pTheta = 0.051245;
      geom.pPhi = -2.952000;
      geom.pAlp1 = 0.035717;
      geom.pAlp2 = 0.035707;
      geom.pDy1 = 1.127657;
      geom.pDx1 = 1.223453;
      geom.pDx2 = 1.207255;
      geom.pDy2 = 1.261082;
      geom.pDx3 = 1.362203;
      geom.pDx4 = 1.344094;
      geom.centralX = -6.423333;
      geom.centralY = 105.571449;
      geom.centralZ = -74.004297;
      geom.pRotationAngleX = -2.191185;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 723 based Row/Col = 72/3
      geom_tower geom;
      geom.id = 723;
      geom.pDz = 6.827935;
      geom.pTheta = 0.031704;
      geom.pPhi = -2.831887;
      geom.pAlp1 = 0.021400;
      geom.pAlp2 = 0.021395;
      geom.pDy1 = 1.127657;
      geom.pDx1 = 1.221439;
      geom.pDx2 = 1.205322;
      geom.pDy2 = 1.261082;
      geom.pDx3 = 1.359961;
      geom.pDx4 = 1.341942;
      geom.centralX = -3.851915;
      geom.centralY = 105.571449;
      geom.centralZ = -74.004297;
      geom.pRotationAngleX = -2.191185;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 724 based Row/Col = 72/4
      geom_tower geom;
      geom.id = 724;
      geom.pDz = 6.827935;
      geom.pTheta = 0.013955;
      geom.pPhi = -2.376460;
      geom.pAlp1 = 0.007128;
      geom.pAlp2 = 0.007127;
      geom.pDy1 = 1.127657;
      geom.pDx1 = 1.220434;
      geom.pDx2 = 1.204357;
      geom.pDy2 = 1.261082;
      geom.pDx3 = 1.358842;
      geom.pDx4 = 1.340867;
      geom.centralX = -1.283625;
      geom.centralY = 105.571449;
      geom.centralZ = -74.004297;
      geom.pRotationAngleX = -2.191185;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 725 based Row/Col = 72/5
      geom_tower geom;
      geom.id = 725;
      geom.pDz = 6.827935;
      geom.pTheta = 0.013955;
      geom.pPhi = -0.765132;
      geom.pAlp1 = -0.007128;
      geom.pAlp2 = -0.007127;
      geom.pDy1 = 1.127657;
      geom.pDx1 = 1.220434;
      geom.pDx2 = 1.204357;
      geom.pDy2 = 1.261082;
      geom.pDx3 = 1.358842;
      geom.pDx4 = 1.340867;
      geom.centralX = 1.283625;
      geom.centralY = 105.571449;
      geom.centralZ = -74.004297;
      geom.pRotationAngleX = -2.191185;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 726 based Row/Col = 72/6
      geom_tower geom;
      geom.id = 726;
      geom.pDz = 6.827935;
      geom.pTheta = 0.031704;
      geom.pPhi = -0.309706;
      geom.pAlp1 = -0.021400;
      geom.pAlp2 = -0.021395;
      geom.pDy1 = 1.127657;
      geom.pDx1 = 1.221439;
      geom.pDx2 = 1.205322;
      geom.pDy2 = 1.261082;
      geom.pDx3 = 1.359961;
      geom.pDx4 = 1.341942;
      geom.centralX = 3.851915;
      geom.centralY = 105.571449;
      geom.centralZ = -74.004297;
      geom.pRotationAngleX = -2.191185;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 727 based Row/Col = 72/7
      geom_tower geom;
      geom.id = 727;
      geom.pDz = 6.827935;
      geom.pTheta = 0.051245;
      geom.pPhi = -0.189593;
      geom.pAlp1 = -0.035717;
      geom.pAlp2 = -0.035707;
      geom.pDy1 = 1.127657;
      geom.pDx1 = 1.223453;
      geom.pDx2 = 1.207255;
      geom.pDy2 = 1.261082;
      geom.pDx3 = 1.362203;
      geom.pDx4 = 1.344094;
      geom.centralX = 6.423333;
      geom.centralY = 105.571449;
      geom.centralZ = -74.004297;
      geom.pRotationAngleX = -2.191185;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 728 based Row/Col = 72/8
      geom_tower geom;
      geom.id = 728;
      geom.pDz = 6.827935;
      geom.pTheta = 0.071115;
      geom.pPhi = -0.136112;
      geom.pAlp1 = -0.050108;
      geom.pAlp2 = -0.050094;
      geom.pDy1 = 1.127657;
      geom.pDx1 = 1.226483;
      geom.pDx2 = 1.210163;
      geom.pDy2 = 1.261082;
      geom.pDx3 = 1.365576;
      geom.pDx4 = 1.347331;
      geom.centralX = 8.999972;
      geom.centralY = 105.571449;
      geom.centralZ = -74.004297;
      geom.pRotationAngleX = -2.191185;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1341 based Row/Col = 134/1
      geom_tower geom;
      geom.id = 1341;
      geom.pDz = 7.050152;
      geom.pTheta = 0.066350;
      geom.pPhi = -3.001610;
      geom.pAlp1 = -0.055028;
      geom.pAlp2 = -0.055013;
      geom.pDy1 = 1.125572;
      geom.pDx1 = 1.193452;
      geom.pDx2 = 1.211320;
      geom.pDy2 = 1.254731;
      geom.pDx3 = 1.325341;
      geom.pDx4 = 1.345254;
      geom.centralX = -8.874324;
      geom.centralY = 104.107425;
      geom.centralZ = 87.947881;
      geom.pRotationAngleX = -0.878335;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1342 based Row/Col = 134/2
      geom_tower geom;
      geom.id = 1342;
      geom.pDz = 7.050152;
      geom.pTheta = 0.047835;
      geom.pPhi = -2.946693;
      geom.pAlp1 = -0.039240;
      geom.pAlp2 = -0.039230;
      geom.pDy1 = 1.125572;
      geom.pDx1 = 1.190961;
      geom.pDx2 = 1.208714;
      geom.pDy2 = 1.254731;
      geom.pDx3 = 1.322577;
      geom.pDx4 = 1.342361;
      geom.centralX = -6.334329;
      geom.centralY = 104.107425;
      geom.centralZ = 87.947881;
      geom.pRotationAngleX = -0.878335;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1343 based Row/Col = 134/3
      geom_tower geom;
      geom.id = 1343;
      geom.pDz = 7.050152;
      geom.pTheta = 0.029644;
      geom.pPhi = -2.823599;
      geom.pAlp1 = -0.023518;
      geom.pAlp2 = -0.023512;
      geom.pDy1 = 1.125572;
      geom.pDx1 = 1.189305;
      geom.pDx2 = 1.206981;
      geom.pDy2 = 1.254731;
      geom.pDx3 = 1.320738;
      geom.pDx4 = 1.340437;
      geom.centralX = -3.798810;
      geom.centralY = 104.107425;
      geom.centralZ = 87.947881;
      geom.pRotationAngleX = -0.878335;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1344 based Row/Col = 134/4
      geom_tower geom;
      geom.id = 1344;
      geom.pDz = 7.050152;
      geom.pTheta = 0.013192;
      geom.pPhi = -2.362369;
      geom.pAlp1 = -0.007835;
      geom.pAlp2 = -0.007833;
      geom.pDy1 = 1.125572;
      geom.pDx1 = 1.188478;
      geom.pDx2 = 1.206116;
      geom.pDy2 = 1.254731;
      geom.pDx3 = 1.319820;
      geom.pDx4 = 1.339476;
      geom.centralX = -1.265973;
      geom.centralY = 104.107425;
      geom.centralZ = 87.947881;
      geom.pRotationAngleX = -0.878335;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1345 based Row/Col = 134/5
      geom_tower geom;
      geom.id = 1345;
      geom.pDz = 7.050152;
      geom.pTheta = 0.013192;
      geom.pPhi = -0.779223;
      geom.pAlp1 = 0.007835;
      geom.pAlp2 = 0.007833;
      geom.pDy1 = 1.125572;
      geom.pDx1 = 1.188478;
      geom.pDx2 = 1.206116;
      geom.pDy2 = 1.254731;
      geom.pDx3 = 1.319820;
      geom.pDx4 = 1.339476;
      geom.centralX = 1.265973;
      geom.centralY = 104.107425;
      geom.centralZ = 87.947881;
      geom.pRotationAngleX = -0.878335;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1346 based Row/Col = 134/6
      geom_tower geom;
      geom.id = 1346;
      geom.pDz = 7.050152;
      geom.pTheta = 0.029644;
      geom.pPhi = -0.317994;
      geom.pAlp1 = 0.023518;
      geom.pAlp2 = 0.023512;
      geom.pDy1 = 1.125572;
      geom.pDx1 = 1.189305;
      geom.pDx2 = 1.206981;
      geom.pDy2 = 1.254731;
      geom.pDx3 = 1.320738;
      geom.pDx4 = 1.340437;
      geom.centralX = 3.798810;
      geom.centralY = 104.107425;
      geom.centralZ = 87.947881;
      geom.pRotationAngleX = -0.878335;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1347 based Row/Col = 134/7
      geom_tower geom;
      geom.id = 1347;
      geom.pDz = 7.050152;
      geom.pTheta = 0.047835;
      geom.pPhi = -0.194900;
      geom.pAlp1 = 0.039240;
      geom.pAlp2 = 0.039230;
      geom.pDy1 = 1.125572;
      geom.pDx1 = 1.190961;
      geom.pDx2 = 1.208714;
      geom.pDy2 = 1.254731;
      geom.pDx3 = 1.322577;
      geom.pDx4 = 1.342361;
      geom.centralX = 6.334329;
      geom.centralY = 104.107425;
      geom.centralZ = 87.947881;
      geom.pRotationAngleX = -0.878335;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1348 based Row/Col = 134/8
      geom_tower geom;
      geom.id = 1348;
      geom.pDz = 7.050152;
      geom.pTheta = 0.066350;
      geom.pPhi = -0.139983;
      geom.pAlp1 = 0.055028;
      geom.pAlp2 = 0.055013;
      geom.pDy1 = 1.125572;
      geom.pDx1 = 1.193452;
      geom.pDx2 = 1.211320;
      geom.pDy2 = 1.254731;
      geom.pDx3 = 1.325341;
      geom.pDx4 = 1.345254;
      geom.centralX = 8.874324;
      geom.centralY = 104.107425;
      geom.centralZ = 87.947881;
      geom.pRotationAngleX = -0.878335;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1331 based Row/Col = 133/1
      geom_tower geom;
      geom.id = 1331;
      geom.pDz = 7.050105;
      geom.pTheta = 0.067313;
      geom.pPhi = 3.006958;
      geom.pAlp1 = -0.055027;
      geom.pAlp2 = -0.055012;
      geom.pDy1 = 1.125572;
      geom.pDx1 = 1.211320;
      geom.pDx2 = 1.229192;
      geom.pDy2 = 1.254731;
      geom.pDx3 = 1.345254;
      geom.pDx4 = 1.365172;
      geom.centralX = -9.005420;
      geom.centralY = 105.630282;
      geom.centralZ = 86.111969;
      geom.pRotationAngleX = -0.878335;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1332 based Row/Col = 133/2
      geom_tower geom;
      geom.id = 1332;
      geom.pDz = 7.050105;
      geom.pTheta = 0.048496;
      geom.pPhi = 2.954049;
      geom.pAlp1 = -0.039236;
      geom.pAlp2 = -0.039226;
      geom.pDy1 = 1.125572;
      geom.pDx1 = 1.208714;
      geom.pDx2 = 1.226468;
      geom.pDy2 = 1.254731;
      geom.pDx3 = 1.342361;
      geom.pDx4 = 1.362146;
      geom.centralX = -6.427763;
      geom.centralY = 105.630282;
      geom.centralZ = 86.111969;
      geom.pRotationAngleX = -0.878335;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1333 based Row/Col = 133/3
      geom_tower geom;
      geom.id = 1333;
      geom.pDz = 7.050105;
      geom.pTheta = 0.029984;
      geom.pPhi = 2.835119;
      geom.pAlp1 = -0.023515;
      geom.pAlp2 = -0.023508;
      geom.pDy1 = 1.125572;
      geom.pDx1 = 1.206981;
      geom.pDx2 = 1.224656;
      geom.pDy2 = 1.254731;
      geom.pDx3 = 1.340437;
      geom.pDx4 = 1.360134;
      geom.centralX = -3.854789;
      geom.centralY = 105.630282;
      geom.centralZ = 86.111969;
      geom.pRotationAngleX = -0.878335;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1334 based Row/Col = 133/4
      geom_tower geom;
      geom.id = 1334;
      geom.pDz = 7.050105;
      geom.pTheta = 0.013141;
      geom.pPhi = 2.382061;
      geom.pAlp1 = -0.007834;
      geom.pAlp2 = -0.007832;
      geom.pDy1 = 1.125572;
      geom.pDx1 = 1.206116;
      geom.pDx2 = 1.223751;
      geom.pDy2 = 1.254731;
      geom.pDx3 = 1.339476;
      geom.pDx4 = 1.359130;
      geom.centralX = -1.284618;
      geom.centralY = 105.630282;
      geom.centralZ = 86.111969;
      geom.pRotationAngleX = -0.878335;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1335 based Row/Col = 133/5
      geom_tower geom;
      geom.id = 1335;
      geom.pDz = 7.050105;
      geom.pTheta = 0.013141;
      geom.pPhi = 0.759532;
      geom.pAlp1 = 0.007834;
      geom.pAlp2 = 0.007832;
      geom.pDy1 = 1.125572;
      geom.pDx1 = 1.206116;
      geom.pDx2 = 1.223751;
      geom.pDy2 = 1.254731;
      geom.pDx3 = 1.339476;
      geom.pDx4 = 1.359130;
      geom.centralX = 1.284618;
      geom.centralY = 105.630282;
      geom.centralZ = 86.111969;
      geom.pRotationAngleX = -0.878335;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1336 based Row/Col = 133/6
      geom_tower geom;
      geom.id = 1336;
      geom.pDz = 7.050105;
      geom.pTheta = 0.029984;
      geom.pPhi = 0.306474;
      geom.pAlp1 = 0.023515;
      geom.pAlp2 = 0.023508;
      geom.pDy1 = 1.125572;
      geom.pDx1 = 1.206981;
      geom.pDx2 = 1.224656;
      geom.pDy2 = 1.254731;
      geom.pDx3 = 1.340437;
      geom.pDx4 = 1.360134;
      geom.centralX = 3.854789;
      geom.centralY = 105.630282;
      geom.centralZ = 86.111969;
      geom.pRotationAngleX = -0.878335;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1337 based Row/Col = 133/7
      geom_tower geom;
      geom.id = 1337;
      geom.pDz = 7.050105;
      geom.pTheta = 0.048496;
      geom.pPhi = 0.187544;
      geom.pAlp1 = 0.039236;
      geom.pAlp2 = 0.039226;
      geom.pDy1 = 1.125572;
      geom.pDx1 = 1.208714;
      geom.pDx2 = 1.226468;
      geom.pDy2 = 1.254731;
      geom.pDx3 = 1.342361;
      geom.pDx4 = 1.362146;
      geom.centralX = 6.427763;
      geom.centralY = 105.630282;
      geom.centralZ = 86.111969;
      geom.pRotationAngleX = -0.878335;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1338 based Row/Col = 133/8
      geom_tower geom;
      geom.id = 1338;
      geom.pDz = 7.050105;
      geom.pTheta = 0.067313;
      geom.pPhi = 0.134635;
      geom.pAlp1 = 0.055027;
      geom.pAlp2 = 0.055012;
      geom.pDy1 = 1.125572;
      geom.pDx1 = 1.211320;
      geom.pDx2 = 1.229192;
      geom.pDy2 = 1.254731;
      geom.pDx3 = 1.345254;
      geom.pDx4 = 1.365172;
      geom.centralX = 9.005420;
      geom.centralY = 105.630282;
      geom.centralZ = 86.111969;
      geom.pRotationAngleX = -0.878335;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 671 based Row/Col = 67/1
      geom_tower geom;
      geom.id = 671;
      geom.pDz = 7.050152;
      geom.pTheta = 0.066350;
      geom.pPhi = 3.001610;
      geom.pAlp1 = 0.055028;
      geom.pAlp2 = 0.055013;
      geom.pDy1 = 1.125572;
      geom.pDx1 = 1.211320;
      geom.pDx2 = 1.193452;
      geom.pDy2 = 1.254731;
      geom.pDx3 = 1.345254;
      geom.pDx4 = 1.325341;
      geom.centralX = -8.874324;
      geom.centralY = 104.107425;
      geom.centralZ = -87.947881;
      geom.pRotationAngleX = -2.263258;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 672 based Row/Col = 67/2
      geom_tower geom;
      geom.id = 672;
      geom.pDz = 7.050152;
      geom.pTheta = 0.047835;
      geom.pPhi = 2.946693;
      geom.pAlp1 = 0.039240;
      geom.pAlp2 = 0.039230;
      geom.pDy1 = 1.125572;
      geom.pDx1 = 1.208714;
      geom.pDx2 = 1.190961;
      geom.pDy2 = 1.254731;
      geom.pDx3 = 1.342361;
      geom.pDx4 = 1.322577;
      geom.centralX = -6.334329;
      geom.centralY = 104.107425;
      geom.centralZ = -87.947881;
      geom.pRotationAngleX = -2.263258;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 673 based Row/Col = 67/3
      geom_tower geom;
      geom.id = 673;
      geom.pDz = 7.050152;
      geom.pTheta = 0.029644;
      geom.pPhi = 2.823599;
      geom.pAlp1 = 0.023518;
      geom.pAlp2 = 0.023512;
      geom.pDy1 = 1.125572;
      geom.pDx1 = 1.206981;
      geom.pDx2 = 1.189305;
      geom.pDy2 = 1.254731;
      geom.pDx3 = 1.340437;
      geom.pDx4 = 1.320738;
      geom.centralX = -3.798810;
      geom.centralY = 104.107425;
      geom.centralZ = -87.947881;
      geom.pRotationAngleX = -2.263258;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 674 based Row/Col = 67/4
      geom_tower geom;
      geom.id = 674;
      geom.pDz = 7.050152;
      geom.pTheta = 0.013192;
      geom.pPhi = 2.362369;
      geom.pAlp1 = 0.007835;
      geom.pAlp2 = 0.007833;
      geom.pDy1 = 1.125572;
      geom.pDx1 = 1.206116;
      geom.pDx2 = 1.188478;
      geom.pDy2 = 1.254731;
      geom.pDx3 = 1.339476;
      geom.pDx4 = 1.319820;
      geom.centralX = -1.265973;
      geom.centralY = 104.107425;
      geom.centralZ = -87.947881;
      geom.pRotationAngleX = -2.263258;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 675 based Row/Col = 67/5
      geom_tower geom;
      geom.id = 675;
      geom.pDz = 7.050152;
      geom.pTheta = 0.013192;
      geom.pPhi = 0.779223;
      geom.pAlp1 = -0.007835;
      geom.pAlp2 = -0.007833;
      geom.pDy1 = 1.125572;
      geom.pDx1 = 1.206116;
      geom.pDx2 = 1.188478;
      geom.pDy2 = 1.254731;
      geom.pDx3 = 1.339476;
      geom.pDx4 = 1.319820;
      geom.centralX = 1.265973;
      geom.centralY = 104.107425;
      geom.centralZ = -87.947881;
      geom.pRotationAngleX = -2.263258;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 676 based Row/Col = 67/6
      geom_tower geom;
      geom.id = 676;
      geom.pDz = 7.050152;
      geom.pTheta = 0.029644;
      geom.pPhi = 0.317994;
      geom.pAlp1 = -0.023518;
      geom.pAlp2 = -0.023512;
      geom.pDy1 = 1.125572;
      geom.pDx1 = 1.206981;
      geom.pDx2 = 1.189305;
      geom.pDy2 = 1.254731;
      geom.pDx3 = 1.340437;
      geom.pDx4 = 1.320738;
      geom.centralX = 3.798810;
      geom.centralY = 104.107425;
      geom.centralZ = -87.947881;
      geom.pRotationAngleX = -2.263258;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 677 based Row/Col = 67/7
      geom_tower geom;
      geom.id = 677;
      geom.pDz = 7.050152;
      geom.pTheta = 0.047835;
      geom.pPhi = 0.194900;
      geom.pAlp1 = -0.039240;
      geom.pAlp2 = -0.039230;
      geom.pDy1 = 1.125572;
      geom.pDx1 = 1.208714;
      geom.pDx2 = 1.190961;
      geom.pDy2 = 1.254731;
      geom.pDx3 = 1.342361;
      geom.pDx4 = 1.322577;
      geom.centralX = 6.334329;
      geom.centralY = 104.107425;
      geom.centralZ = -87.947881;
      geom.pRotationAngleX = -2.263258;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 678 based Row/Col = 67/8
      geom_tower geom;
      geom.id = 678;
      geom.pDz = 7.050152;
      geom.pTheta = 0.066350;
      geom.pPhi = 0.139983;
      geom.pAlp1 = -0.055028;
      geom.pAlp2 = -0.055013;
      geom.pDy1 = 1.125572;
      geom.pDx1 = 1.211320;
      geom.pDx2 = 1.193452;
      geom.pDy2 = 1.254731;
      geom.pDx3 = 1.345254;
      geom.pDx4 = 1.325341;
      geom.centralX = 8.874324;
      geom.centralY = 104.107425;
      geom.centralZ = -87.947881;
      geom.pRotationAngleX = -2.263258;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 681 based Row/Col = 68/1
      geom_tower geom;
      geom.id = 681;
      geom.pDz = 7.050105;
      geom.pTheta = 0.067313;
      geom.pPhi = -3.006958;
      geom.pAlp1 = 0.055027;
      geom.pAlp2 = 0.055012;
      geom.pDy1 = 1.125572;
      geom.pDx1 = 1.229192;
      geom.pDx2 = 1.211320;
      geom.pDy2 = 1.254731;
      geom.pDx3 = 1.365172;
      geom.pDx4 = 1.345254;
      geom.centralX = -9.005420;
      geom.centralY = 105.630282;
      geom.centralZ = -86.111969;
      geom.pRotationAngleX = -2.263258;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 682 based Row/Col = 68/2
      geom_tower geom;
      geom.id = 682;
      geom.pDz = 7.050105;
      geom.pTheta = 0.048496;
      geom.pPhi = -2.954049;
      geom.pAlp1 = 0.039236;
      geom.pAlp2 = 0.039226;
      geom.pDy1 = 1.125572;
      geom.pDx1 = 1.226468;
      geom.pDx2 = 1.208714;
      geom.pDy2 = 1.254731;
      geom.pDx3 = 1.362146;
      geom.pDx4 = 1.342361;
      geom.centralX = -6.427763;
      geom.centralY = 105.630282;
      geom.centralZ = -86.111969;
      geom.pRotationAngleX = -2.263258;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 683 based Row/Col = 68/3
      geom_tower geom;
      geom.id = 683;
      geom.pDz = 7.050105;
      geom.pTheta = 0.029984;
      geom.pPhi = -2.835119;
      geom.pAlp1 = 0.023515;
      geom.pAlp2 = 0.023508;
      geom.pDy1 = 1.125572;
      geom.pDx1 = 1.224656;
      geom.pDx2 = 1.206981;
      geom.pDy2 = 1.254731;
      geom.pDx3 = 1.360134;
      geom.pDx4 = 1.340437;
      geom.centralX = -3.854789;
      geom.centralY = 105.630282;
      geom.centralZ = -86.111969;
      geom.pRotationAngleX = -2.263258;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 684 based Row/Col = 68/4
      geom_tower geom;
      geom.id = 684;
      geom.pDz = 7.050105;
      geom.pTheta = 0.013141;
      geom.pPhi = -2.382061;
      geom.pAlp1 = 0.007834;
      geom.pAlp2 = 0.007832;
      geom.pDy1 = 1.125572;
      geom.pDx1 = 1.223751;
      geom.pDx2 = 1.206116;
      geom.pDy2 = 1.254731;
      geom.pDx3 = 1.359130;
      geom.pDx4 = 1.339476;
      geom.centralX = -1.284618;
      geom.centralY = 105.630282;
      geom.centralZ = -86.111969;
      geom.pRotationAngleX = -2.263258;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 685 based Row/Col = 68/5
      geom_tower geom;
      geom.id = 685;
      geom.pDz = 7.050105;
      geom.pTheta = 0.013141;
      geom.pPhi = -0.759532;
      geom.pAlp1 = -0.007834;
      geom.pAlp2 = -0.007832;
      geom.pDy1 = 1.125572;
      geom.pDx1 = 1.223751;
      geom.pDx2 = 1.206116;
      geom.pDy2 = 1.254731;
      geom.pDx3 = 1.359130;
      geom.pDx4 = 1.339476;
      geom.centralX = 1.284618;
      geom.centralY = 105.630282;
      geom.centralZ = -86.111969;
      geom.pRotationAngleX = -2.263258;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 686 based Row/Col = 68/6
      geom_tower geom;
      geom.id = 686;
      geom.pDz = 7.050105;
      geom.pTheta = 0.029984;
      geom.pPhi = -0.306474;
      geom.pAlp1 = -0.023515;
      geom.pAlp2 = -0.023508;
      geom.pDy1 = 1.125572;
      geom.pDx1 = 1.224656;
      geom.pDx2 = 1.206981;
      geom.pDy2 = 1.254731;
      geom.pDx3 = 1.360134;
      geom.pDx4 = 1.340437;
      geom.centralX = 3.854789;
      geom.centralY = 105.630282;
      geom.centralZ = -86.111969;
      geom.pRotationAngleX = -2.263258;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 687 based Row/Col = 68/7
      geom_tower geom;
      geom.id = 687;
      geom.pDz = 7.050105;
      geom.pTheta = 0.048496;
      geom.pPhi = -0.187544;
      geom.pAlp1 = -0.039236;
      geom.pAlp2 = -0.039226;
      geom.pDy1 = 1.125572;
      geom.pDx1 = 1.226468;
      geom.pDx2 = 1.208714;
      geom.pDy2 = 1.254731;
      geom.pDx3 = 1.362146;
      geom.pDx4 = 1.342361;
      geom.centralX = 6.427763;
      geom.centralY = 105.630282;
      geom.centralZ = -86.111969;
      geom.pRotationAngleX = -2.263258;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 688 based Row/Col = 68/8
      geom_tower geom;
      geom.id = 688;
      geom.pDz = 7.050105;
      geom.pTheta = 0.067313;
      geom.pPhi = -0.134635;
      geom.pAlp1 = -0.055027;
      geom.pAlp2 = -0.055012;
      geom.pDy1 = 1.125572;
      geom.pDx1 = 1.229192;
      geom.pDx2 = 1.211320;
      geom.pDy2 = 1.254731;
      geom.pDx3 = 1.365172;
      geom.pDx4 = 1.345254;
      geom.centralX = 9.005420;
      geom.centralY = 105.630282;
      geom.centralZ = -86.111969;
      geom.pRotationAngleX = -2.263258;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1381 based Row/Col = 138/1
      geom_tower geom;
      geom.id = 1381;
      geom.pDz = 7.329474;
      geom.pTheta = 0.062437;
      geom.pPhi = -3.002302;
      geom.pAlp1 = -0.059395;
      geom.pAlp2 = -0.059385;
      geom.pDy1 = 1.123036;
      geom.pDx1 = 1.193132;
      geom.pDx2 = 1.212359;
      geom.pDy2 = 1.248304;
      geom.pDx3 = 1.322039;
      geom.pDx4 = 1.343406;
      geom.centralX = -8.869453;
      geom.centralY = 104.046719;
      geom.centralZ = 100.627910;
      geom.pRotationAngleX = -0.810475;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1382 based Row/Col = 138/2
      geom_tower geom;
      geom.id = 1382;
      geom.pDz = 7.329474;
      geom.pTheta = 0.045009;
      geom.pPhi = -2.947661;
      geom.pAlp1 = -0.042369;
      geom.pAlp2 = -0.042361;
      geom.pDy1 = 1.123036;
      geom.pDx1 = 1.190930;
      geom.pDx2 = 1.210047;
      geom.pDy2 = 1.248304;
      geom.pDx3 = 1.319599;
      geom.pDx4 = 1.340845;
      geom.centralX = -6.331364;
      geom.centralY = 104.046719;
      geom.centralZ = 100.627910;
      geom.pRotationAngleX = -0.810475;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1383 based Row/Col = 138/3
      geom_tower geom;
      geom.id = 1383;
      geom.pDz = 7.329474;
      geom.pTheta = 0.027884;
      geom.pPhi = -2.825128;
      geom.pAlp1 = -0.025399;
      geom.pAlp2 = -0.025394;
      geom.pDy1 = 1.123036;
      geom.pDx1 = 1.189464;
      geom.pDx2 = 1.208508;
      geom.pDy2 = 1.248304;
      geom.pDx3 = 1.317976;
      geom.pDx4 = 1.339140;
      geom.centralX = -3.797236;
      geom.centralY = 104.046719;
      geom.centralZ = 100.627910;
      geom.pRotationAngleX = -0.810475;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1384 based Row/Col = 138/4
      geom_tower geom;
      geom.id = 1384;
      geom.pDz = 7.329474;
      geom.pTheta = 0.012384;
      geom.pPhi = -2.364962;
      geom.pAlp1 = -0.008462;
      geom.pAlp2 = -0.008461;
      geom.pDy1 = 1.123036;
      geom.pDx1 = 1.188733;
      geom.pDx2 = 1.207740;
      geom.pDy2 = 1.248304;
      geom.pDx3 = 1.317165;
      geom.pDx4 = 1.338289;
      geom.centralX = -1.265482;
      geom.centralY = 104.046719;
      geom.centralZ = 100.627910;
      geom.pRotationAngleX = -0.810475;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1385 based Row/Col = 138/5
      geom_tower geom;
      geom.id = 1385;
      geom.pDz = 7.329474;
      geom.pTheta = 0.012384;
      geom.pPhi = -0.776630;
      geom.pAlp1 = 0.008462;
      geom.pAlp2 = 0.008461;
      geom.pDy1 = 1.123036;
      geom.pDx1 = 1.188733;
      geom.pDx2 = 1.207740;
      geom.pDy2 = 1.248304;
      geom.pDx3 = 1.317165;
      geom.pDx4 = 1.338289;
      geom.centralX = 1.265482;
      geom.centralY = 104.046719;
      geom.centralZ = 100.627910;
      geom.pRotationAngleX = -0.810475;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1386 based Row/Col = 138/6
      geom_tower geom;
      geom.id = 1386;
      geom.pDz = 7.329474;
      geom.pTheta = 0.027884;
      geom.pPhi = -0.316465;
      geom.pAlp1 = 0.025399;
      geom.pAlp2 = 0.025394;
      geom.pDy1 = 1.123036;
      geom.pDx1 = 1.189464;
      geom.pDx2 = 1.208508;
      geom.pDy2 = 1.248304;
      geom.pDx3 = 1.317976;
      geom.pDx4 = 1.339140;
      geom.centralX = 3.797236;
      geom.centralY = 104.046719;
      geom.centralZ = 100.627910;
      geom.pRotationAngleX = -0.810475;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1387 based Row/Col = 138/7
      geom_tower geom;
      geom.id = 1387;
      geom.pDz = 7.329474;
      geom.pTheta = 0.045009;
      geom.pPhi = -0.193932;
      geom.pAlp1 = 0.042369;
      geom.pAlp2 = 0.042361;
      geom.pDy1 = 1.123036;
      geom.pDx1 = 1.190930;
      geom.pDx2 = 1.210047;
      geom.pDy2 = 1.248304;
      geom.pDx3 = 1.319599;
      geom.pDx4 = 1.340845;
      geom.centralX = 6.331364;
      geom.centralY = 104.046719;
      geom.centralZ = 100.627910;
      geom.pRotationAngleX = -0.810475;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1388 based Row/Col = 138/8
      geom_tower geom;
      geom.id = 1388;
      geom.pDz = 7.329474;
      geom.pTheta = 0.062437;
      geom.pPhi = -0.139290;
      geom.pAlp1 = 0.059395;
      geom.pAlp2 = 0.059385;
      geom.pDy1 = 1.123036;
      geom.pDx1 = 1.193132;
      geom.pDx2 = 1.212359;
      geom.pDy2 = 1.248304;
      geom.pDx3 = 1.322039;
      geom.pDx4 = 1.343406;
      geom.centralX = 8.869453;
      geom.centralY = 104.046719;
      geom.centralZ = 100.627910;
      geom.pRotationAngleX = -0.810475;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1371 based Row/Col = 137/1
      geom_tower geom;
      geom.id = 1371;
      geom.pDz = 7.329517;
      geom.pTheta = 0.063401;
      geom.pPhi = 3.008717;
      geom.pAlp1 = -0.059394;
      geom.pAlp2 = -0.059383;
      geom.pDy1 = 1.123036;
      geom.pDx1 = 1.212359;
      geom.pDx2 = 1.231590;
      geom.pDy2 = 1.248304;
      geom.pDx3 = 1.343406;
      geom.pDx4 = 1.364778;
      geom.centralX = -9.010450;
      geom.centralY = 105.684384;
      geom.centralZ = 98.905981;
      geom.pRotationAngleX = -0.810475;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1372 based Row/Col = 137/2
      geom_tower geom;
      geom.id = 1372;
      geom.pDz = 7.329517;
      geom.pTheta = 0.045667;
      geom.pPhi = 2.956487;
      geom.pAlp1 = -0.042365;
      geom.pAlp2 = -0.042357;
      geom.pDy1 = 1.123036;
      geom.pDx1 = 1.210047;
      geom.pDx2 = 1.229164;
      geom.pDy2 = 1.248304;
      geom.pDx3 = 1.340845;
      geom.pDx4 = 1.362091;
      geom.centralX = -6.431880;
      geom.centralY = 105.684384;
      geom.centralZ = 98.905981;
      geom.pRotationAngleX = -0.810475;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1373 based Row/Col = 137/3
      geom_tower geom;
      geom.id = 1373;
      geom.pDz = 7.329517;
      geom.pTheta = 0.028214;
      geom.pPhi = 2.838967;
      geom.pAlp1 = -0.025395;
      geom.pAlp2 = -0.025391;
      geom.pDy1 = 1.123036;
      geom.pDx1 = 1.208508;
      geom.pDx2 = 1.227551;
      geom.pDy2 = 1.248304;
      geom.pDx3 = 1.339140;
      geom.pDx4 = 1.360303;
      geom.centralX = -3.857468;
      geom.centralY = 105.684384;
      geom.centralZ = 98.905981;
      geom.pRotationAngleX = -0.810475;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1374 based Row/Col = 137/4
      geom_tower geom;
      geom.id = 1374;
      geom.pDz = 7.329517;
      geom.pTheta = 0.012302;
      geom.pPhi = 2.388788;
      geom.pAlp1 = -0.008461;
      geom.pAlp2 = -0.008460;
      geom.pDy1 = 1.123036;
      geom.pDx1 = 1.207740;
      geom.pDx2 = 1.226745;
      geom.pDy2 = 1.248304;
      geom.pDx3 = 1.338289;
      geom.pDx4 = 1.359410;
      geom.centralX = -1.285546;
      geom.centralY = 105.684384;
      geom.centralZ = 98.905981;
      geom.pRotationAngleX = -0.810475;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1375 based Row/Col = 137/5
      geom_tower geom;
      geom.id = 1375;
      geom.pDz = 7.329517;
      geom.pTheta = 0.012302;
      geom.pPhi = 0.752805;
      geom.pAlp1 = 0.008461;
      geom.pAlp2 = 0.008460;
      geom.pDy1 = 1.123036;
      geom.pDx1 = 1.207740;
      geom.pDx2 = 1.226745;
      geom.pDy2 = 1.248304;
      geom.pDx3 = 1.338289;
      geom.pDx4 = 1.359410;
      geom.centralX = 1.285546;
      geom.centralY = 105.684384;
      geom.centralZ = 98.905981;
      geom.pRotationAngleX = -0.810475;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1376 based Row/Col = 137/6
      geom_tower geom;
      geom.id = 1376;
      geom.pDz = 7.329517;
      geom.pTheta = 0.028214;
      geom.pPhi = 0.302626;
      geom.pAlp1 = 0.025395;
      geom.pAlp2 = 0.025391;
      geom.pDy1 = 1.123036;
      geom.pDx1 = 1.208508;
      geom.pDx2 = 1.227551;
      geom.pDy2 = 1.248304;
      geom.pDx3 = 1.339140;
      geom.pDx4 = 1.360303;
      geom.centralX = 3.857468;
      geom.centralY = 105.684384;
      geom.centralZ = 98.905981;
      geom.pRotationAngleX = -0.810475;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1377 based Row/Col = 137/7
      geom_tower geom;
      geom.id = 1377;
      geom.pDz = 7.329517;
      geom.pTheta = 0.045667;
      geom.pPhi = 0.185105;
      geom.pAlp1 = 0.042365;
      geom.pAlp2 = 0.042357;
      geom.pDy1 = 1.123036;
      geom.pDx1 = 1.210047;
      geom.pDx2 = 1.229164;
      geom.pDy2 = 1.248304;
      geom.pDx3 = 1.340845;
      geom.pDx4 = 1.362091;
      geom.centralX = 6.431880;
      geom.centralY = 105.684384;
      geom.centralZ = 98.905981;
      geom.pRotationAngleX = -0.810475;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1378 based Row/Col = 137/8
      geom_tower geom;
      geom.id = 1378;
      geom.pDz = 7.329517;
      geom.pTheta = 0.063401;
      geom.pPhi = 0.132875;
      geom.pAlp1 = 0.059394;
      geom.pAlp2 = 0.059383;
      geom.pDy1 = 1.123036;
      geom.pDx1 = 1.212359;
      geom.pDx2 = 1.231590;
      geom.pDy2 = 1.248304;
      geom.pDx3 = 1.343406;
      geom.pDx4 = 1.364778;
      geom.centralX = 9.010450;
      geom.centralY = 105.684384;
      geom.centralZ = 98.905981;
      geom.pRotationAngleX = -0.810475;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 631 based Row/Col = 63/1
      geom_tower geom;
      geom.id = 631;
      geom.pDz = 7.329474;
      geom.pTheta = 0.062437;
      geom.pPhi = 3.002302;
      geom.pAlp1 = 0.059395;
      geom.pAlp2 = 0.059385;
      geom.pDy1 = 1.123036;
      geom.pDx1 = 1.212359;
      geom.pDx2 = 1.193132;
      geom.pDy2 = 1.248304;
      geom.pDx3 = 1.343406;
      geom.pDx4 = 1.322039;
      geom.centralX = -8.869453;
      geom.centralY = 104.046719;
      geom.centralZ = -100.627910;
      geom.pRotationAngleX = -2.331118;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 632 based Row/Col = 63/2
      geom_tower geom;
      geom.id = 632;
      geom.pDz = 7.329474;
      geom.pTheta = 0.045009;
      geom.pPhi = 2.947661;
      geom.pAlp1 = 0.042369;
      geom.pAlp2 = 0.042361;
      geom.pDy1 = 1.123036;
      geom.pDx1 = 1.210047;
      geom.pDx2 = 1.190930;
      geom.pDy2 = 1.248304;
      geom.pDx3 = 1.340845;
      geom.pDx4 = 1.319599;
      geom.centralX = -6.331364;
      geom.centralY = 104.046719;
      geom.centralZ = -100.627910;
      geom.pRotationAngleX = -2.331118;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 633 based Row/Col = 63/3
      geom_tower geom;
      geom.id = 633;
      geom.pDz = 7.329474;
      geom.pTheta = 0.027884;
      geom.pPhi = 2.825128;
      geom.pAlp1 = 0.025399;
      geom.pAlp2 = 0.025394;
      geom.pDy1 = 1.123036;
      geom.pDx1 = 1.208508;
      geom.pDx2 = 1.189464;
      geom.pDy2 = 1.248304;
      geom.pDx3 = 1.339140;
      geom.pDx4 = 1.317976;
      geom.centralX = -3.797236;
      geom.centralY = 104.046719;
      geom.centralZ = -100.627910;
      geom.pRotationAngleX = -2.331118;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 634 based Row/Col = 63/4
      geom_tower geom;
      geom.id = 634;
      geom.pDz = 7.329474;
      geom.pTheta = 0.012384;
      geom.pPhi = 2.364962;
      geom.pAlp1 = 0.008462;
      geom.pAlp2 = 0.008461;
      geom.pDy1 = 1.123036;
      geom.pDx1 = 1.207740;
      geom.pDx2 = 1.188733;
      geom.pDy2 = 1.248304;
      geom.pDx3 = 1.338289;
      geom.pDx4 = 1.317165;
      geom.centralX = -1.265482;
      geom.centralY = 104.046719;
      geom.centralZ = -100.627910;
      geom.pRotationAngleX = -2.331118;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 635 based Row/Col = 63/5
      geom_tower geom;
      geom.id = 635;
      geom.pDz = 7.329474;
      geom.pTheta = 0.012384;
      geom.pPhi = 0.776630;
      geom.pAlp1 = -0.008462;
      geom.pAlp2 = -0.008461;
      geom.pDy1 = 1.123036;
      geom.pDx1 = 1.207740;
      geom.pDx2 = 1.188733;
      geom.pDy2 = 1.248304;
      geom.pDx3 = 1.338289;
      geom.pDx4 = 1.317165;
      geom.centralX = 1.265482;
      geom.centralY = 104.046719;
      geom.centralZ = -100.627910;
      geom.pRotationAngleX = -2.331118;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 636 based Row/Col = 63/6
      geom_tower geom;
      geom.id = 636;
      geom.pDz = 7.329474;
      geom.pTheta = 0.027884;
      geom.pPhi = 0.316465;
      geom.pAlp1 = -0.025399;
      geom.pAlp2 = -0.025394;
      geom.pDy1 = 1.123036;
      geom.pDx1 = 1.208508;
      geom.pDx2 = 1.189464;
      geom.pDy2 = 1.248304;
      geom.pDx3 = 1.339140;
      geom.pDx4 = 1.317976;
      geom.centralX = 3.797236;
      geom.centralY = 104.046719;
      geom.centralZ = -100.627910;
      geom.pRotationAngleX = -2.331118;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 637 based Row/Col = 63/7
      geom_tower geom;
      geom.id = 637;
      geom.pDz = 7.329474;
      geom.pTheta = 0.045009;
      geom.pPhi = 0.193932;
      geom.pAlp1 = -0.042369;
      geom.pAlp2 = -0.042361;
      geom.pDy1 = 1.123036;
      geom.pDx1 = 1.210047;
      geom.pDx2 = 1.190930;
      geom.pDy2 = 1.248304;
      geom.pDx3 = 1.340845;
      geom.pDx4 = 1.319599;
      geom.centralX = 6.331364;
      geom.centralY = 104.046719;
      geom.centralZ = -100.627910;
      geom.pRotationAngleX = -2.331118;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 638 based Row/Col = 63/8
      geom_tower geom;
      geom.id = 638;
      geom.pDz = 7.329474;
      geom.pTheta = 0.062437;
      geom.pPhi = 0.139290;
      geom.pAlp1 = -0.059395;
      geom.pAlp2 = -0.059385;
      geom.pDy1 = 1.123036;
      geom.pDx1 = 1.212359;
      geom.pDx2 = 1.193132;
      geom.pDy2 = 1.248304;
      geom.pDx3 = 1.343406;
      geom.pDx4 = 1.322039;
      geom.centralX = 8.869453;
      geom.centralY = 104.046719;
      geom.centralZ = -100.627910;
      geom.pRotationAngleX = -2.331118;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 641 based Row/Col = 64/1
      geom_tower geom;
      geom.id = 641;
      geom.pDz = 7.329517;
      geom.pTheta = 0.063401;
      geom.pPhi = -3.008717;
      geom.pAlp1 = 0.059394;
      geom.pAlp2 = 0.059383;
      geom.pDy1 = 1.123036;
      geom.pDx1 = 1.231590;
      geom.pDx2 = 1.212359;
      geom.pDy2 = 1.248304;
      geom.pDx3 = 1.364778;
      geom.pDx4 = 1.343406;
      geom.centralX = -9.010450;
      geom.centralY = 105.684384;
      geom.centralZ = -98.905981;
      geom.pRotationAngleX = -2.331118;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 642 based Row/Col = 64/2
      geom_tower geom;
      geom.id = 642;
      geom.pDz = 7.329517;
      geom.pTheta = 0.045667;
      geom.pPhi = -2.956487;
      geom.pAlp1 = 0.042365;
      geom.pAlp2 = 0.042357;
      geom.pDy1 = 1.123036;
      geom.pDx1 = 1.229164;
      geom.pDx2 = 1.210047;
      geom.pDy2 = 1.248304;
      geom.pDx3 = 1.362091;
      geom.pDx4 = 1.340845;
      geom.centralX = -6.431880;
      geom.centralY = 105.684384;
      geom.centralZ = -98.905981;
      geom.pRotationAngleX = -2.331118;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 643 based Row/Col = 64/3
      geom_tower geom;
      geom.id = 643;
      geom.pDz = 7.329517;
      geom.pTheta = 0.028214;
      geom.pPhi = -2.838967;
      geom.pAlp1 = 0.025395;
      geom.pAlp2 = 0.025391;
      geom.pDy1 = 1.123036;
      geom.pDx1 = 1.227551;
      geom.pDx2 = 1.208508;
      geom.pDy2 = 1.248304;
      geom.pDx3 = 1.360303;
      geom.pDx4 = 1.339140;
      geom.centralX = -3.857468;
      geom.centralY = 105.684384;
      geom.centralZ = -98.905981;
      geom.pRotationAngleX = -2.331118;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 644 based Row/Col = 64/4
      geom_tower geom;
      geom.id = 644;
      geom.pDz = 7.329517;
      geom.pTheta = 0.012302;
      geom.pPhi = -2.388788;
      geom.pAlp1 = 0.008461;
      geom.pAlp2 = 0.008460;
      geom.pDy1 = 1.123036;
      geom.pDx1 = 1.226745;
      geom.pDx2 = 1.207740;
      geom.pDy2 = 1.248304;
      geom.pDx3 = 1.359410;
      geom.pDx4 = 1.338289;
      geom.centralX = -1.285546;
      geom.centralY = 105.684384;
      geom.centralZ = -98.905981;
      geom.pRotationAngleX = -2.331118;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 645 based Row/Col = 64/5
      geom_tower geom;
      geom.id = 645;
      geom.pDz = 7.329517;
      geom.pTheta = 0.012302;
      geom.pPhi = -0.752805;
      geom.pAlp1 = -0.008461;
      geom.pAlp2 = -0.008460;
      geom.pDy1 = 1.123036;
      geom.pDx1 = 1.226745;
      geom.pDx2 = 1.207740;
      geom.pDy2 = 1.248304;
      geom.pDx3 = 1.359410;
      geom.pDx4 = 1.338289;
      geom.centralX = 1.285546;
      geom.centralY = 105.684384;
      geom.centralZ = -98.905981;
      geom.pRotationAngleX = -2.331118;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 646 based Row/Col = 64/6
      geom_tower geom;
      geom.id = 646;
      geom.pDz = 7.329517;
      geom.pTheta = 0.028214;
      geom.pPhi = -0.302626;
      geom.pAlp1 = -0.025395;
      geom.pAlp2 = -0.025391;
      geom.pDy1 = 1.123036;
      geom.pDx1 = 1.227551;
      geom.pDx2 = 1.208508;
      geom.pDy2 = 1.248304;
      geom.pDx3 = 1.360303;
      geom.pDx4 = 1.339140;
      geom.centralX = 3.857468;
      geom.centralY = 105.684384;
      geom.centralZ = -98.905981;
      geom.pRotationAngleX = -2.331118;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 647 based Row/Col = 64/7
      geom_tower geom;
      geom.id = 647;
      geom.pDz = 7.329517;
      geom.pTheta = 0.045667;
      geom.pPhi = -0.185105;
      geom.pAlp1 = -0.042365;
      geom.pAlp2 = -0.042357;
      geom.pDy1 = 1.123036;
      geom.pDx1 = 1.229164;
      geom.pDx2 = 1.210047;
      geom.pDy2 = 1.248304;
      geom.pDx3 = 1.362091;
      geom.pDx4 = 1.340845;
      geom.centralX = 6.431880;
      geom.centralY = 105.684384;
      geom.centralZ = -98.905981;
      geom.pRotationAngleX = -2.331118;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 648 based Row/Col = 64/8
      geom_tower geom;
      geom.id = 648;
      geom.pDz = 7.329517;
      geom.pTheta = 0.063401;
      geom.pPhi = -0.132875;
      geom.pAlp1 = -0.059394;
      geom.pAlp2 = -0.059383;
      geom.pDy1 = 1.123036;
      geom.pDx1 = 1.231590;
      geom.pDx2 = 1.212359;
      geom.pDy2 = 1.248304;
      geom.pDx3 = 1.364778;
      geom.pDx4 = 1.343406;
      geom.centralX = 9.010450;
      geom.centralY = 105.684384;
      geom.centralZ = -98.905981;
      geom.pRotationAngleX = -2.331118;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1421 based Row/Col = 142/1
      geom_tower geom;
      geom.id = 1421;
      geom.pDz = 7.685130;
      geom.pTheta = 0.058530;
      geom.pPhi = -3.003071;
      geom.pAlp1 = -0.063228;
      geom.pAlp2 = -0.063214;
      geom.pDy1 = 1.121258;
      geom.pDx1 = 1.192528;
      geom.pDx2 = 1.212946;
      geom.pDy2 = 1.243210;
      geom.pDx3 = 1.319133;
      geom.pDx4 = 1.341767;
      geom.centralX = -8.864168;
      geom.centralY = 103.981472;
      geom.centralZ = 114.073654;
      geom.pRotationAngleX = -0.747148;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1422 based Row/Col = 142/2
      geom_tower geom;
      geom.id = 1422;
      geom.pDz = 7.685130;
      geom.pTheta = 0.042189;
      geom.pPhi = -2.948733;
      geom.pAlp1 = -0.045117;
      geom.pAlp2 = -0.045107;
      geom.pDy1 = 1.121258;
      geom.pDx1 = 1.190595;
      geom.pDx2 = 1.210911;
      geom.pDy2 = 1.243210;
      geom.pDx3 = 1.316995;
      geom.pDx4 = 1.339516;
      geom.centralX = -6.328070;
      geom.centralY = 103.981472;
      geom.centralZ = 114.073654;
      geom.pRotationAngleX = -0.747148;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1423 based Row/Col = 142/3
      geom_tower geom;
      geom.id = 1423;
      geom.pDz = 7.685130;
      geom.pTheta = 0.026128;
      geom.pPhi = -2.826820;
      geom.pAlp1 = -0.027052;
      geom.pAlp2 = -0.027046;
      geom.pDy1 = 1.121258;
      geom.pDx1 = 1.189308;
      geom.pDx2 = 1.209557;
      geom.pDy2 = 1.243210;
      geom.pDx3 = 1.315572;
      geom.pDx4 = 1.338018;
      geom.centralX = -3.795452;
      geom.centralY = 103.981472;
      geom.centralZ = 114.073654;
      geom.pRotationAngleX = -0.747148;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1424 based Row/Col = 142/4
      geom_tower geom;
      geom.id = 1424;
      geom.pDz = 7.685130;
      geom.pTheta = 0.011578;
      geom.pPhi = -2.367842;
      geom.pAlp1 = -0.009014;
      geom.pAlp2 = -0.009012;
      geom.pDy1 = 1.121258;
      geom.pDx1 = 1.188666;
      geom.pDx2 = 1.208880;
      geom.pDy2 = 1.243210;
      geom.pDx3 = 1.314861;
      geom.pDx4 = 1.337270;
      geom.centralX = -1.264919;
      geom.centralY = 103.981472;
      geom.centralZ = 114.073654;
      geom.pRotationAngleX = -0.747148;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1425 based Row/Col = 142/5
      geom_tower geom;
      geom.id = 1425;
      geom.pDz = 7.685130;
      geom.pTheta = 0.011578;
      geom.pPhi = -0.773751;
      geom.pAlp1 = 0.009014;
      geom.pAlp2 = 0.009012;
      geom.pDy1 = 1.121258;
      geom.pDx1 = 1.188666;
      geom.pDx2 = 1.208880;
      geom.pDy2 = 1.243210;
      geom.pDx3 = 1.314861;
      geom.pDx4 = 1.337270;
      geom.centralX = 1.264919;
      geom.centralY = 103.981472;
      geom.centralZ = 114.073654;
      geom.pRotationAngleX = -0.747148;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1426 based Row/Col = 142/6
      geom_tower geom;
      geom.id = 1426;
      geom.pDz = 7.685130;
      geom.pTheta = 0.026128;
      geom.pPhi = -0.314772;
      geom.pAlp1 = 0.027052;
      geom.pAlp2 = 0.027046;
      geom.pDy1 = 1.121258;
      geom.pDx1 = 1.189308;
      geom.pDx2 = 1.209557;
      geom.pDy2 = 1.243210;
      geom.pDx3 = 1.315572;
      geom.pDx4 = 1.338018;
      geom.centralX = 3.795452;
      geom.centralY = 103.981472;
      geom.centralZ = 114.073654;
      geom.pRotationAngleX = -0.747148;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1427 based Row/Col = 142/7
      geom_tower geom;
      geom.id = 1427;
      geom.pDz = 7.685130;
      geom.pTheta = 0.042189;
      geom.pPhi = -0.192860;
      geom.pAlp1 = 0.045117;
      geom.pAlp2 = 0.045107;
      geom.pDy1 = 1.121258;
      geom.pDx1 = 1.190595;
      geom.pDx2 = 1.210911;
      geom.pDy2 = 1.243210;
      geom.pDx3 = 1.316995;
      geom.pDx4 = 1.339516;
      geom.centralX = 6.328070;
      geom.centralY = 103.981472;
      geom.centralZ = 114.073654;
      geom.pRotationAngleX = -0.747148;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1428 based Row/Col = 142/8
      geom_tower geom;
      geom.id = 1428;
      geom.pDz = 7.685130;
      geom.pTheta = 0.058530;
      geom.pPhi = -0.138521;
      geom.pAlp1 = 0.063228;
      geom.pAlp2 = 0.063214;
      geom.pDy1 = 1.121258;
      geom.pDx1 = 1.192528;
      geom.pDx2 = 1.212946;
      geom.pDy2 = 1.243210;
      geom.pDx3 = 1.319133;
      geom.pDx4 = 1.341767;
      geom.centralX = 8.864168;
      geom.centralY = 103.981472;
      geom.centralZ = 114.073654;
      geom.pRotationAngleX = -0.747148;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1411 based Row/Col = 141/1
      geom_tower geom;
      geom.id = 1411;
      geom.pDz = 7.685137;
      geom.pTheta = 0.059478;
      geom.pPhi = 3.010612;
      geom.pAlp1 = -0.063226;
      geom.pAlp2 = -0.063212;
      geom.pDy1 = 1.121258;
      geom.pDx1 = 1.212946;
      geom.pDx2 = 1.233369;
      geom.pDy2 = 1.243210;
      geom.pDx3 = 1.341767;
      geom.pDx4 = 1.364406;
      geom.centralX = -9.013848;
      geom.centralY = 105.719785;
      geom.centralZ = 112.463485;
      geom.pRotationAngleX = -0.747148;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1412 based Row/Col = 141/2
      geom_tower geom;
      geom.id = 1412;
      geom.pDz = 7.685137;
      geom.pTheta = 0.042832;
      geom.pPhi = 2.959111;
      geom.pAlp1 = -0.045113;
      geom.pAlp2 = -0.045103;
      geom.pDy1 = 1.121258;
      geom.pDx1 = 1.210911;
      geom.pDx2 = 1.231228;
      geom.pDy2 = 1.243210;
      geom.pDx3 = 1.339516;
      geom.pDx4 = 1.362038;
      geom.centralX = -6.434803;
      geom.centralY = 105.719785;
      geom.centralZ = 112.463485;
      geom.pRotationAngleX = -0.747148;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1413 based Row/Col = 141/3
      geom_tower geom;
      geom.id = 1413;
      geom.pDz = 7.685137;
      geom.pTheta = 0.026441;
      geom.pPhi = 2.843111;
      geom.pAlp1 = -0.027048;
      geom.pAlp2 = -0.027043;
      geom.pDy1 = 1.121258;
      geom.pDx1 = 1.209557;
      geom.pDx2 = 1.229804;
      geom.pDy2 = 1.243210;
      geom.pDx3 = 1.338018;
      geom.pDx4 = 1.360462;
      geom.centralX = -3.859419;
      geom.centralY = 105.719785;
      geom.centralZ = 112.463485;
      geom.pRotationAngleX = -0.747148;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1414 based Row/Col = 141/4
      geom_tower geom;
      geom.id = 1414;
      geom.pDz = 7.685137;
      geom.pTheta = 0.011465;
      geom.pPhi = 2.396109;
      geom.pAlp1 = -0.009013;
      geom.pAlp2 = -0.009011;
      geom.pDy1 = 1.121258;
      geom.pDx1 = 1.208880;
      geom.pDx2 = 1.229092;
      geom.pDy2 = 1.243210;
      geom.pDx3 = 1.337270;
      geom.pDx4 = 1.359676;
      geom.centralX = -1.286230;
      geom.centralY = 105.719785;
      geom.centralZ = 112.463485;
      geom.pRotationAngleX = -0.747148;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1415 based Row/Col = 141/5
      geom_tower geom;
      geom.id = 1415;
      geom.pDz = 7.685137;
      geom.pTheta = 0.011465;
      geom.pPhi = 0.745484;
      geom.pAlp1 = 0.009013;
      geom.pAlp2 = 0.009011;
      geom.pDy1 = 1.121258;
      geom.pDx1 = 1.208880;
      geom.pDx2 = 1.229092;
      geom.pDy2 = 1.243210;
      geom.pDx3 = 1.337270;
      geom.pDx4 = 1.359676;
      geom.centralX = 1.286230;
      geom.centralY = 105.719785;
      geom.centralZ = 112.463485;
      geom.pRotationAngleX = -0.747148;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1416 based Row/Col = 141/6
      geom_tower geom;
      geom.id = 1416;
      geom.pDz = 7.685137;
      geom.pTheta = 0.026441;
      geom.pPhi = 0.298482;
      geom.pAlp1 = 0.027048;
      geom.pAlp2 = 0.027043;
      geom.pDy1 = 1.121258;
      geom.pDx1 = 1.209557;
      geom.pDx2 = 1.229804;
      geom.pDy2 = 1.243210;
      geom.pDx3 = 1.338018;
      geom.pDx4 = 1.360462;
      geom.centralX = 3.859419;
      geom.centralY = 105.719785;
      geom.centralZ = 112.463485;
      geom.pRotationAngleX = -0.747148;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1417 based Row/Col = 141/7
      geom_tower geom;
      geom.id = 1417;
      geom.pDz = 7.685137;
      geom.pTheta = 0.042832;
      geom.pPhi = 0.182482;
      geom.pAlp1 = 0.045113;
      geom.pAlp2 = 0.045103;
      geom.pDy1 = 1.121258;
      geom.pDx1 = 1.210911;
      geom.pDx2 = 1.231228;
      geom.pDy2 = 1.243210;
      geom.pDx3 = 1.339516;
      geom.pDx4 = 1.362038;
      geom.centralX = 6.434803;
      geom.centralY = 105.719785;
      geom.centralZ = 112.463485;
      geom.pRotationAngleX = -0.747148;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1418 based Row/Col = 141/8
      geom_tower geom;
      geom.id = 1418;
      geom.pDz = 7.685137;
      geom.pTheta = 0.059478;
      geom.pPhi = 0.130981;
      geom.pAlp1 = 0.063226;
      geom.pAlp2 = 0.063212;
      geom.pDy1 = 1.121258;
      geom.pDx1 = 1.212946;
      geom.pDx2 = 1.233369;
      geom.pDy2 = 1.243210;
      geom.pDx3 = 1.341767;
      geom.pDx4 = 1.364406;
      geom.centralX = 9.013848;
      geom.centralY = 105.719785;
      geom.centralZ = 112.463485;
      geom.pRotationAngleX = -0.747148;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 591 based Row/Col = 59/1
      geom_tower geom;
      geom.id = 591;
      geom.pDz = 7.685130;
      geom.pTheta = 0.058530;
      geom.pPhi = 3.003071;
      geom.pAlp1 = 0.063228;
      geom.pAlp2 = 0.063214;
      geom.pDy1 = 1.121258;
      geom.pDx1 = 1.212946;
      geom.pDx2 = 1.192528;
      geom.pDy2 = 1.243210;
      geom.pDx3 = 1.341767;
      geom.pDx4 = 1.319133;
      geom.centralX = -8.864168;
      geom.centralY = 103.981472;
      geom.centralZ = -114.073654;
      geom.pRotationAngleX = -2.394445;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 592 based Row/Col = 59/2
      geom_tower geom;
      geom.id = 592;
      geom.pDz = 7.685130;
      geom.pTheta = 0.042189;
      geom.pPhi = 2.948733;
      geom.pAlp1 = 0.045117;
      geom.pAlp2 = 0.045107;
      geom.pDy1 = 1.121258;
      geom.pDx1 = 1.210911;
      geom.pDx2 = 1.190595;
      geom.pDy2 = 1.243210;
      geom.pDx3 = 1.339516;
      geom.pDx4 = 1.316995;
      geom.centralX = -6.328070;
      geom.centralY = 103.981472;
      geom.centralZ = -114.073654;
      geom.pRotationAngleX = -2.394445;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 593 based Row/Col = 59/3
      geom_tower geom;
      geom.id = 593;
      geom.pDz = 7.685130;
      geom.pTheta = 0.026128;
      geom.pPhi = 2.826820;
      geom.pAlp1 = 0.027052;
      geom.pAlp2 = 0.027046;
      geom.pDy1 = 1.121258;
      geom.pDx1 = 1.209557;
      geom.pDx2 = 1.189308;
      geom.pDy2 = 1.243210;
      geom.pDx3 = 1.338018;
      geom.pDx4 = 1.315572;
      geom.centralX = -3.795452;
      geom.centralY = 103.981472;
      geom.centralZ = -114.073654;
      geom.pRotationAngleX = -2.394445;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 594 based Row/Col = 59/4
      geom_tower geom;
      geom.id = 594;
      geom.pDz = 7.685130;
      geom.pTheta = 0.011578;
      geom.pPhi = 2.367842;
      geom.pAlp1 = 0.009014;
      geom.pAlp2 = 0.009012;
      geom.pDy1 = 1.121258;
      geom.pDx1 = 1.208880;
      geom.pDx2 = 1.188666;
      geom.pDy2 = 1.243210;
      geom.pDx3 = 1.337270;
      geom.pDx4 = 1.314861;
      geom.centralX = -1.264919;
      geom.centralY = 103.981472;
      geom.centralZ = -114.073654;
      geom.pRotationAngleX = -2.394445;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 595 based Row/Col = 59/5
      geom_tower geom;
      geom.id = 595;
      geom.pDz = 7.685130;
      geom.pTheta = 0.011578;
      geom.pPhi = 0.773751;
      geom.pAlp1 = -0.009014;
      geom.pAlp2 = -0.009012;
      geom.pDy1 = 1.121258;
      geom.pDx1 = 1.208880;
      geom.pDx2 = 1.188666;
      geom.pDy2 = 1.243210;
      geom.pDx3 = 1.337270;
      geom.pDx4 = 1.314861;
      geom.centralX = 1.264919;
      geom.centralY = 103.981472;
      geom.centralZ = -114.073654;
      geom.pRotationAngleX = -2.394445;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 596 based Row/Col = 59/6
      geom_tower geom;
      geom.id = 596;
      geom.pDz = 7.685130;
      geom.pTheta = 0.026128;
      geom.pPhi = 0.314772;
      geom.pAlp1 = -0.027052;
      geom.pAlp2 = -0.027046;
      geom.pDy1 = 1.121258;
      geom.pDx1 = 1.209557;
      geom.pDx2 = 1.189308;
      geom.pDy2 = 1.243210;
      geom.pDx3 = 1.338018;
      geom.pDx4 = 1.315572;
      geom.centralX = 3.795452;
      geom.centralY = 103.981472;
      geom.centralZ = -114.073654;
      geom.pRotationAngleX = -2.394445;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 597 based Row/Col = 59/7
      geom_tower geom;
      geom.id = 597;
      geom.pDz = 7.685130;
      geom.pTheta = 0.042189;
      geom.pPhi = 0.192860;
      geom.pAlp1 = -0.045117;
      geom.pAlp2 = -0.045107;
      geom.pDy1 = 1.121258;
      geom.pDx1 = 1.210911;
      geom.pDx2 = 1.190595;
      geom.pDy2 = 1.243210;
      geom.pDx3 = 1.339516;
      geom.pDx4 = 1.316995;
      geom.centralX = 6.328070;
      geom.centralY = 103.981472;
      geom.centralZ = -114.073654;
      geom.pRotationAngleX = -2.394445;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 598 based Row/Col = 59/8
      geom_tower geom;
      geom.id = 598;
      geom.pDz = 7.685130;
      geom.pTheta = 0.058530;
      geom.pPhi = 0.138521;
      geom.pAlp1 = -0.063228;
      geom.pAlp2 = -0.063214;
      geom.pDy1 = 1.121258;
      geom.pDx1 = 1.212946;
      geom.pDx2 = 1.192528;
      geom.pDy2 = 1.243210;
      geom.pDx3 = 1.341767;
      geom.pDx4 = 1.319133;
      geom.centralX = 8.864168;
      geom.centralY = 103.981472;
      geom.centralZ = -114.073654;
      geom.pRotationAngleX = -2.394445;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 601 based Row/Col = 60/1
      geom_tower geom;
      geom.id = 601;
      geom.pDz = 7.685137;
      geom.pTheta = 0.059478;
      geom.pPhi = -3.010612;
      geom.pAlp1 = 0.063226;
      geom.pAlp2 = 0.063212;
      geom.pDy1 = 1.121258;
      geom.pDx1 = 1.233369;
      geom.pDx2 = 1.212946;
      geom.pDy2 = 1.243210;
      geom.pDx3 = 1.364406;
      geom.pDx4 = 1.341767;
      geom.centralX = -9.013848;
      geom.centralY = 105.719785;
      geom.centralZ = -112.463485;
      geom.pRotationAngleX = -2.394445;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 602 based Row/Col = 60/2
      geom_tower geom;
      geom.id = 602;
      geom.pDz = 7.685137;
      geom.pTheta = 0.042832;
      geom.pPhi = -2.959111;
      geom.pAlp1 = 0.045113;
      geom.pAlp2 = 0.045103;
      geom.pDy1 = 1.121258;
      geom.pDx1 = 1.231228;
      geom.pDx2 = 1.210911;
      geom.pDy2 = 1.243210;
      geom.pDx3 = 1.362038;
      geom.pDx4 = 1.339516;
      geom.centralX = -6.434803;
      geom.centralY = 105.719785;
      geom.centralZ = -112.463485;
      geom.pRotationAngleX = -2.394445;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 603 based Row/Col = 60/3
      geom_tower geom;
      geom.id = 603;
      geom.pDz = 7.685137;
      geom.pTheta = 0.026441;
      geom.pPhi = -2.843111;
      geom.pAlp1 = 0.027048;
      geom.pAlp2 = 0.027043;
      geom.pDy1 = 1.121258;
      geom.pDx1 = 1.229804;
      geom.pDx2 = 1.209557;
      geom.pDy2 = 1.243210;
      geom.pDx3 = 1.360462;
      geom.pDx4 = 1.338018;
      geom.centralX = -3.859419;
      geom.centralY = 105.719785;
      geom.centralZ = -112.463485;
      geom.pRotationAngleX = -2.394445;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 604 based Row/Col = 60/4
      geom_tower geom;
      geom.id = 604;
      geom.pDz = 7.685137;
      geom.pTheta = 0.011465;
      geom.pPhi = -2.396109;
      geom.pAlp1 = 0.009013;
      geom.pAlp2 = 0.009011;
      geom.pDy1 = 1.121258;
      geom.pDx1 = 1.229092;
      geom.pDx2 = 1.208880;
      geom.pDy2 = 1.243210;
      geom.pDx3 = 1.359676;
      geom.pDx4 = 1.337270;
      geom.centralX = -1.286230;
      geom.centralY = 105.719785;
      geom.centralZ = -112.463485;
      geom.pRotationAngleX = -2.394445;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 605 based Row/Col = 60/5
      geom_tower geom;
      geom.id = 605;
      geom.pDz = 7.685137;
      geom.pTheta = 0.011465;
      geom.pPhi = -0.745484;
      geom.pAlp1 = -0.009013;
      geom.pAlp2 = -0.009011;
      geom.pDy1 = 1.121258;
      geom.pDx1 = 1.229092;
      geom.pDx2 = 1.208880;
      geom.pDy2 = 1.243210;
      geom.pDx3 = 1.359676;
      geom.pDx4 = 1.337270;
      geom.centralX = 1.286230;
      geom.centralY = 105.719785;
      geom.centralZ = -112.463485;
      geom.pRotationAngleX = -2.394445;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 606 based Row/Col = 60/6
      geom_tower geom;
      geom.id = 606;
      geom.pDz = 7.685137;
      geom.pTheta = 0.026441;
      geom.pPhi = -0.298482;
      geom.pAlp1 = -0.027048;
      geom.pAlp2 = -0.027043;
      geom.pDy1 = 1.121258;
      geom.pDx1 = 1.229804;
      geom.pDx2 = 1.209557;
      geom.pDy2 = 1.243210;
      geom.pDx3 = 1.360462;
      geom.pDx4 = 1.338018;
      geom.centralX = 3.859419;
      geom.centralY = 105.719785;
      geom.centralZ = -112.463485;
      geom.pRotationAngleX = -2.394445;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 607 based Row/Col = 60/7
      geom_tower geom;
      geom.id = 607;
      geom.pDz = 7.685137;
      geom.pTheta = 0.042832;
      geom.pPhi = -0.182482;
      geom.pAlp1 = -0.045113;
      geom.pAlp2 = -0.045103;
      geom.pDy1 = 1.121258;
      geom.pDx1 = 1.231228;
      geom.pDx2 = 1.210911;
      geom.pDy2 = 1.243210;
      geom.pDx3 = 1.362038;
      geom.pDx4 = 1.339516;
      geom.centralX = 6.434803;
      geom.centralY = 105.719785;
      geom.centralZ = -112.463485;
      geom.pRotationAngleX = -2.394445;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 608 based Row/Col = 60/8
      geom_tower geom;
      geom.id = 608;
      geom.pDz = 7.685137;
      geom.pTheta = 0.059478;
      geom.pPhi = -0.130981;
      geom.pAlp1 = -0.063226;
      geom.pAlp2 = -0.063212;
      geom.pDy1 = 1.121258;
      geom.pDx1 = 1.233369;
      geom.pDx2 = 1.212946;
      geom.pDy2 = 1.243210;
      geom.pDx3 = 1.364406;
      geom.pDx4 = 1.341767;
      geom.centralX = 9.013848;
      geom.centralY = 105.719785;
      geom.centralZ = -112.463485;
      geom.pRotationAngleX = -2.394445;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1461 based Row/Col = 146/1
      geom_tower geom;
      geom.id = 1461;
      geom.pDz = 8.104373;
      geom.pTheta = 0.054667;
      geom.pPhi = -3.000388;
      geom.pAlp1 = -0.066593;
      geom.pAlp2 = -0.066581;
      geom.pDy1 = 1.116523;
      geom.pDx1 = 1.192076;
      geom.pDx2 = 1.213477;
      geom.pDy2 = 1.238115;
      geom.pDx3 = 1.316602;
      geom.pDx4 = 1.340328;
      geom.centralX = -8.859871;
      geom.centralY = 103.927974;
      geom.centralZ = 128.392841;
      geom.pRotationAngleX = -0.687682;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1462 based Row/Col = 146/2
      geom_tower geom;
      geom.id = 1462;
      geom.pDz = 8.104373;
      geom.pTheta = 0.039418;
      geom.pPhi = -2.945056;
      geom.pAlp1 = -0.047531;
      geom.pAlp2 = -0.047522;
      geom.pDy1 = 1.116523;
      geom.pDx1 = 1.190394;
      geom.pDx2 = 1.211700;
      geom.pDy2 = 1.238115;
      geom.pDx3 = 1.314744;
      geom.pDx4 = 1.338366;
      geom.centralX = -6.325449;
      geom.centralY = 103.927974;
      geom.centralZ = 128.392841;
      geom.pRotationAngleX = -0.687682;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1463 based Row/Col = 146/3
      geom_tower geom;
      geom.id = 1463;
      geom.pDz = 8.104373;
      geom.pTheta = 0.024440;
      geom.pPhi = -2.821088;
      geom.pAlp1 = -0.028504;
      geom.pAlp2 = -0.028499;
      geom.pDy1 = 1.116523;
      geom.pDx1 = 1.189274;
      geom.pDx2 = 1.210518;
      geom.pDy2 = 1.238115;
      geom.pDx3 = 1.313507;
      geom.pDx4 = 1.337060;
      geom.centralX = -3.794059;
      geom.centralY = 103.927974;
      geom.centralZ = 128.392841;
      geom.pRotationAngleX = -0.687682;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1464 based Row/Col = 146/4
      geom_tower geom;
      geom.id = 1464;
      geom.pDz = 8.104373;
      geom.pTheta = 0.010913;
      geom.pPhi = -2.358194;
      geom.pAlp1 = -0.009499;
      geom.pAlp2 = -0.009497;
      geom.pDy1 = 1.116523;
      geom.pDx1 = 1.188714;
      geom.pDx2 = 1.209927;
      geom.pDy2 = 1.238115;
      geom.pDx3 = 1.312889;
      geom.pDx4 = 1.336408;
      geom.centralX = -1.264485;
      geom.centralY = 103.927974;
      geom.centralZ = 128.392841;
      geom.pRotationAngleX = -0.687682;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1465 based Row/Col = 146/5
      geom_tower geom;
      geom.id = 1465;
      geom.pDz = 8.104373;
      geom.pTheta = 0.010913;
      geom.pPhi = -0.783399;
      geom.pAlp1 = 0.009499;
      geom.pAlp2 = 0.009497;
      geom.pDy1 = 1.116523;
      geom.pDx1 = 1.188714;
      geom.pDx2 = 1.209927;
      geom.pDy2 = 1.238115;
      geom.pDx3 = 1.312889;
      geom.pDx4 = 1.336408;
      geom.centralX = 1.264485;
      geom.centralY = 103.927974;
      geom.centralZ = 128.392841;
      geom.pRotationAngleX = -0.687682;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1466 based Row/Col = 146/6
      geom_tower geom;
      geom.id = 1466;
      geom.pDz = 8.104373;
      geom.pTheta = 0.024440;
      geom.pPhi = -0.320505;
      geom.pAlp1 = 0.028504;
      geom.pAlp2 = 0.028499;
      geom.pDy1 = 1.116523;
      geom.pDx1 = 1.189274;
      geom.pDx2 = 1.210518;
      geom.pDy2 = 1.238115;
      geom.pDx3 = 1.313507;
      geom.pDx4 = 1.337060;
      geom.centralX = 3.794059;
      geom.centralY = 103.927974;
      geom.centralZ = 128.392841;
      geom.pRotationAngleX = -0.687682;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1467 based Row/Col = 146/7
      geom_tower geom;
      geom.id = 1467;
      geom.pDz = 8.104373;
      geom.pTheta = 0.039418;
      geom.pPhi = -0.196536;
      geom.pAlp1 = 0.047531;
      geom.pAlp2 = 0.047522;
      geom.pDy1 = 1.116523;
      geom.pDx1 = 1.190394;
      geom.pDx2 = 1.211700;
      geom.pDy2 = 1.238115;
      geom.pDx3 = 1.314744;
      geom.pDx4 = 1.338366;
      geom.centralX = 6.325449;
      geom.centralY = 103.927974;
      geom.centralZ = 128.392841;
      geom.pRotationAngleX = -0.687682;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1468 based Row/Col = 146/8
      geom_tower geom;
      geom.id = 1468;
      geom.pDz = 8.104373;
      geom.pTheta = 0.054667;
      geom.pPhi = -0.141205;
      geom.pAlp1 = 0.066593;
      geom.pAlp2 = 0.066581;
      geom.pDy1 = 1.116523;
      geom.pDx1 = 1.192076;
      geom.pDx2 = 1.213477;
      geom.pDy2 = 1.238115;
      geom.pDx3 = 1.316602;
      geom.pDx4 = 1.340328;
      geom.centralX = 8.859871;
      geom.centralY = 103.927974;
      geom.centralZ = 128.392841;
      geom.pRotationAngleX = -0.687682;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1451 based Row/Col = 145/1
      geom_tower geom;
      geom.id = 1451;
      geom.pDz = 8.104410;
      geom.pTheta = 0.055599;
      geom.pPhi = 3.010014;
      geom.pAlp1 = -0.066592;
      geom.pAlp2 = -0.066579;
      geom.pDy1 = 1.116523;
      geom.pDx1 = 1.213477;
      geom.pDx2 = 1.234882;
      geom.pDy2 = 1.238115;
      geom.pDx3 = 1.340328;
      geom.pDx4 = 1.364059;
      geom.centralX = -9.016888;
      geom.centralY = 105.751313;
      geom.centralZ = 126.895066;
      geom.pRotationAngleX = -0.687682;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1452 based Row/Col = 145/2
      geom_tower geom;
      geom.id = 1452;
      geom.pDz = 8.104410;
      geom.pTheta = 0.040041;
      geom.pPhi = 2.958301;
      geom.pAlp1 = -0.047527;
      geom.pAlp2 = -0.047518;
      geom.pDy1 = 1.116523;
      geom.pDx1 = 1.211700;
      geom.pDx2 = 1.233007;
      geom.pDy2 = 1.238115;
      geom.pDx3 = 1.338366;
      geom.pDx4 = 1.361989;
      geom.centralX = -6.437436;
      geom.centralY = 105.751313;
      geom.centralZ = 126.895066;
      geom.pRotationAngleX = -0.687682;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1453 based Row/Col = 145/3
      geom_tower geom;
      geom.id = 1453;
      geom.pDz = 8.104410;
      geom.pTheta = 0.024724;
      geom.pPhi = 2.841849;
      geom.pAlp1 = -0.028501;
      geom.pAlp2 = -0.028496;
      geom.pDy1 = 1.116523;
      geom.pDx1 = 1.210518;
      geom.pDx2 = 1.231760;
      geom.pDy2 = 1.238115;
      geom.pDx3 = 1.337060;
      geom.pDx4 = 1.360611;
      geom.centralX = -3.861184;
      geom.centralY = 105.751313;
      geom.centralZ = 126.895066;
      geom.pRotationAngleX = -0.687682;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1454 based Row/Col = 145/4
      geom_tower geom;
      geom.id = 1454;
      geom.pDz = 8.104410;
      geom.pTheta = 0.010739;
      geom.pPhi = 2.393887;
      geom.pAlp1 = -0.009498;
      geom.pAlp2 = -0.009496;
      geom.pDy1 = 1.116523;
      geom.pDx1 = 1.209927;
      geom.pDx2 = 1.231136;
      geom.pDy2 = 1.238115;
      geom.pDx3 = 1.336408;
      geom.pDx4 = 1.359923;
      geom.centralX = -1.286848;
      geom.centralY = 105.751313;
      geom.centralZ = 126.895066;
      geom.pRotationAngleX = -0.687682;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1455 based Row/Col = 145/5
      geom_tower geom;
      geom.id = 1455;
      geom.pDz = 8.104410;
      geom.pTheta = 0.010739;
      geom.pPhi = 0.747706;
      geom.pAlp1 = 0.009498;
      geom.pAlp2 = 0.009496;
      geom.pDy1 = 1.116523;
      geom.pDx1 = 1.209927;
      geom.pDx2 = 1.231136;
      geom.pDy2 = 1.238115;
      geom.pDx3 = 1.336408;
      geom.pDx4 = 1.359923;
      geom.centralX = 1.286848;
      geom.centralY = 105.751313;
      geom.centralZ = 126.895066;
      geom.pRotationAngleX = -0.687682;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1456 based Row/Col = 145/6
      geom_tower geom;
      geom.id = 1456;
      geom.pDz = 8.104410;
      geom.pTheta = 0.024724;
      geom.pPhi = 0.299744;
      geom.pAlp1 = 0.028501;
      geom.pAlp2 = 0.028496;
      geom.pDy1 = 1.116523;
      geom.pDx1 = 1.210518;
      geom.pDx2 = 1.231760;
      geom.pDy2 = 1.238115;
      geom.pDx3 = 1.337060;
      geom.pDx4 = 1.360611;
      geom.centralX = 3.861184;
      geom.centralY = 105.751313;
      geom.centralZ = 126.895066;
      geom.pRotationAngleX = -0.687682;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1457 based Row/Col = 145/7
      geom_tower geom;
      geom.id = 1457;
      geom.pDz = 8.104410;
      geom.pTheta = 0.040041;
      geom.pPhi = 0.183291;
      geom.pAlp1 = 0.047527;
      geom.pAlp2 = 0.047518;
      geom.pDy1 = 1.116523;
      geom.pDx1 = 1.211700;
      geom.pDx2 = 1.233007;
      geom.pDy2 = 1.238115;
      geom.pDx3 = 1.338366;
      geom.pDx4 = 1.361989;
      geom.centralX = 6.437436;
      geom.centralY = 105.751313;
      geom.centralZ = 126.895066;
      geom.pRotationAngleX = -0.687682;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1458 based Row/Col = 145/8
      geom_tower geom;
      geom.id = 1458;
      geom.pDz = 8.104410;
      geom.pTheta = 0.055599;
      geom.pPhi = 0.131578;
      geom.pAlp1 = 0.066592;
      geom.pAlp2 = 0.066579;
      geom.pDy1 = 1.116523;
      geom.pDx1 = 1.213477;
      geom.pDx2 = 1.234882;
      geom.pDy2 = 1.238115;
      geom.pDx3 = 1.340328;
      geom.pDx4 = 1.364059;
      geom.centralX = 9.016888;
      geom.centralY = 105.751313;
      geom.centralZ = 126.895066;
      geom.pRotationAngleX = -0.687682;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 551 based Row/Col = 55/1
      geom_tower geom;
      geom.id = 551;
      geom.pDz = 8.104373;
      geom.pTheta = 0.054667;
      geom.pPhi = 3.000388;
      geom.pAlp1 = 0.066593;
      geom.pAlp2 = 0.066581;
      geom.pDy1 = 1.116523;
      geom.pDx1 = 1.213477;
      geom.pDx2 = 1.192076;
      geom.pDy2 = 1.238115;
      geom.pDx3 = 1.340328;
      geom.pDx4 = 1.316602;
      geom.centralX = -8.859871;
      geom.centralY = 103.927974;
      geom.centralZ = -128.392841;
      geom.pRotationAngleX = -2.453911;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 552 based Row/Col = 55/2
      geom_tower geom;
      geom.id = 552;
      geom.pDz = 8.104373;
      geom.pTheta = 0.039418;
      geom.pPhi = 2.945056;
      geom.pAlp1 = 0.047531;
      geom.pAlp2 = 0.047522;
      geom.pDy1 = 1.116523;
      geom.pDx1 = 1.211700;
      geom.pDx2 = 1.190394;
      geom.pDy2 = 1.238115;
      geom.pDx3 = 1.338366;
      geom.pDx4 = 1.314744;
      geom.centralX = -6.325449;
      geom.centralY = 103.927974;
      geom.centralZ = -128.392841;
      geom.pRotationAngleX = -2.453911;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 553 based Row/Col = 55/3
      geom_tower geom;
      geom.id = 553;
      geom.pDz = 8.104373;
      geom.pTheta = 0.024440;
      geom.pPhi = 2.821088;
      geom.pAlp1 = 0.028504;
      geom.pAlp2 = 0.028499;
      geom.pDy1 = 1.116523;
      geom.pDx1 = 1.210518;
      geom.pDx2 = 1.189274;
      geom.pDy2 = 1.238115;
      geom.pDx3 = 1.337060;
      geom.pDx4 = 1.313507;
      geom.centralX = -3.794059;
      geom.centralY = 103.927974;
      geom.centralZ = -128.392841;
      geom.pRotationAngleX = -2.453911;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 554 based Row/Col = 55/4
      geom_tower geom;
      geom.id = 554;
      geom.pDz = 8.104373;
      geom.pTheta = 0.010913;
      geom.pPhi = 2.358194;
      geom.pAlp1 = 0.009499;
      geom.pAlp2 = 0.009497;
      geom.pDy1 = 1.116523;
      geom.pDx1 = 1.209927;
      geom.pDx2 = 1.188714;
      geom.pDy2 = 1.238115;
      geom.pDx3 = 1.336408;
      geom.pDx4 = 1.312889;
      geom.centralX = -1.264485;
      geom.centralY = 103.927974;
      geom.centralZ = -128.392841;
      geom.pRotationAngleX = -2.453911;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 555 based Row/Col = 55/5
      geom_tower geom;
      geom.id = 555;
      geom.pDz = 8.104373;
      geom.pTheta = 0.010913;
      geom.pPhi = 0.783399;
      geom.pAlp1 = -0.009499;
      geom.pAlp2 = -0.009497;
      geom.pDy1 = 1.116523;
      geom.pDx1 = 1.209927;
      geom.pDx2 = 1.188714;
      geom.pDy2 = 1.238115;
      geom.pDx3 = 1.336408;
      geom.pDx4 = 1.312889;
      geom.centralX = 1.264485;
      geom.centralY = 103.927974;
      geom.centralZ = -128.392841;
      geom.pRotationAngleX = -2.453911;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 556 based Row/Col = 55/6
      geom_tower geom;
      geom.id = 556;
      geom.pDz = 8.104373;
      geom.pTheta = 0.024440;
      geom.pPhi = 0.320505;
      geom.pAlp1 = -0.028504;
      geom.pAlp2 = -0.028499;
      geom.pDy1 = 1.116523;
      geom.pDx1 = 1.210518;
      geom.pDx2 = 1.189274;
      geom.pDy2 = 1.238115;
      geom.pDx3 = 1.337060;
      geom.pDx4 = 1.313507;
      geom.centralX = 3.794059;
      geom.centralY = 103.927974;
      geom.centralZ = -128.392841;
      geom.pRotationAngleX = -2.453911;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 557 based Row/Col = 55/7
      geom_tower geom;
      geom.id = 557;
      geom.pDz = 8.104373;
      geom.pTheta = 0.039418;
      geom.pPhi = 0.196536;
      geom.pAlp1 = -0.047531;
      geom.pAlp2 = -0.047522;
      geom.pDy1 = 1.116523;
      geom.pDx1 = 1.211700;
      geom.pDx2 = 1.190394;
      geom.pDy2 = 1.238115;
      geom.pDx3 = 1.338366;
      geom.pDx4 = 1.314744;
      geom.centralX = 6.325449;
      geom.centralY = 103.927974;
      geom.centralZ = -128.392841;
      geom.pRotationAngleX = -2.453911;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 558 based Row/Col = 55/8
      geom_tower geom;
      geom.id = 558;
      geom.pDz = 8.104373;
      geom.pTheta = 0.054667;
      geom.pPhi = 0.141205;
      geom.pAlp1 = -0.066593;
      geom.pAlp2 = -0.066581;
      geom.pDy1 = 1.116523;
      geom.pDx1 = 1.213477;
      geom.pDx2 = 1.192076;
      geom.pDy2 = 1.238115;
      geom.pDx3 = 1.340328;
      geom.pDx4 = 1.316602;
      geom.centralX = 8.859871;
      geom.centralY = 103.927974;
      geom.centralZ = -128.392841;
      geom.pRotationAngleX = -2.453911;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 561 based Row/Col = 56/1
      geom_tower geom;
      geom.id = 561;
      geom.pDz = 8.104410;
      geom.pTheta = 0.055599;
      geom.pPhi = -3.010014;
      geom.pAlp1 = 0.066592;
      geom.pAlp2 = 0.066579;
      geom.pDy1 = 1.116523;
      geom.pDx1 = 1.234882;
      geom.pDx2 = 1.213477;
      geom.pDy2 = 1.238115;
      geom.pDx3 = 1.364059;
      geom.pDx4 = 1.340328;
      geom.centralX = -9.016888;
      geom.centralY = 105.751313;
      geom.centralZ = -126.895066;
      geom.pRotationAngleX = -2.453911;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 562 based Row/Col = 56/2
      geom_tower geom;
      geom.id = 562;
      geom.pDz = 8.104410;
      geom.pTheta = 0.040041;
      geom.pPhi = -2.958301;
      geom.pAlp1 = 0.047527;
      geom.pAlp2 = 0.047518;
      geom.pDy1 = 1.116523;
      geom.pDx1 = 1.233007;
      geom.pDx2 = 1.211700;
      geom.pDy2 = 1.238115;
      geom.pDx3 = 1.361989;
      geom.pDx4 = 1.338366;
      geom.centralX = -6.437436;
      geom.centralY = 105.751313;
      geom.centralZ = -126.895066;
      geom.pRotationAngleX = -2.453911;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 563 based Row/Col = 56/3
      geom_tower geom;
      geom.id = 563;
      geom.pDz = 8.104410;
      geom.pTheta = 0.024724;
      geom.pPhi = -2.841849;
      geom.pAlp1 = 0.028501;
      geom.pAlp2 = 0.028496;
      geom.pDy1 = 1.116523;
      geom.pDx1 = 1.231760;
      geom.pDx2 = 1.210518;
      geom.pDy2 = 1.238115;
      geom.pDx3 = 1.360611;
      geom.pDx4 = 1.337060;
      geom.centralX = -3.861184;
      geom.centralY = 105.751313;
      geom.centralZ = -126.895066;
      geom.pRotationAngleX = -2.453911;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 564 based Row/Col = 56/4
      geom_tower geom;
      geom.id = 564;
      geom.pDz = 8.104410;
      geom.pTheta = 0.010739;
      geom.pPhi = -2.393887;
      geom.pAlp1 = 0.009498;
      geom.pAlp2 = 0.009496;
      geom.pDy1 = 1.116523;
      geom.pDx1 = 1.231136;
      geom.pDx2 = 1.209927;
      geom.pDy2 = 1.238115;
      geom.pDx3 = 1.359923;
      geom.pDx4 = 1.336408;
      geom.centralX = -1.286848;
      geom.centralY = 105.751313;
      geom.centralZ = -126.895066;
      geom.pRotationAngleX = -2.453911;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 565 based Row/Col = 56/5
      geom_tower geom;
      geom.id = 565;
      geom.pDz = 8.104410;
      geom.pTheta = 0.010739;
      geom.pPhi = -0.747706;
      geom.pAlp1 = -0.009498;
      geom.pAlp2 = -0.009496;
      geom.pDy1 = 1.116523;
      geom.pDx1 = 1.231136;
      geom.pDx2 = 1.209927;
      geom.pDy2 = 1.238115;
      geom.pDx3 = 1.359923;
      geom.pDx4 = 1.336408;
      geom.centralX = 1.286848;
      geom.centralY = 105.751313;
      geom.centralZ = -126.895066;
      geom.pRotationAngleX = -2.453911;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 566 based Row/Col = 56/6
      geom_tower geom;
      geom.id = 566;
      geom.pDz = 8.104410;
      geom.pTheta = 0.024724;
      geom.pPhi = -0.299744;
      geom.pAlp1 = -0.028501;
      geom.pAlp2 = -0.028496;
      geom.pDy1 = 1.116523;
      geom.pDx1 = 1.231760;
      geom.pDx2 = 1.210518;
      geom.pDy2 = 1.238115;
      geom.pDx3 = 1.360611;
      geom.pDx4 = 1.337060;
      geom.centralX = 3.861184;
      geom.centralY = 105.751313;
      geom.centralZ = -126.895066;
      geom.pRotationAngleX = -2.453911;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 567 based Row/Col = 56/7
      geom_tower geom;
      geom.id = 567;
      geom.pDz = 8.104410;
      geom.pTheta = 0.040041;
      geom.pPhi = -0.183291;
      geom.pAlp1 = -0.047527;
      geom.pAlp2 = -0.047518;
      geom.pDy1 = 1.116523;
      geom.pDx1 = 1.233007;
      geom.pDx2 = 1.211700;
      geom.pDy2 = 1.238115;
      geom.pDx3 = 1.361989;
      geom.pDx4 = 1.338366;
      geom.centralX = 6.437436;
      geom.centralY = 105.751313;
      geom.centralZ = -126.895066;
      geom.pRotationAngleX = -2.453911;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 568 based Row/Col = 56/8
      geom_tower geom;
      geom.id = 568;
      geom.pDz = 8.104410;
      geom.pTheta = 0.055599;
      geom.pPhi = -0.131578;
      geom.pAlp1 = -0.066592;
      geom.pAlp2 = -0.066579;
      geom.pDy1 = 1.116523;
      geom.pDx1 = 1.234882;
      geom.pDx2 = 1.213477;
      geom.pDy2 = 1.238115;
      geom.pDx3 = 1.364059;
      geom.pDx4 = 1.340328;
      geom.centralX = 9.016888;
      geom.centralY = 105.751313;
      geom.centralZ = -126.895066;
      geom.pRotationAngleX = -2.453911;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1041 based Row/Col = 104/1
      geom_tower geom;
      geom.id = 1041;
      geom.pDz = 6.751986;
      geom.pTheta = 0.086342;
      geom.pPhi = -3.014158;
      geom.pAlp1 = -0.005984;
      geom.pAlp2 = -0.005981;
      geom.pDy1 = 1.124739;
      geom.pDx1 = 1.196698;
      geom.pDx2 = 1.198647;
      geom.pDy2 = 1.273242;
      geom.pDx3 = 1.362991;
      geom.pDx4 = 1.365197;
      geom.centralX = -8.940723;
      geom.centralY = 104.903842;
      geom.centralZ = 8.491506;
      geom.pRotationAngleX = -1.501347;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1042 based Row/Col = 104/2
      geom_tower geom;
      geom.id = 1042;
      geom.pDz = 6.751986;
      geom.pTheta = 0.062151;
      geom.pPhi = -2.963891;
      geom.pAlp1 = -0.004259;
      geom.pAlp2 = -0.004257;
      geom.pDy1 = 1.124739;
      geom.pDx1 = 1.192396;
      geom.pDx2 = 1.194326;
      geom.pDy2 = 1.273242;
      geom.pDx3 = 1.358093;
      geom.pDx4 = 1.360277;
      geom.centralX = -6.378566;
      geom.centralY = 104.903842;
      geom.centralZ = 8.491506;
      geom.pRotationAngleX = -1.501347;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1043 based Row/Col = 104/3
      geom_tower geom;
      geom.id = 1043;
      geom.pDz = 6.751986;
      geom.pTheta = 0.038315;
      geom.pPhi = -2.850534;
      geom.pAlp1 = -0.002550;
      geom.pAlp2 = -0.002548;
      geom.pDy1 = 1.124739;
      geom.pDx1 = 1.189540;
      geom.pDx2 = 1.191456;
      geom.pDy2 = 1.273242;
      geom.pDx3 = 1.354841;
      geom.pDx4 = 1.357009;
      geom.centralX = -3.824082;
      geom.centralY = 104.903842;
      geom.centralZ = 8.491506;
      geom.pRotationAngleX = -1.501347;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1044 based Row/Col = 104/4
      geom_tower geom;
      geom.id = 1044;
      geom.pDz = 6.751986;
      geom.pTheta = 0.016452;
      geom.pPhi = -2.409299;
      geom.pAlp1 = -0.000849;
      geom.pAlp2 = -0.000848;
      geom.pDy1 = 1.124739;
      geom.pDx1 = 1.188116;
      geom.pDx2 = 1.190025;
      geom.pDy2 = 1.273242;
      geom.pDx3 = 1.353219;
      geom.pDx4 = 1.355380;
      geom.centralX = -1.274185;
      geom.centralY = 104.903842;
      geom.centralZ = 8.491506;
      geom.pRotationAngleX = -1.501347;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1045 based Row/Col = 104/5
      geom_tower geom;
      geom.id = 1045;
      geom.pDz = 6.751986;
      geom.pTheta = 0.016452;
      geom.pPhi = -0.732294;
      geom.pAlp1 = 0.000849;
      geom.pAlp2 = 0.000848;
      geom.pDy1 = 1.124739;
      geom.pDx1 = 1.188116;
      geom.pDx2 = 1.190025;
      geom.pDy2 = 1.273242;
      geom.pDx3 = 1.353219;
      geom.pDx4 = 1.355380;
      geom.centralX = 1.274185;
      geom.centralY = 104.903842;
      geom.centralZ = 8.491506;
      geom.pRotationAngleX = -1.501347;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1046 based Row/Col = 104/6
      geom_tower geom;
      geom.id = 1046;
      geom.pDz = 6.751986;
      geom.pTheta = 0.038315;
      geom.pPhi = -0.291058;
      geom.pAlp1 = 0.002550;
      geom.pAlp2 = 0.002548;
      geom.pDy1 = 1.124739;
      geom.pDx1 = 1.189540;
      geom.pDx2 = 1.191456;
      geom.pDy2 = 1.273242;
      geom.pDx3 = 1.354841;
      geom.pDx4 = 1.357009;
      geom.centralX = 3.824082;
      geom.centralY = 104.903842;
      geom.centralZ = 8.491506;
      geom.pRotationAngleX = -1.501347;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1047 based Row/Col = 104/7
      geom_tower geom;
      geom.id = 1047;
      geom.pDz = 6.751986;
      geom.pTheta = 0.062151;
      geom.pPhi = -0.177701;
      geom.pAlp1 = 0.004259;
      geom.pAlp2 = 0.004257;
      geom.pDy1 = 1.124739;
      geom.pDx1 = 1.192396;
      geom.pDx2 = 1.194326;
      geom.pDy2 = 1.273242;
      geom.pDx3 = 1.358093;
      geom.pDx4 = 1.360277;
      geom.centralX = 6.378566;
      geom.centralY = 104.903842;
      geom.centralZ = 8.491506;
      geom.pRotationAngleX = -1.501347;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1048 based Row/Col = 104/8
      geom_tower geom;
      geom.id = 1048;
      geom.pDz = 6.751986;
      geom.pTheta = 0.086342;
      geom.pPhi = -0.127434;
      geom.pAlp1 = 0.005984;
      geom.pAlp2 = 0.005981;
      geom.pDy1 = 1.124739;
      geom.pDx1 = 1.196698;
      geom.pDx2 = 1.198647;
      geom.pDy2 = 1.273242;
      geom.pDx3 = 1.362991;
      geom.pDx4 = 1.365197;
      geom.centralX = 8.940723;
      geom.centralY = 104.903842;
      geom.centralZ = 8.491506;
      geom.pRotationAngleX = -1.501347;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1031 based Row/Col = 103/1
      geom_tower geom;
      geom.id = 1031;
      geom.pDz = 6.751968;
      geom.pTheta = 0.086470;
      geom.pPhi = 3.014432;
      geom.pAlp1 = -0.005984;
      geom.pAlp2 = -0.005981;
      geom.pDy1 = 1.124739;
      geom.pDx1 = 1.198647;
      geom.pDx2 = 1.200597;
      geom.pDy2 = 1.273242;
      geom.pDx3 = 1.365197;
      geom.pDx4 = 1.367404;
      geom.centralX = -8.955069;
      geom.centralY = 105.070593;
      geom.centralZ = 6.094317;
      geom.pRotationAngleX = -1.501347;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1032 based Row/Col = 103/2
      geom_tower geom;
      geom.id = 1032;
      geom.pDz = 6.751968;
      geom.pTheta = 0.062242;
      geom.pPhi = 2.964268;
      geom.pAlp1 = -0.004259;
      geom.pAlp2 = -0.004257;
      geom.pDy1 = 1.124739;
      geom.pDx1 = 1.194326;
      geom.pDx2 = 1.196256;
      geom.pDy2 = 1.273242;
      geom.pDx3 = 1.360277;
      geom.pDx4 = 1.362460;
      geom.centralX = -6.388778;
      geom.centralY = 105.070593;
      geom.centralZ = 6.094317;
      geom.pRotationAngleX = -1.501347;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1033 based Row/Col = 103/3
      geom_tower geom;
      geom.id = 1033;
      geom.pDz = 6.751968;
      geom.pTheta = 0.038366;
      geom.pPhi = 2.851129;
      geom.pAlp1 = -0.002550;
      geom.pAlp2 = -0.002548;
      geom.pDy1 = 1.124739;
      geom.pDx1 = 1.191456;
      geom.pDx2 = 1.193373;
      geom.pDy2 = 1.273242;
      geom.pDx3 = 1.357009;
      geom.pDx4 = 1.359178;
      geom.centralX = -3.830194;
      geom.centralY = 105.070593;
      geom.centralZ = 6.094317;
      geom.pRotationAngleX = -1.501347;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1034 based Row/Col = 103/4
      geom_tower geom;
      geom.id = 1034;
      geom.pDz = 6.751968;
      geom.pTheta = 0.016461;
      geom.pPhi = 2.410375;
      geom.pAlp1 = -0.000849;
      geom.pAlp2 = -0.000848;
      geom.pDy1 = 1.124739;
      geom.pDx1 = 1.190025;
      geom.pDx2 = 1.191935;
      geom.pDy2 = 1.273242;
      geom.pDx3 = 1.355380;
      geom.pDx4 = 1.357540;
      geom.centralX = -1.276220;
      geom.centralY = 105.070593;
      geom.centralZ = 6.094317;
      geom.pRotationAngleX = -1.501347;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1035 based Row/Col = 103/5
      geom_tower geom;
      geom.id = 1035;
      geom.pDz = 6.751968;
      geom.pTheta = 0.016461;
      geom.pPhi = 0.731218;
      geom.pAlp1 = 0.000849;
      geom.pAlp2 = 0.000848;
      geom.pDy1 = 1.124739;
      geom.pDx1 = 1.190025;
      geom.pDx2 = 1.191935;
      geom.pDy2 = 1.273242;
      geom.pDx3 = 1.355380;
      geom.pDx4 = 1.357540;
      geom.centralX = 1.276220;
      geom.centralY = 105.070593;
      geom.centralZ = 6.094317;
      geom.pRotationAngleX = -1.501347;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1036 based Row/Col = 103/6
      geom_tower geom;
      geom.id = 1036;
      geom.pDz = 6.751968;
      geom.pTheta = 0.038366;
      geom.pPhi = 0.290464;
      geom.pAlp1 = 0.002550;
      geom.pAlp2 = 0.002548;
      geom.pDy1 = 1.124739;
      geom.pDx1 = 1.191456;
      geom.pDx2 = 1.193373;
      geom.pDy2 = 1.273242;
      geom.pDx3 = 1.357009;
      geom.pDx4 = 1.359178;
      geom.centralX = 3.830194;
      geom.centralY = 105.070593;
      geom.centralZ = 6.094317;
      geom.pRotationAngleX = -1.501347;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1037 based Row/Col = 103/7
      geom_tower geom;
      geom.id = 1037;
      geom.pDz = 6.751968;
      geom.pTheta = 0.062242;
      geom.pPhi = 0.177325;
      geom.pAlp1 = 0.004259;
      geom.pAlp2 = 0.004257;
      geom.pDy1 = 1.124739;
      geom.pDx1 = 1.194326;
      geom.pDx2 = 1.196256;
      geom.pDy2 = 1.273242;
      geom.pDx3 = 1.360277;
      geom.pDx4 = 1.362460;
      geom.centralX = 6.388778;
      geom.centralY = 105.070593;
      geom.centralZ = 6.094317;
      geom.pRotationAngleX = -1.501347;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1038 based Row/Col = 103/8
      geom_tower geom;
      geom.id = 1038;
      geom.pDz = 6.751968;
      geom.pTheta = 0.086470;
      geom.pPhi = 0.127161;
      geom.pAlp1 = 0.005984;
      geom.pAlp2 = 0.005981;
      geom.pDy1 = 1.124739;
      geom.pDx1 = 1.198647;
      geom.pDx2 = 1.200597;
      geom.pDy2 = 1.273242;
      geom.pDx3 = 1.365197;
      geom.pDx4 = 1.367404;
      geom.centralX = 8.955069;
      geom.centralY = 105.070593;
      geom.centralZ = 6.094317;
      geom.pRotationAngleX = -1.501347;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 971 based Row/Col = 97/1
      geom_tower geom;
      geom.id = 971;
      geom.pDz = 6.751986;
      geom.pTheta = 0.086342;
      geom.pPhi = 3.014158;
      geom.pAlp1 = 0.005984;
      geom.pAlp2 = 0.005981;
      geom.pDy1 = 1.124739;
      geom.pDx1 = 1.198647;
      geom.pDx2 = 1.196698;
      geom.pDy2 = 1.273242;
      geom.pDx3 = 1.365197;
      geom.pDx4 = 1.362991;
      geom.centralX = -8.940723;
      geom.centralY = 104.903842;
      geom.centralZ = -8.491506;
      geom.pRotationAngleX = -1.640246;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 972 based Row/Col = 97/2
      geom_tower geom;
      geom.id = 972;
      geom.pDz = 6.751986;
      geom.pTheta = 0.062151;
      geom.pPhi = 2.963891;
      geom.pAlp1 = 0.004259;
      geom.pAlp2 = 0.004257;
      geom.pDy1 = 1.124739;
      geom.pDx1 = 1.194326;
      geom.pDx2 = 1.192396;
      geom.pDy2 = 1.273242;
      geom.pDx3 = 1.360277;
      geom.pDx4 = 1.358093;
      geom.centralX = -6.378566;
      geom.centralY = 104.903842;
      geom.centralZ = -8.491506;
      geom.pRotationAngleX = -1.640246;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 973 based Row/Col = 97/3
      geom_tower geom;
      geom.id = 973;
      geom.pDz = 6.751986;
      geom.pTheta = 0.038315;
      geom.pPhi = 2.850534;
      geom.pAlp1 = 0.002550;
      geom.pAlp2 = 0.002548;
      geom.pDy1 = 1.124739;
      geom.pDx1 = 1.191456;
      geom.pDx2 = 1.189540;
      geom.pDy2 = 1.273242;
      geom.pDx3 = 1.357009;
      geom.pDx4 = 1.354841;
      geom.centralX = -3.824082;
      geom.centralY = 104.903842;
      geom.centralZ = -8.491506;
      geom.pRotationAngleX = -1.640246;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 974 based Row/Col = 97/4
      geom_tower geom;
      geom.id = 974;
      geom.pDz = 6.751986;
      geom.pTheta = 0.016452;
      geom.pPhi = 2.409299;
      geom.pAlp1 = 0.000849;
      geom.pAlp2 = 0.000848;
      geom.pDy1 = 1.124739;
      geom.pDx1 = 1.190025;
      geom.pDx2 = 1.188116;
      geom.pDy2 = 1.273242;
      geom.pDx3 = 1.355380;
      geom.pDx4 = 1.353219;
      geom.centralX = -1.274185;
      geom.centralY = 104.903842;
      geom.centralZ = -8.491506;
      geom.pRotationAngleX = -1.640246;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 975 based Row/Col = 97/5
      geom_tower geom;
      geom.id = 975;
      geom.pDz = 6.751986;
      geom.pTheta = 0.016452;
      geom.pPhi = 0.732294;
      geom.pAlp1 = -0.000849;
      geom.pAlp2 = -0.000848;
      geom.pDy1 = 1.124739;
      geom.pDx1 = 1.190025;
      geom.pDx2 = 1.188116;
      geom.pDy2 = 1.273242;
      geom.pDx3 = 1.355380;
      geom.pDx4 = 1.353219;
      geom.centralX = 1.274185;
      geom.centralY = 104.903842;
      geom.centralZ = -8.491506;
      geom.pRotationAngleX = -1.640246;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 976 based Row/Col = 97/6
      geom_tower geom;
      geom.id = 976;
      geom.pDz = 6.751986;
      geom.pTheta = 0.038315;
      geom.pPhi = 0.291058;
      geom.pAlp1 = -0.002550;
      geom.pAlp2 = -0.002548;
      geom.pDy1 = 1.124739;
      geom.pDx1 = 1.191456;
      geom.pDx2 = 1.189540;
      geom.pDy2 = 1.273242;
      geom.pDx3 = 1.357009;
      geom.pDx4 = 1.354841;
      geom.centralX = 3.824082;
      geom.centralY = 104.903842;
      geom.centralZ = -8.491506;
      geom.pRotationAngleX = -1.640246;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 977 based Row/Col = 97/7
      geom_tower geom;
      geom.id = 977;
      geom.pDz = 6.751986;
      geom.pTheta = 0.062151;
      geom.pPhi = 0.177701;
      geom.pAlp1 = -0.004259;
      geom.pAlp2 = -0.004257;
      geom.pDy1 = 1.124739;
      geom.pDx1 = 1.194326;
      geom.pDx2 = 1.192396;
      geom.pDy2 = 1.273242;
      geom.pDx3 = 1.360277;
      geom.pDx4 = 1.358093;
      geom.centralX = 6.378566;
      geom.centralY = 104.903842;
      geom.centralZ = -8.491506;
      geom.pRotationAngleX = -1.640246;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 978 based Row/Col = 97/8
      geom_tower geom;
      geom.id = 978;
      geom.pDz = 6.751986;
      geom.pTheta = 0.086342;
      geom.pPhi = 0.127434;
      geom.pAlp1 = -0.005984;
      geom.pAlp2 = -0.005981;
      geom.pDy1 = 1.124739;
      geom.pDx1 = 1.198647;
      geom.pDx2 = 1.196698;
      geom.pDy2 = 1.273242;
      geom.pDx3 = 1.365197;
      geom.pDx4 = 1.362991;
      geom.centralX = 8.940723;
      geom.centralY = 104.903842;
      geom.centralZ = -8.491506;
      geom.pRotationAngleX = -1.640246;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 981 based Row/Col = 98/1
      geom_tower geom;
      geom.id = 981;
      geom.pDz = 6.751968;
      geom.pTheta = 0.086470;
      geom.pPhi = -3.014432;
      geom.pAlp1 = 0.005984;
      geom.pAlp2 = 0.005981;
      geom.pDy1 = 1.124739;
      geom.pDx1 = 1.200597;
      geom.pDx2 = 1.198647;
      geom.pDy2 = 1.273242;
      geom.pDx3 = 1.367404;
      geom.pDx4 = 1.365197;
      geom.centralX = -8.955069;
      geom.centralY = 105.070593;
      geom.centralZ = -6.094317;
      geom.pRotationAngleX = -1.640246;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 982 based Row/Col = 98/2
      geom_tower geom;
      geom.id = 982;
      geom.pDz = 6.751968;
      geom.pTheta = 0.062242;
      geom.pPhi = -2.964268;
      geom.pAlp1 = 0.004259;
      geom.pAlp2 = 0.004257;
      geom.pDy1 = 1.124739;
      geom.pDx1 = 1.196256;
      geom.pDx2 = 1.194326;
      geom.pDy2 = 1.273242;
      geom.pDx3 = 1.362460;
      geom.pDx4 = 1.360277;
      geom.centralX = -6.388778;
      geom.centralY = 105.070593;
      geom.centralZ = -6.094317;
      geom.pRotationAngleX = -1.640246;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 983 based Row/Col = 98/3
      geom_tower geom;
      geom.id = 983;
      geom.pDz = 6.751968;
      geom.pTheta = 0.038366;
      geom.pPhi = -2.851129;
      geom.pAlp1 = 0.002550;
      geom.pAlp2 = 0.002548;
      geom.pDy1 = 1.124739;
      geom.pDx1 = 1.193373;
      geom.pDx2 = 1.191456;
      geom.pDy2 = 1.273242;
      geom.pDx3 = 1.359178;
      geom.pDx4 = 1.357009;
      geom.centralX = -3.830194;
      geom.centralY = 105.070593;
      geom.centralZ = -6.094317;
      geom.pRotationAngleX = -1.640246;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 984 based Row/Col = 98/4
      geom_tower geom;
      geom.id = 984;
      geom.pDz = 6.751968;
      geom.pTheta = 0.016461;
      geom.pPhi = -2.410375;
      geom.pAlp1 = 0.000849;
      geom.pAlp2 = 0.000848;
      geom.pDy1 = 1.124739;
      geom.pDx1 = 1.191935;
      geom.pDx2 = 1.190025;
      geom.pDy2 = 1.273242;
      geom.pDx3 = 1.357540;
      geom.pDx4 = 1.355380;
      geom.centralX = -1.276220;
      geom.centralY = 105.070593;
      geom.centralZ = -6.094317;
      geom.pRotationAngleX = -1.640246;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 985 based Row/Col = 98/5
      geom_tower geom;
      geom.id = 985;
      geom.pDz = 6.751968;
      geom.pTheta = 0.016461;
      geom.pPhi = -0.731218;
      geom.pAlp1 = -0.000849;
      geom.pAlp2 = -0.000848;
      geom.pDy1 = 1.124739;
      geom.pDx1 = 1.191935;
      geom.pDx2 = 1.190025;
      geom.pDy2 = 1.273242;
      geom.pDx3 = 1.357540;
      geom.pDx4 = 1.355380;
      geom.centralX = 1.276220;
      geom.centralY = 105.070593;
      geom.centralZ = -6.094317;
      geom.pRotationAngleX = -1.640246;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 986 based Row/Col = 98/6
      geom_tower geom;
      geom.id = 986;
      geom.pDz = 6.751968;
      geom.pTheta = 0.038366;
      geom.pPhi = -0.290464;
      geom.pAlp1 = -0.002550;
      geom.pAlp2 = -0.002548;
      geom.pDy1 = 1.124739;
      geom.pDx1 = 1.193373;
      geom.pDx2 = 1.191456;
      geom.pDy2 = 1.273242;
      geom.pDx3 = 1.359178;
      geom.pDx4 = 1.357009;
      geom.centralX = 3.830194;
      geom.centralY = 105.070593;
      geom.centralZ = -6.094317;
      geom.pRotationAngleX = -1.640246;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 987 based Row/Col = 98/7
      geom_tower geom;
      geom.id = 987;
      geom.pDz = 6.751968;
      geom.pTheta = 0.062242;
      geom.pPhi = -0.177325;
      geom.pAlp1 = -0.004259;
      geom.pAlp2 = -0.004257;
      geom.pDy1 = 1.124739;
      geom.pDx1 = 1.196256;
      geom.pDx2 = 1.194326;
      geom.pDy2 = 1.273242;
      geom.pDx3 = 1.362460;
      geom.pDx4 = 1.360277;
      geom.centralX = 6.388778;
      geom.centralY = 105.070593;
      geom.centralZ = -6.094317;
      geom.pRotationAngleX = -1.640246;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 988 based Row/Col = 98/8
      geom_tower geom;
      geom.id = 988;
      geom.pDz = 6.751968;
      geom.pTheta = 0.086470;
      geom.pPhi = -0.127161;
      geom.pAlp1 = -0.005984;
      geom.pAlp2 = -0.005981;
      geom.pDy1 = 1.124739;
      geom.pDx1 = 1.200597;
      geom.pDx2 = 1.198647;
      geom.pDy2 = 1.273242;
      geom.pDx3 = 1.367404;
      geom.pDx4 = 1.365197;
      geom.centralX = 8.955069;
      geom.centralY = 105.070593;
      geom.centralZ = -6.094317;
      geom.pRotationAngleX = -1.640246;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1081 based Row/Col = 108/1
      geom_tower geom;
      geom.id = 1081;
      geom.pDz = 6.751980;
      geom.pTheta = 0.085365;
      geom.pPhi = -3.013701;
      geom.pAlp1 = -0.013767;
      geom.pAlp2 = -0.013762;
      geom.pDy1 = 1.125316;
      geom.pDx1 = 1.192951;
      geom.pDx2 = 1.197438;
      geom.pDy2 = 1.272617;
      geom.pDx3 = 1.357161;
      geom.pDx4 = 1.362234;
      geom.centralX = -8.917743;
      geom.centralY = 104.635618;
      geom.centralZ = 18.145887;
      geom.pRotationAngleX = -1.410444;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1082 based Row/Col = 108/2
      geom_tower geom;
      geom.id = 1082;
      geom.pDz = 6.751980;
      geom.pTheta = 0.061452;
      geom.pPhi = -2.963264;
      geom.pAlp1 = -0.009801;
      geom.pAlp2 = -0.009797;
      geom.pDy1 = 1.125316;
      geom.pDx1 = 1.188768;
      geom.pDx2 = 1.193210;
      geom.pDy2 = 1.272617;
      geom.pDx3 = 1.352404;
      geom.pDx4 = 1.357425;
      geom.centralX = -6.362345;
      geom.centralY = 104.635618;
      geom.centralZ = 18.145887;
      geom.pRotationAngleX = -1.410444;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1083 based Row/Col = 108/3
      geom_tower geom;
      geom.id = 1083;
      geom.pDz = 6.751980;
      geom.pTheta = 0.037891;
      geom.pPhi = -2.849549;
      geom.pAlp1 = -0.005867;
      geom.pAlp2 = -0.005865;
      geom.pDy1 = 1.125316;
      geom.pDx1 = 1.185991;
      geom.pDx2 = 1.190402;
      geom.pDy2 = 1.272617;
      geom.pDx3 = 1.349245;
      geom.pDx4 = 1.354232;
      geom.centralX = -3.814425;
      geom.centralY = 104.635618;
      geom.centralZ = 18.145887;
      geom.pRotationAngleX = -1.410444;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1084 based Row/Col = 108/4
      geom_tower geom;
      geom.id = 1084;
      geom.pDz = 6.751980;
      geom.pTheta = 0.016291;
      geom.pPhi = -2.407524;
      geom.pAlp1 = -0.001954;
      geom.pAlp2 = -0.001953;
      geom.pDy1 = 1.125316;
      geom.pDx1 = 1.184605;
      geom.pDx2 = 1.189002;
      geom.pDy2 = 1.272617;
      geom.pDx3 = 1.347669;
      geom.pDx4 = 1.352639;
      geom.centralX = -1.270979;
      geom.centralY = 104.635618;
      geom.centralZ = 18.145887;
      geom.pRotationAngleX = -1.410444;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1085 based Row/Col = 108/5
      geom_tower geom;
      geom.id = 1085;
      geom.pDz = 6.751980;
      geom.pTheta = 0.016291;
      geom.pPhi = -0.734069;
      geom.pAlp1 = 0.001954;
      geom.pAlp2 = 0.001953;
      geom.pDy1 = 1.125316;
      geom.pDx1 = 1.184605;
      geom.pDx2 = 1.189002;
      geom.pDy2 = 1.272617;
      geom.pDx3 = 1.347669;
      geom.pDx4 = 1.352639;
      geom.centralX = 1.270979;
      geom.centralY = 104.635618;
      geom.centralZ = 18.145887;
      geom.pRotationAngleX = -1.410444;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1086 based Row/Col = 108/6
      geom_tower geom;
      geom.id = 1086;
      geom.pDz = 6.751980;
      geom.pTheta = 0.037891;
      geom.pPhi = -0.292044;
      geom.pAlp1 = 0.005867;
      geom.pAlp2 = 0.005865;
      geom.pDy1 = 1.125316;
      geom.pDx1 = 1.185991;
      geom.pDx2 = 1.190402;
      geom.pDy2 = 1.272617;
      geom.pDx3 = 1.349245;
      geom.pDx4 = 1.354232;
      geom.centralX = 3.814425;
      geom.centralY = 104.635618;
      geom.centralZ = 18.145887;
      geom.pRotationAngleX = -1.410444;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1087 based Row/Col = 108/7
      geom_tower geom;
      geom.id = 1087;
      geom.pDz = 6.751980;
      geom.pTheta = 0.061452;
      geom.pPhi = -0.178328;
      geom.pAlp1 = 0.009801;
      geom.pAlp2 = 0.009797;
      geom.pDy1 = 1.125316;
      geom.pDx1 = 1.188768;
      geom.pDx2 = 1.193210;
      geom.pDy2 = 1.272617;
      geom.pDx3 = 1.352404;
      geom.pDx4 = 1.357425;
      geom.centralX = 6.362345;
      geom.centralY = 104.635618;
      geom.centralZ = 18.145887;
      geom.pRotationAngleX = -1.410444;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1088 based Row/Col = 108/8
      geom_tower geom;
      geom.id = 1088;
      geom.pDz = 6.751980;
      geom.pTheta = 0.085365;
      geom.pPhi = -0.127892;
      geom.pAlp1 = 0.013767;
      geom.pAlp2 = 0.013762;
      geom.pDy1 = 1.125316;
      geom.pDx1 = 1.192951;
      geom.pDx2 = 1.197438;
      geom.pDy2 = 1.272617;
      geom.pDx3 = 1.357161;
      geom.pDx4 = 1.362234;
      geom.centralX = 8.917743;
      geom.centralY = 104.635618;
      geom.centralZ = 18.145887;
      geom.pRotationAngleX = -1.410444;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1071 based Row/Col = 107/1
      geom_tower geom;
      geom.id = 1071;
      geom.pDz = 6.751956;
      geom.pTheta = 0.085659;
      geom.pPhi = 3.014292;
      geom.pAlp1 = -0.013767;
      geom.pAlp2 = -0.013762;
      geom.pDy1 = 1.125316;
      geom.pDx1 = 1.197438;
      geom.pDx2 = 1.201925;
      geom.pDy2 = 1.272617;
      geom.pDx3 = 1.362234;
      geom.pDx4 = 1.367307;
      geom.centralX = -8.950751;
      geom.centralY = 105.019285;
      geom.centralZ = 15.773781;
      geom.pRotationAngleX = -1.410444;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1072 based Row/Col = 107/2
      geom_tower geom;
      geom.id = 1072;
      geom.pDz = 6.751956;
      geom.pTheta = 0.061659;
      geom.pPhi = 2.964079;
      geom.pAlp1 = -0.009800;
      geom.pAlp2 = -0.009796;
      geom.pDy1 = 1.125316;
      geom.pDx1 = 1.193210;
      geom.pDx2 = 1.197652;
      geom.pDy2 = 1.272617;
      geom.pDx3 = 1.357425;
      geom.pDx4 = 1.362447;
      geom.centralX = -6.385841;
      geom.centralY = 105.019285;
      geom.centralZ = 15.773781;
      geom.pRotationAngleX = -1.410444;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1073 based Row/Col = 107/3
      geom_tower geom;
      geom.id = 1073;
      geom.pDz = 6.751956;
      geom.pTheta = 0.038009;
      geom.pPhi = 2.850835;
      geom.pAlp1 = -0.005867;
      geom.pAlp2 = -0.005864;
      geom.pDy1 = 1.125316;
      geom.pDx1 = 1.190402;
      geom.pDx2 = 1.194814;
      geom.pDy2 = 1.272617;
      geom.pDx3 = 1.354232;
      geom.pDx4 = 1.359219;
      geom.centralX = -3.828491;
      geom.centralY = 105.019285;
      geom.centralZ = 15.773781;
      geom.pRotationAngleX = -1.410444;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1074 based Row/Col = 107/4
      geom_tower geom;
      geom.id = 1074;
      geom.pDz = 6.751956;
      geom.pTheta = 0.016314;
      geom.pPhi = 2.409846;
      geom.pAlp1 = -0.001953;
      geom.pAlp2 = -0.001953;
      geom.pDy1 = 1.125316;
      geom.pDx1 = 1.189002;
      geom.pDx2 = 1.193398;
      geom.pDy2 = 1.272617;
      geom.pDx3 = 1.352639;
      geom.pDx4 = 1.357609;
      geom.centralX = -1.275662;
      geom.centralY = 105.019285;
      geom.centralZ = 15.773781;
      geom.pRotationAngleX = -1.410444;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1075 based Row/Col = 107/5
      geom_tower geom;
      geom.id = 1075;
      geom.pDz = 6.751956;
      geom.pTheta = 0.016314;
      geom.pPhi = 0.731746;
      geom.pAlp1 = 0.001953;
      geom.pAlp2 = 0.001953;
      geom.pDy1 = 1.125316;
      geom.pDx1 = 1.189002;
      geom.pDx2 = 1.193398;
      geom.pDy2 = 1.272617;
      geom.pDx3 = 1.352639;
      geom.pDx4 = 1.357609;
      geom.centralX = 1.275662;
      geom.centralY = 105.019285;
      geom.centralZ = 15.773781;
      geom.pRotationAngleX = -1.410444;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1076 based Row/Col = 107/6
      geom_tower geom;
      geom.id = 1076;
      geom.pDz = 6.751956;
      geom.pTheta = 0.038009;
      geom.pPhi = 0.290758;
      geom.pAlp1 = 0.005867;
      geom.pAlp2 = 0.005864;
      geom.pDy1 = 1.125316;
      geom.pDx1 = 1.190402;
      geom.pDx2 = 1.194814;
      geom.pDy2 = 1.272617;
      geom.pDx3 = 1.354232;
      geom.pDx4 = 1.359219;
      geom.centralX = 3.828491;
      geom.centralY = 105.019285;
      geom.centralZ = 15.773781;
      geom.pRotationAngleX = -1.410444;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1077 based Row/Col = 107/7
      geom_tower geom;
      geom.id = 1077;
      geom.pDz = 6.751956;
      geom.pTheta = 0.061659;
      geom.pPhi = 0.177513;
      geom.pAlp1 = 0.009800;
      geom.pAlp2 = 0.009796;
      geom.pDy1 = 1.125316;
      geom.pDx1 = 1.193210;
      geom.pDx2 = 1.197652;
      geom.pDy2 = 1.272617;
      geom.pDx3 = 1.357425;
      geom.pDx4 = 1.362447;
      geom.centralX = 6.385841;
      geom.centralY = 105.019285;
      geom.centralZ = 15.773781;
      geom.pRotationAngleX = -1.410444;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1078 based Row/Col = 107/8
      geom_tower geom;
      geom.id = 1078;
      geom.pDz = 6.751956;
      geom.pTheta = 0.085659;
      geom.pPhi = 0.127301;
      geom.pAlp1 = 0.013767;
      geom.pAlp2 = 0.013762;
      geom.pDy1 = 1.125316;
      geom.pDx1 = 1.197438;
      geom.pDx2 = 1.201925;
      geom.pDy2 = 1.272617;
      geom.pDx3 = 1.362234;
      geom.pDx4 = 1.367307;
      geom.centralX = 8.950751;
      geom.centralY = 105.019285;
      geom.centralZ = 15.773781;
      geom.pRotationAngleX = -1.410444;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 931 based Row/Col = 93/1
      geom_tower geom;
      geom.id = 931;
      geom.pDz = 6.751980;
      geom.pTheta = 0.085365;
      geom.pPhi = 3.013701;
      geom.pAlp1 = 0.013767;
      geom.pAlp2 = 0.013762;
      geom.pDy1 = 1.125316;
      geom.pDx1 = 1.197438;
      geom.pDx2 = 1.192951;
      geom.pDy2 = 1.272617;
      geom.pDx3 = 1.362234;
      geom.pDx4 = 1.357161;
      geom.centralX = -8.917743;
      geom.centralY = 104.635618;
      geom.centralZ = -18.145887;
      geom.pRotationAngleX = -1.731149;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 932 based Row/Col = 93/2
      geom_tower geom;
      geom.id = 932;
      geom.pDz = 6.751980;
      geom.pTheta = 0.061452;
      geom.pPhi = 2.963264;
      geom.pAlp1 = 0.009801;
      geom.pAlp2 = 0.009797;
      geom.pDy1 = 1.125316;
      geom.pDx1 = 1.193210;
      geom.pDx2 = 1.188768;
      geom.pDy2 = 1.272617;
      geom.pDx3 = 1.357425;
      geom.pDx4 = 1.352404;
      geom.centralX = -6.362345;
      geom.centralY = 104.635618;
      geom.centralZ = -18.145887;
      geom.pRotationAngleX = -1.731149;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 933 based Row/Col = 93/3
      geom_tower geom;
      geom.id = 933;
      geom.pDz = 6.751980;
      geom.pTheta = 0.037891;
      geom.pPhi = 2.849549;
      geom.pAlp1 = 0.005867;
      geom.pAlp2 = 0.005865;
      geom.pDy1 = 1.125316;
      geom.pDx1 = 1.190402;
      geom.pDx2 = 1.185991;
      geom.pDy2 = 1.272617;
      geom.pDx3 = 1.354232;
      geom.pDx4 = 1.349245;
      geom.centralX = -3.814425;
      geom.centralY = 104.635618;
      geom.centralZ = -18.145887;
      geom.pRotationAngleX = -1.731149;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 934 based Row/Col = 93/4
      geom_tower geom;
      geom.id = 934;
      geom.pDz = 6.751980;
      geom.pTheta = 0.016291;
      geom.pPhi = 2.407524;
      geom.pAlp1 = 0.001954;
      geom.pAlp2 = 0.001953;
      geom.pDy1 = 1.125316;
      geom.pDx1 = 1.189002;
      geom.pDx2 = 1.184605;
      geom.pDy2 = 1.272617;
      geom.pDx3 = 1.352639;
      geom.pDx4 = 1.347669;
      geom.centralX = -1.270979;
      geom.centralY = 104.635618;
      geom.centralZ = -18.145887;
      geom.pRotationAngleX = -1.731149;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 935 based Row/Col = 93/5
      geom_tower geom;
      geom.id = 935;
      geom.pDz = 6.751980;
      geom.pTheta = 0.016291;
      geom.pPhi = 0.734069;
      geom.pAlp1 = -0.001954;
      geom.pAlp2 = -0.001953;
      geom.pDy1 = 1.125316;
      geom.pDx1 = 1.189002;
      geom.pDx2 = 1.184605;
      geom.pDy2 = 1.272617;
      geom.pDx3 = 1.352639;
      geom.pDx4 = 1.347669;
      geom.centralX = 1.270979;
      geom.centralY = 104.635618;
      geom.centralZ = -18.145887;
      geom.pRotationAngleX = -1.731149;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 936 based Row/Col = 93/6
      geom_tower geom;
      geom.id = 936;
      geom.pDz = 6.751980;
      geom.pTheta = 0.037891;
      geom.pPhi = 0.292044;
      geom.pAlp1 = -0.005867;
      geom.pAlp2 = -0.005865;
      geom.pDy1 = 1.125316;
      geom.pDx1 = 1.190402;
      geom.pDx2 = 1.185991;
      geom.pDy2 = 1.272617;
      geom.pDx3 = 1.354232;
      geom.pDx4 = 1.349245;
      geom.centralX = 3.814425;
      geom.centralY = 104.635618;
      geom.centralZ = -18.145887;
      geom.pRotationAngleX = -1.731149;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 937 based Row/Col = 93/7
      geom_tower geom;
      geom.id = 937;
      geom.pDz = 6.751980;
      geom.pTheta = 0.061452;
      geom.pPhi = 0.178328;
      geom.pAlp1 = -0.009801;
      geom.pAlp2 = -0.009797;
      geom.pDy1 = 1.125316;
      geom.pDx1 = 1.193210;
      geom.pDx2 = 1.188768;
      geom.pDy2 = 1.272617;
      geom.pDx3 = 1.357425;
      geom.pDx4 = 1.352404;
      geom.centralX = 6.362345;
      geom.centralY = 104.635618;
      geom.centralZ = -18.145887;
      geom.pRotationAngleX = -1.731149;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 938 based Row/Col = 93/8
      geom_tower geom;
      geom.id = 938;
      geom.pDz = 6.751980;
      geom.pTheta = 0.085365;
      geom.pPhi = 0.127892;
      geom.pAlp1 = -0.013767;
      geom.pAlp2 = -0.013762;
      geom.pDy1 = 1.125316;
      geom.pDx1 = 1.197438;
      geom.pDx2 = 1.192951;
      geom.pDy2 = 1.272617;
      geom.pDx3 = 1.362234;
      geom.pDx4 = 1.357161;
      geom.centralX = 8.917743;
      geom.centralY = 104.635618;
      geom.centralZ = -18.145887;
      geom.pRotationAngleX = -1.731149;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 941 based Row/Col = 94/1
      geom_tower geom;
      geom.id = 941;
      geom.pDz = 6.751956;
      geom.pTheta = 0.085659;
      geom.pPhi = -3.014292;
      geom.pAlp1 = 0.013767;
      geom.pAlp2 = 0.013762;
      geom.pDy1 = 1.125316;
      geom.pDx1 = 1.201925;
      geom.pDx2 = 1.197438;
      geom.pDy2 = 1.272617;
      geom.pDx3 = 1.367307;
      geom.pDx4 = 1.362234;
      geom.centralX = -8.950751;
      geom.centralY = 105.019285;
      geom.centralZ = -15.773781;
      geom.pRotationAngleX = -1.731149;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 942 based Row/Col = 94/2
      geom_tower geom;
      geom.id = 942;
      geom.pDz = 6.751956;
      geom.pTheta = 0.061659;
      geom.pPhi = -2.964079;
      geom.pAlp1 = 0.009800;
      geom.pAlp2 = 0.009796;
      geom.pDy1 = 1.125316;
      geom.pDx1 = 1.197652;
      geom.pDx2 = 1.193210;
      geom.pDy2 = 1.272617;
      geom.pDx3 = 1.362447;
      geom.pDx4 = 1.357425;
      geom.centralX = -6.385841;
      geom.centralY = 105.019285;
      geom.centralZ = -15.773781;
      geom.pRotationAngleX = -1.731149;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 943 based Row/Col = 94/3
      geom_tower geom;
      geom.id = 943;
      geom.pDz = 6.751956;
      geom.pTheta = 0.038009;
      geom.pPhi = -2.850835;
      geom.pAlp1 = 0.005867;
      geom.pAlp2 = 0.005864;
      geom.pDy1 = 1.125316;
      geom.pDx1 = 1.194814;
      geom.pDx2 = 1.190402;
      geom.pDy2 = 1.272617;
      geom.pDx3 = 1.359219;
      geom.pDx4 = 1.354232;
      geom.centralX = -3.828491;
      geom.centralY = 105.019285;
      geom.centralZ = -15.773781;
      geom.pRotationAngleX = -1.731149;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 944 based Row/Col = 94/4
      geom_tower geom;
      geom.id = 944;
      geom.pDz = 6.751956;
      geom.pTheta = 0.016314;
      geom.pPhi = -2.409846;
      geom.pAlp1 = 0.001953;
      geom.pAlp2 = 0.001953;
      geom.pDy1 = 1.125316;
      geom.pDx1 = 1.193398;
      geom.pDx2 = 1.189002;
      geom.pDy2 = 1.272617;
      geom.pDx3 = 1.357609;
      geom.pDx4 = 1.352639;
      geom.centralX = -1.275662;
      geom.centralY = 105.019285;
      geom.centralZ = -15.773781;
      geom.pRotationAngleX = -1.731149;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 945 based Row/Col = 94/5
      geom_tower geom;
      geom.id = 945;
      geom.pDz = 6.751956;
      geom.pTheta = 0.016314;
      geom.pPhi = -0.731746;
      geom.pAlp1 = -0.001953;
      geom.pAlp2 = -0.001953;
      geom.pDy1 = 1.125316;
      geom.pDx1 = 1.193398;
      geom.pDx2 = 1.189002;
      geom.pDy2 = 1.272617;
      geom.pDx3 = 1.357609;
      geom.pDx4 = 1.352639;
      geom.centralX = 1.275662;
      geom.centralY = 105.019285;
      geom.centralZ = -15.773781;
      geom.pRotationAngleX = -1.731149;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 946 based Row/Col = 94/6
      geom_tower geom;
      geom.id = 946;
      geom.pDz = 6.751956;
      geom.pTheta = 0.038009;
      geom.pPhi = -0.290758;
      geom.pAlp1 = -0.005867;
      geom.pAlp2 = -0.005864;
      geom.pDy1 = 1.125316;
      geom.pDx1 = 1.194814;
      geom.pDx2 = 1.190402;
      geom.pDy2 = 1.272617;
      geom.pDx3 = 1.359219;
      geom.pDx4 = 1.354232;
      geom.centralX = 3.828491;
      geom.centralY = 105.019285;
      geom.centralZ = -15.773781;
      geom.pRotationAngleX = -1.731149;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 947 based Row/Col = 94/7
      geom_tower geom;
      geom.id = 947;
      geom.pDz = 6.751956;
      geom.pTheta = 0.061659;
      geom.pPhi = -0.177513;
      geom.pAlp1 = -0.009800;
      geom.pAlp2 = -0.009796;
      geom.pDy1 = 1.125316;
      geom.pDx1 = 1.197652;
      geom.pDx2 = 1.193210;
      geom.pDy2 = 1.272617;
      geom.pDx3 = 1.362447;
      geom.pDx4 = 1.357425;
      geom.centralX = 6.385841;
      geom.centralY = 105.019285;
      geom.centralZ = -15.773781;
      geom.pRotationAngleX = -1.731149;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 948 based Row/Col = 94/8
      geom_tower geom;
      geom.id = 948;
      geom.pDz = 6.751956;
      geom.pTheta = 0.085659;
      geom.pPhi = -0.127301;
      geom.pAlp1 = -0.013767;
      geom.pAlp2 = -0.013762;
      geom.pDy1 = 1.125316;
      geom.pDx1 = 1.201925;
      geom.pDx2 = 1.197438;
      geom.pDy2 = 1.272617;
      geom.pDx3 = 1.367307;
      geom.pDx4 = 1.362234;
      geom.centralX = 8.950751;
      geom.centralY = 105.019285;
      geom.centralZ = -15.773781;
      geom.pRotationAngleX = -1.731149;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1121 based Row/Col = 112/1
      geom_tower geom;
      geom.id = 1121;
      geom.pDz = 6.751955;
      geom.pTheta = 0.083707;
      geom.pPhi = -3.014083;
      geom.pAlp1 = -0.021345;
      geom.pAlp2 = -0.021342;
      geom.pDy1 = 1.126302;
      geom.pDx1 = 1.190598;
      geom.pDx2 = 1.197557;
      geom.pDy2 = 1.270041;
      geom.pDx3 = 1.351425;
      geom.pDx4 = 1.359272;
      geom.centralX = -8.900303;
      geom.centralY = 104.430798;
      geom.centralZ = 27.928983;
      geom.pRotationAngleX = -1.320603;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1122 based Row/Col = 112/2
      geom_tower geom;
      geom.id = 1122;
      geom.pDz = 6.751955;
      geom.pTheta = 0.060255;
      geom.pPhi = -2.963800;
      geom.pAlp1 = -0.015198;
      geom.pAlp2 = -0.015196;
      geom.pDy1 = 1.126302;
      geom.pDx1 = 1.186591;
      geom.pDx2 = 1.193483;
      geom.pDy2 = 1.270041;
      geom.pDx3 = 1.346879;
      geom.pDx4 = 1.354649;
      geom.centralX = -6.350189;
      geom.centralY = 104.430798;
      geom.centralZ = 27.928983;
      geom.pRotationAngleX = -1.320603;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1123 based Row/Col = 112/3
      geom_tower geom;
      geom.id = 1123;
      geom.pDz = 6.751955;
      geom.pTheta = 0.037147;
      geom.pPhi = -2.850403;
      geom.pAlp1 = -0.009100;
      geom.pAlp2 = -0.009098;
      geom.pDy1 = 1.126302;
      geom.pDx1 = 1.183930;
      geom.pDx2 = 1.190778;
      geom.pDy2 = 1.270041;
      geom.pDx3 = 1.343859;
      geom.pDx4 = 1.351579;
      geom.centralX = -3.807252;
      geom.centralY = 104.430798;
      geom.centralZ = 27.928983;
      geom.pRotationAngleX = -1.320603;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1124 based Row/Col = 112/4
      geom_tower geom;
      geom.id = 1124;
      geom.pDz = 6.751955;
      geom.pTheta = 0.015953;
      geom.pPhi = -2.409074;
      geom.pAlp1 = -0.003030;
      geom.pAlp2 = -0.003029;
      geom.pDy1 = 1.126302;
      geom.pDx1 = 1.182603;
      geom.pDx2 = 1.189428;
      geom.pDy2 = 1.270041;
      geom.pDx3 = 1.342352;
      geom.pDx4 = 1.350047;
      geom.centralX = -1.268608;
      geom.centralY = 104.430798;
      geom.centralZ = 27.928983;
      geom.pRotationAngleX = -1.320603;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1125 based Row/Col = 112/5
      geom_tower geom;
      geom.id = 1125;
      geom.pDz = 6.751955;
      geom.pTheta = 0.015953;
      geom.pPhi = -0.732518;
      geom.pAlp1 = 0.003030;
      geom.pAlp2 = 0.003029;
      geom.pDy1 = 1.126302;
      geom.pDx1 = 1.182603;
      geom.pDx2 = 1.189428;
      geom.pDy2 = 1.270041;
      geom.pDx3 = 1.342352;
      geom.pDx4 = 1.350047;
      geom.centralX = 1.268608;
      geom.centralY = 104.430798;
      geom.centralZ = 27.928983;
      geom.pRotationAngleX = -1.320603;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1126 based Row/Col = 112/6
      geom_tower geom;
      geom.id = 1126;
      geom.pDz = 6.751955;
      geom.pTheta = 0.037147;
      geom.pPhi = -0.291189;
      geom.pAlp1 = 0.009100;
      geom.pAlp2 = 0.009098;
      geom.pDy1 = 1.126302;
      geom.pDx1 = 1.183930;
      geom.pDx2 = 1.190778;
      geom.pDy2 = 1.270041;
      geom.pDx3 = 1.343859;
      geom.pDx4 = 1.351579;
      geom.centralX = 3.807252;
      geom.centralY = 104.430798;
      geom.centralZ = 27.928983;
      geom.pRotationAngleX = -1.320603;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1127 based Row/Col = 112/7
      geom_tower geom;
      geom.id = 1127;
      geom.pDz = 6.751955;
      geom.pTheta = 0.060255;
      geom.pPhi = -0.177793;
      geom.pAlp1 = 0.015198;
      geom.pAlp2 = 0.015196;
      geom.pDy1 = 1.126302;
      geom.pDx1 = 1.186591;
      geom.pDx2 = 1.193483;
      geom.pDy2 = 1.270041;
      geom.pDx3 = 1.346879;
      geom.pDx4 = 1.354649;
      geom.centralX = 6.350189;
      geom.centralY = 104.430798;
      geom.centralZ = 27.928983;
      geom.pRotationAngleX = -1.320603;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1128 based Row/Col = 112/8
      geom_tower geom;
      geom.id = 1128;
      geom.pDz = 6.751955;
      geom.pTheta = 0.083707;
      geom.pPhi = -0.127510;
      geom.pAlp1 = 0.021345;
      geom.pAlp2 = 0.021342;
      geom.pDy1 = 1.126302;
      geom.pDx1 = 1.190598;
      geom.pDx2 = 1.197557;
      geom.pDy2 = 1.270041;
      geom.pDx3 = 1.351425;
      geom.pDx4 = 1.359272;
      geom.centralX = 8.900303;
      geom.centralY = 104.430798;
      geom.centralZ = 27.928983;
      geom.pRotationAngleX = -1.320603;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1111 based Row/Col = 111/1
      geom_tower geom;
      geom.id = 1111;
      geom.pDz = 6.751991;
      geom.pTheta = 0.084148;
      geom.pPhi = 3.015360;
      geom.pAlp1 = -0.021345;
      geom.pAlp2 = -0.021342;
      geom.pDy1 = 1.126302;
      geom.pDx1 = 1.197557;
      geom.pDx2 = 1.204518;
      geom.pDy2 = 1.270041;
      geom.pDx3 = 1.359272;
      geom.pDx4 = 1.367120;
      geom.centralX = -8.951457;
      geom.centralY = 105.025349;
      geom.centralZ = 25.602406;
      geom.pRotationAngleX = -1.320603;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1112 based Row/Col = 111/2
      geom_tower geom;
      geom.id = 1112;
      geom.pDz = 6.751991;
      geom.pTheta = 0.060563;
      geom.pPhi = 2.965561;
      geom.pAlp1 = -0.015197;
      geom.pAlp2 = -0.015195;
      geom.pDy1 = 1.126302;
      geom.pDx1 = 1.193483;
      geom.pDx2 = 1.200376;
      geom.pDy2 = 1.270041;
      geom.pDx3 = 1.354649;
      geom.pDx4 = 1.362420;
      geom.centralX = -6.386608;
      geom.centralY = 105.025349;
      geom.centralZ = 25.602406;
      geom.pRotationAngleX = -1.320603;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1113 based Row/Col = 111/3
      geom_tower geom;
      geom.id = 1113;
      geom.pDz = 6.751991;
      geom.pTheta = 0.037318;
      geom.pPhi = 2.853184;
      geom.pAlp1 = -0.009099;
      geom.pAlp2 = -0.009097;
      geom.pDy1 = 1.126302;
      geom.pDx1 = 1.190778;
      geom.pDx2 = 1.197625;
      geom.pDy2 = 1.270041;
      geom.pDx3 = 1.351579;
      geom.pDx4 = 1.359299;
      geom.centralX = -3.829055;
      geom.centralY = 105.025349;
      geom.centralZ = 25.602406;
      geom.pRotationAngleX = -1.320603;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1114 based Row/Col = 111/4
      geom_tower geom;
      geom.id = 1114;
      geom.pDz = 6.751991;
      geom.pTheta = 0.015967;
      geom.pPhi = 2.414118;
      geom.pAlp1 = -0.003030;
      geom.pAlp2 = -0.003029;
      geom.pDy1 = 1.126302;
      geom.pDx1 = 1.189428;
      geom.pDx2 = 1.196253;
      geom.pDy2 = 1.270041;
      geom.pDx3 = 1.350048;
      geom.pDx4 = 1.357742;
      geom.centralX = -1.275868;
      geom.centralY = 105.025349;
      geom.centralZ = 25.602406;
      geom.pRotationAngleX = -1.320603;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1115 based Row/Col = 111/5
      geom_tower geom;
      geom.id = 1115;
      geom.pDz = 6.751991;
      geom.pTheta = 0.015967;
      geom.pPhi = 0.727474;
      geom.pAlp1 = 0.003030;
      geom.pAlp2 = 0.003029;
      geom.pDy1 = 1.126302;
      geom.pDx1 = 1.189428;
      geom.pDx2 = 1.196253;
      geom.pDy2 = 1.270041;
      geom.pDx3 = 1.350048;
      geom.pDx4 = 1.357742;
      geom.centralX = 1.275868;
      geom.centralY = 105.025349;
      geom.centralZ = 25.602406;
      geom.pRotationAngleX = -1.320603;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1116 based Row/Col = 111/6
      geom_tower geom;
      geom.id = 1116;
      geom.pDz = 6.751991;
      geom.pTheta = 0.037318;
      geom.pPhi = 0.288409;
      geom.pAlp1 = 0.009099;
      geom.pAlp2 = 0.009097;
      geom.pDy1 = 1.126302;
      geom.pDx1 = 1.190778;
      geom.pDx2 = 1.197625;
      geom.pDy2 = 1.270041;
      geom.pDx3 = 1.351579;
      geom.pDx4 = 1.359299;
      geom.centralX = 3.829055;
      geom.centralY = 105.025349;
      geom.centralZ = 25.602406;
      geom.pRotationAngleX = -1.320603;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1117 based Row/Col = 111/7
      geom_tower geom;
      geom.id = 1117;
      geom.pDz = 6.751991;
      geom.pTheta = 0.060563;
      geom.pPhi = 0.176032;
      geom.pAlp1 = 0.015197;
      geom.pAlp2 = 0.015195;
      geom.pDy1 = 1.126302;
      geom.pDx1 = 1.193483;
      geom.pDx2 = 1.200376;
      geom.pDy2 = 1.270041;
      geom.pDx3 = 1.354649;
      geom.pDx4 = 1.362420;
      geom.centralX = 6.386608;
      geom.centralY = 105.025349;
      geom.centralZ = 25.602406;
      geom.pRotationAngleX = -1.320603;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1118 based Row/Col = 111/8
      geom_tower geom;
      geom.id = 1118;
      geom.pDz = 6.751991;
      geom.pTheta = 0.084148;
      geom.pPhi = 0.126233;
      geom.pAlp1 = 0.021345;
      geom.pAlp2 = 0.021342;
      geom.pDy1 = 1.126302;
      geom.pDx1 = 1.197557;
      geom.pDx2 = 1.204518;
      geom.pDy2 = 1.270041;
      geom.pDx3 = 1.359272;
      geom.pDx4 = 1.367120;
      geom.centralX = 8.951457;
      geom.centralY = 105.025349;
      geom.centralZ = 25.602406;
      geom.pRotationAngleX = -1.320603;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 891 based Row/Col = 89/1
      geom_tower geom;
      geom.id = 891;
      geom.pDz = 6.751955;
      geom.pTheta = 0.083707;
      geom.pPhi = 3.014083;
      geom.pAlp1 = 0.021345;
      geom.pAlp2 = 0.021342;
      geom.pDy1 = 1.126302;
      geom.pDx1 = 1.197557;
      geom.pDx2 = 1.190598;
      geom.pDy2 = 1.270041;
      geom.pDx3 = 1.359272;
      geom.pDx4 = 1.351425;
      geom.centralX = -8.900303;
      geom.centralY = 104.430798;
      geom.centralZ = -27.928983;
      geom.pRotationAngleX = -1.820989;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 892 based Row/Col = 89/2
      geom_tower geom;
      geom.id = 892;
      geom.pDz = 6.751955;
      geom.pTheta = 0.060255;
      geom.pPhi = 2.963800;
      geom.pAlp1 = 0.015198;
      geom.pAlp2 = 0.015196;
      geom.pDy1 = 1.126302;
      geom.pDx1 = 1.193483;
      geom.pDx2 = 1.186591;
      geom.pDy2 = 1.270041;
      geom.pDx3 = 1.354649;
      geom.pDx4 = 1.346879;
      geom.centralX = -6.350189;
      geom.centralY = 104.430798;
      geom.centralZ = -27.928983;
      geom.pRotationAngleX = -1.820989;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 893 based Row/Col = 89/3
      geom_tower geom;
      geom.id = 893;
      geom.pDz = 6.751955;
      geom.pTheta = 0.037147;
      geom.pPhi = 2.850403;
      geom.pAlp1 = 0.009100;
      geom.pAlp2 = 0.009098;
      geom.pDy1 = 1.126302;
      geom.pDx1 = 1.190778;
      geom.pDx2 = 1.183930;
      geom.pDy2 = 1.270041;
      geom.pDx3 = 1.351579;
      geom.pDx4 = 1.343859;
      geom.centralX = -3.807252;
      geom.centralY = 104.430798;
      geom.centralZ = -27.928983;
      geom.pRotationAngleX = -1.820989;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 894 based Row/Col = 89/4
      geom_tower geom;
      geom.id = 894;
      geom.pDz = 6.751955;
      geom.pTheta = 0.015953;
      geom.pPhi = 2.409074;
      geom.pAlp1 = 0.003030;
      geom.pAlp2 = 0.003029;
      geom.pDy1 = 1.126302;
      geom.pDx1 = 1.189428;
      geom.pDx2 = 1.182603;
      geom.pDy2 = 1.270041;
      geom.pDx3 = 1.350047;
      geom.pDx4 = 1.342352;
      geom.centralX = -1.268608;
      geom.centralY = 104.430798;
      geom.centralZ = -27.928983;
      geom.pRotationAngleX = -1.820989;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 895 based Row/Col = 89/5
      geom_tower geom;
      geom.id = 895;
      geom.pDz = 6.751955;
      geom.pTheta = 0.015953;
      geom.pPhi = 0.732518;
      geom.pAlp1 = -0.003030;
      geom.pAlp2 = -0.003029;
      geom.pDy1 = 1.126302;
      geom.pDx1 = 1.189428;
      geom.pDx2 = 1.182603;
      geom.pDy2 = 1.270041;
      geom.pDx3 = 1.350047;
      geom.pDx4 = 1.342352;
      geom.centralX = 1.268608;
      geom.centralY = 104.430798;
      geom.centralZ = -27.928983;
      geom.pRotationAngleX = -1.820989;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 896 based Row/Col = 89/6
      geom_tower geom;
      geom.id = 896;
      geom.pDz = 6.751955;
      geom.pTheta = 0.037147;
      geom.pPhi = 0.291189;
      geom.pAlp1 = -0.009100;
      geom.pAlp2 = -0.009098;
      geom.pDy1 = 1.126302;
      geom.pDx1 = 1.190778;
      geom.pDx2 = 1.183930;
      geom.pDy2 = 1.270041;
      geom.pDx3 = 1.351579;
      geom.pDx4 = 1.343859;
      geom.centralX = 3.807252;
      geom.centralY = 104.430798;
      geom.centralZ = -27.928983;
      geom.pRotationAngleX = -1.820989;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 897 based Row/Col = 89/7
      geom_tower geom;
      geom.id = 897;
      geom.pDz = 6.751955;
      geom.pTheta = 0.060255;
      geom.pPhi = 0.177793;
      geom.pAlp1 = -0.015198;
      geom.pAlp2 = -0.015196;
      geom.pDy1 = 1.126302;
      geom.pDx1 = 1.193483;
      geom.pDx2 = 1.186591;
      geom.pDy2 = 1.270041;
      geom.pDx3 = 1.354649;
      geom.pDx4 = 1.346879;
      geom.centralX = 6.350189;
      geom.centralY = 104.430798;
      geom.centralZ = -27.928983;
      geom.pRotationAngleX = -1.820989;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 898 based Row/Col = 89/8
      geom_tower geom;
      geom.id = 898;
      geom.pDz = 6.751955;
      geom.pTheta = 0.083707;
      geom.pPhi = 0.127510;
      geom.pAlp1 = -0.021345;
      geom.pAlp2 = -0.021342;
      geom.pDy1 = 1.126302;
      geom.pDx1 = 1.197557;
      geom.pDx2 = 1.190598;
      geom.pDy2 = 1.270041;
      geom.pDx3 = 1.359272;
      geom.pDx4 = 1.351425;
      geom.centralX = 8.900303;
      geom.centralY = 104.430798;
      geom.centralZ = -27.928983;
      geom.pRotationAngleX = -1.820989;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 901 based Row/Col = 90/1
      geom_tower geom;
      geom.id = 901;
      geom.pDz = 6.751991;
      geom.pTheta = 0.084148;
      geom.pPhi = -3.015360;
      geom.pAlp1 = 0.021345;
      geom.pAlp2 = 0.021342;
      geom.pDy1 = 1.126302;
      geom.pDx1 = 1.204518;
      geom.pDx2 = 1.197557;
      geom.pDy2 = 1.270041;
      geom.pDx3 = 1.367120;
      geom.pDx4 = 1.359272;
      geom.centralX = -8.951457;
      geom.centralY = 105.025349;
      geom.centralZ = -25.602406;
      geom.pRotationAngleX = -1.820989;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 902 based Row/Col = 90/2
      geom_tower geom;
      geom.id = 902;
      geom.pDz = 6.751991;
      geom.pTheta = 0.060563;
      geom.pPhi = -2.965561;
      geom.pAlp1 = 0.015197;
      geom.pAlp2 = 0.015195;
      geom.pDy1 = 1.126302;
      geom.pDx1 = 1.200376;
      geom.pDx2 = 1.193483;
      geom.pDy2 = 1.270041;
      geom.pDx3 = 1.362420;
      geom.pDx4 = 1.354649;
      geom.centralX = -6.386608;
      geom.centralY = 105.025349;
      geom.centralZ = -25.602406;
      geom.pRotationAngleX = -1.820989;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 903 based Row/Col = 90/3
      geom_tower geom;
      geom.id = 903;
      geom.pDz = 6.751991;
      geom.pTheta = 0.037318;
      geom.pPhi = -2.853184;
      geom.pAlp1 = 0.009099;
      geom.pAlp2 = 0.009097;
      geom.pDy1 = 1.126302;
      geom.pDx1 = 1.197625;
      geom.pDx2 = 1.190778;
      geom.pDy2 = 1.270041;
      geom.pDx3 = 1.359299;
      geom.pDx4 = 1.351579;
      geom.centralX = -3.829055;
      geom.centralY = 105.025349;
      geom.centralZ = -25.602406;
      geom.pRotationAngleX = -1.820989;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 904 based Row/Col = 90/4
      geom_tower geom;
      geom.id = 904;
      geom.pDz = 6.751991;
      geom.pTheta = 0.015967;
      geom.pPhi = -2.414118;
      geom.pAlp1 = 0.003030;
      geom.pAlp2 = 0.003029;
      geom.pDy1 = 1.126302;
      geom.pDx1 = 1.196253;
      geom.pDx2 = 1.189428;
      geom.pDy2 = 1.270041;
      geom.pDx3 = 1.357742;
      geom.pDx4 = 1.350048;
      geom.centralX = -1.275868;
      geom.centralY = 105.025349;
      geom.centralZ = -25.602406;
      geom.pRotationAngleX = -1.820989;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 905 based Row/Col = 90/5
      geom_tower geom;
      geom.id = 905;
      geom.pDz = 6.751991;
      geom.pTheta = 0.015967;
      geom.pPhi = -0.727474;
      geom.pAlp1 = -0.003030;
      geom.pAlp2 = -0.003029;
      geom.pDy1 = 1.126302;
      geom.pDx1 = 1.196253;
      geom.pDx2 = 1.189428;
      geom.pDy2 = 1.270041;
      geom.pDx3 = 1.357742;
      geom.pDx4 = 1.350048;
      geom.centralX = 1.275868;
      geom.centralY = 105.025349;
      geom.centralZ = -25.602406;
      geom.pRotationAngleX = -1.820989;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 906 based Row/Col = 90/6
      geom_tower geom;
      geom.id = 906;
      geom.pDz = 6.751991;
      geom.pTheta = 0.037318;
      geom.pPhi = -0.288409;
      geom.pAlp1 = -0.009099;
      geom.pAlp2 = -0.009097;
      geom.pDy1 = 1.126302;
      geom.pDx1 = 1.197625;
      geom.pDx2 = 1.190778;
      geom.pDy2 = 1.270041;
      geom.pDx3 = 1.359299;
      geom.pDx4 = 1.351579;
      geom.centralX = 3.829055;
      geom.centralY = 105.025349;
      geom.centralZ = -25.602406;
      geom.pRotationAngleX = -1.820989;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 907 based Row/Col = 90/7
      geom_tower geom;
      geom.id = 907;
      geom.pDz = 6.751991;
      geom.pTheta = 0.060563;
      geom.pPhi = -0.176032;
      geom.pAlp1 = -0.015197;
      geom.pAlp2 = -0.015195;
      geom.pDy1 = 1.126302;
      geom.pDx1 = 1.200376;
      geom.pDx2 = 1.193483;
      geom.pDy2 = 1.270041;
      geom.pDx3 = 1.362420;
      geom.pDx4 = 1.354649;
      geom.centralX = 6.386608;
      geom.centralY = 105.025349;
      geom.centralZ = -25.602406;
      geom.pRotationAngleX = -1.820989;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 908 based Row/Col = 90/8
      geom_tower geom;
      geom.id = 908;
      geom.pDz = 6.751991;
      geom.pTheta = 0.084148;
      geom.pPhi = -0.126233;
      geom.pAlp1 = -0.021345;
      geom.pAlp2 = -0.021342;
      geom.pDy1 = 1.126302;
      geom.pDx1 = 1.204518;
      geom.pDx2 = 1.197557;
      geom.pDy2 = 1.270041;
      geom.pDx3 = 1.367120;
      geom.pDx4 = 1.359272;
      geom.centralX = 8.951457;
      geom.centralY = 105.025349;
      geom.centralZ = -25.602406;
      geom.pRotationAngleX = -1.820989;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1161 based Row/Col = 116/1
      geom_tower geom;
      geom.id = 1161;
      geom.pDz = 6.751923;
      geom.pTheta = 0.081450;
      geom.pPhi = -3.013805;
      geom.pAlp1 = -0.028586;
      geom.pAlp2 = -0.028578;
      geom.pDy1 = 1.127255;
      geom.pDx1 = 1.189615;
      geom.pDx2 = 1.198938;
      geom.pDy2 = 1.267455;
      geom.pDx3 = 1.345894;
      geom.pDx4 = 1.356375;
      geom.centralX = -8.888448;
      geom.centralY = 104.290051;
      geom.centralZ = 37.931439;
      geom.pRotationAngleX = -1.232836;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1162 based Row/Col = 116/2
      geom_tower geom;
      geom.id = 1162;
      geom.pDz = 6.751923;
      geom.pTheta = 0.058632;
      geom.pPhi = -2.963426;
      geom.pAlp1 = -0.020358;
      geom.pAlp2 = -0.020353;
      geom.pDy1 = 1.127255;
      geom.pDx1 = 1.185832;
      geom.pDx2 = 1.195070;
      geom.pDy2 = 1.267455;
      geom.pDx3 = 1.341616;
      geom.pDx4 = 1.351999;
      geom.centralX = -6.342113;
      geom.centralY = 104.290051;
      geom.centralZ = 37.931439;
      geom.pRotationAngleX = -1.232836;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1163 based Row/Col = 116/3
      geom_tower geom;
      geom.id = 1163;
      geom.pDz = 6.751923;
      geom.pTheta = 0.036150;
      geom.pPhi = -2.849824;
      geom.pAlp1 = -0.012191;
      geom.pAlp2 = -0.012187;
      geom.pDy1 = 1.127255;
      geom.pDx1 = 1.183319;
      geom.pDx2 = 1.192500;
      geom.pDy2 = 1.267455;
      geom.pDx3 = 1.338773;
      geom.pDx4 = 1.349093;
      geom.centralX = -3.802562;
      geom.centralY = 104.290051;
      geom.centralZ = 37.931439;
      geom.pRotationAngleX = -1.232836;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1164 based Row/Col = 116/4
      geom_tower geom;
      geom.id = 1164;
      geom.pDz = 6.751923;
      geom.pTheta = 0.015537;
      geom.pPhi = -2.408038;
      geom.pAlp1 = -0.004060;
      geom.pAlp2 = -0.004058;
      geom.pDy1 = 1.127255;
      geom.pDx1 = 1.182065;
      geom.pDx2 = 1.191218;
      geom.pDy2 = 1.267455;
      geom.pDx3 = 1.337355;
      geom.pDx4 = 1.347643;
      geom.centralX = -1.267070;
      geom.centralY = 104.290051;
      geom.centralZ = 37.931439;
      geom.pRotationAngleX = -1.232836;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1165 based Row/Col = 116/5
      geom_tower geom;
      geom.id = 1165;
      geom.pDz = 6.751923;
      geom.pTheta = 0.015537;
      geom.pPhi = -0.733555;
      geom.pAlp1 = 0.004060;
      geom.pAlp2 = 0.004058;
      geom.pDy1 = 1.127255;
      geom.pDx1 = 1.182065;
      geom.pDx2 = 1.191218;
      geom.pDy2 = 1.267455;
      geom.pDx3 = 1.337355;
      geom.pDx4 = 1.347643;
      geom.centralX = 1.267070;
      geom.centralY = 104.290051;
      geom.centralZ = 37.931439;
      geom.pRotationAngleX = -1.232836;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1166 based Row/Col = 116/6
      geom_tower geom;
      geom.id = 1166;
      geom.pDz = 6.751923;
      geom.pTheta = 0.036150;
      geom.pPhi = -0.291769;
      geom.pAlp1 = 0.012191;
      geom.pAlp2 = 0.012187;
      geom.pDy1 = 1.127255;
      geom.pDx1 = 1.183319;
      geom.pDx2 = 1.192500;
      geom.pDy2 = 1.267455;
      geom.pDx3 = 1.338773;
      geom.pDx4 = 1.349093;
      geom.centralX = 3.802562;
      geom.centralY = 104.290051;
      geom.centralZ = 37.931439;
      geom.pRotationAngleX = -1.232836;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1167 based Row/Col = 116/7
      geom_tower geom;
      geom.id = 1167;
      geom.pDz = 6.751923;
      geom.pTheta = 0.058632;
      geom.pPhi = -0.178166;
      geom.pAlp1 = 0.020358;
      geom.pAlp2 = 0.020353;
      geom.pDy1 = 1.127255;
      geom.pDx1 = 1.185832;
      geom.pDx2 = 1.195070;
      geom.pDy2 = 1.267455;
      geom.pDx3 = 1.341616;
      geom.pDx4 = 1.351999;
      geom.centralX = 6.342113;
      geom.centralY = 104.290051;
      geom.centralZ = 37.931439;
      geom.pRotationAngleX = -1.232836;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1168 based Row/Col = 116/8
      geom_tower geom;
      geom.id = 1168;
      geom.pDz = 6.751923;
      geom.pTheta = 0.081450;
      geom.pPhi = -0.127788;
      geom.pAlp1 = 0.028586;
      geom.pAlp2 = 0.028578;
      geom.pDy1 = 1.127255;
      geom.pDx1 = 1.189615;
      geom.pDx2 = 1.198938;
      geom.pDy2 = 1.267455;
      geom.pDx3 = 1.345894;
      geom.pDx4 = 1.356375;
      geom.centralX = 8.888448;
      geom.centralY = 104.290051;
      geom.centralZ = 37.931439;
      geom.pRotationAngleX = -1.232836;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1151 based Row/Col = 115/1
      geom_tower geom;
      geom.id = 1151;
      geom.pDz = 6.751919;
      geom.pTheta = 0.082028;
      geom.pPhi = 3.015224;
      geom.pAlp1 = -0.028585;
      geom.pAlp2 = -0.028578;
      geom.pDy1 = 1.127255;
      geom.pDx1 = 1.198938;
      geom.pDx2 = 1.208263;
      geom.pDy2 = 1.267455;
      geom.pDx3 = 1.356375;
      geom.pDx4 = 1.366857;
      geom.centralX = -8.956910;
      geom.centralY = 105.085706;
      geom.centralZ = 35.667474;
      geom.pRotationAngleX = -1.232836;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1152 based Row/Col = 115/2
      geom_tower geom;
      geom.id = 1152;
      geom.pDz = 6.751919;
      geom.pTheta = 0.059038;
      geom.pPhi = 2.965382;
      geom.pAlp1 = -0.020357;
      geom.pAlp2 = -0.020351;
      geom.pDy1 = 1.127255;
      geom.pDx1 = 1.195070;
      geom.pDx2 = 1.204307;
      geom.pDy2 = 1.267455;
      geom.pDx3 = 1.351999;
      geom.pDx4 = 1.362383;
      geom.centralX = -6.390862;
      geom.centralY = 105.085706;
      geom.centralZ = 35.667474;
      geom.pRotationAngleX = -1.232836;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1153 based Row/Col = 115/3
      geom_tower geom;
      geom.id = 1153;
      geom.pDz = 6.751919;
      geom.pTheta = 0.036380;
      geom.pPhi = 2.852912;
      geom.pAlp1 = -0.012190;
      geom.pAlp2 = -0.012186;
      geom.pDy1 = 1.127255;
      geom.pDx1 = 1.192500;
      geom.pDx2 = 1.201680;
      geom.pDy2 = 1.267455;
      geom.pDx3 = 1.349093;
      geom.pDx4 = 1.359412;
      geom.centralX = -3.831751;
      geom.centralY = 105.085706;
      geom.centralZ = 35.667474;
      geom.pRotationAngleX = -1.232836;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1154 based Row/Col = 115/4
      geom_tower geom;
      geom.id = 1154;
      geom.pDz = 6.751919;
      geom.pTheta = 0.015572;
      geom.pPhi = 2.413632;
      geom.pAlp1 = -0.004059;
      geom.pAlp2 = -0.004058;
      geom.pDy1 = 1.127255;
      geom.pDx1 = 1.191218;
      geom.pDx2 = 1.200369;
      geom.pDy2 = 1.267455;
      geom.pDx3 = 1.347643;
      geom.pDx4 = 1.357930;
      geom.centralX = -1.276790;
      geom.centralY = 105.085706;
      geom.centralZ = 35.667474;
      geom.pRotationAngleX = -1.232836;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1155 based Row/Col = 115/5
      geom_tower geom;
      geom.id = 1155;
      geom.pDz = 6.751919;
      geom.pTheta = 0.015572;
      geom.pPhi = 0.727960;
      geom.pAlp1 = 0.004059;
      geom.pAlp2 = 0.004058;
      geom.pDy1 = 1.127255;
      geom.pDx1 = 1.191218;
      geom.pDx2 = 1.200369;
      geom.pDy2 = 1.267455;
      geom.pDx3 = 1.347643;
      geom.pDx4 = 1.357930;
      geom.centralX = 1.276790;
      geom.centralY = 105.085706;
      geom.centralZ = 35.667474;
      geom.pRotationAngleX = -1.232836;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1156 based Row/Col = 115/6
      geom_tower geom;
      geom.id = 1156;
      geom.pDz = 6.751919;
      geom.pTheta = 0.036380;
      geom.pPhi = 0.288681;
      geom.pAlp1 = 0.012190;
      geom.pAlp2 = 0.012186;
      geom.pDy1 = 1.127255;
      geom.pDx1 = 1.192500;
      geom.pDx2 = 1.201680;
      geom.pDy2 = 1.267455;
      geom.pDx3 = 1.349093;
      geom.pDx4 = 1.359412;
      geom.centralX = 3.831751;
      geom.centralY = 105.085706;
      geom.centralZ = 35.667474;
      geom.pRotationAngleX = -1.232836;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1157 based Row/Col = 115/7
      geom_tower geom;
      geom.id = 1157;
      geom.pDz = 6.751919;
      geom.pTheta = 0.059038;
      geom.pPhi = 0.176211;
      geom.pAlp1 = 0.020357;
      geom.pAlp2 = 0.020351;
      geom.pDy1 = 1.127255;
      geom.pDx1 = 1.195070;
      geom.pDx2 = 1.204307;
      geom.pDy2 = 1.267455;
      geom.pDx3 = 1.351999;
      geom.pDx4 = 1.362383;
      geom.centralX = 6.390862;
      geom.centralY = 105.085706;
      geom.centralZ = 35.667474;
      geom.pRotationAngleX = -1.232836;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1158 based Row/Col = 115/8
      geom_tower geom;
      geom.id = 1158;
      geom.pDz = 6.751919;
      geom.pTheta = 0.082028;
      geom.pPhi = 0.126369;
      geom.pAlp1 = 0.028585;
      geom.pAlp2 = 0.028578;
      geom.pDy1 = 1.127255;
      geom.pDx1 = 1.198938;
      geom.pDx2 = 1.208263;
      geom.pDy2 = 1.267455;
      geom.pDx3 = 1.356375;
      geom.pDx4 = 1.366857;
      geom.centralX = 8.956910;
      geom.centralY = 105.085706;
      geom.centralZ = 35.667474;
      geom.pRotationAngleX = -1.232836;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 851 based Row/Col = 85/1
      geom_tower geom;
      geom.id = 851;
      geom.pDz = 6.751923;
      geom.pTheta = 0.081450;
      geom.pPhi = 3.013805;
      geom.pAlp1 = 0.028586;
      geom.pAlp2 = 0.028578;
      geom.pDy1 = 1.127255;
      geom.pDx1 = 1.198938;
      geom.pDx2 = 1.189615;
      geom.pDy2 = 1.267455;
      geom.pDx3 = 1.356375;
      geom.pDx4 = 1.345894;
      geom.centralX = -8.888448;
      geom.centralY = 104.290051;
      geom.centralZ = -37.931439;
      geom.pRotationAngleX = -1.908756;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 852 based Row/Col = 85/2
      geom_tower geom;
      geom.id = 852;
      geom.pDz = 6.751923;
      geom.pTheta = 0.058632;
      geom.pPhi = 2.963426;
      geom.pAlp1 = 0.020358;
      geom.pAlp2 = 0.020353;
      geom.pDy1 = 1.127255;
      geom.pDx1 = 1.195070;
      geom.pDx2 = 1.185832;
      geom.pDy2 = 1.267455;
      geom.pDx3 = 1.351999;
      geom.pDx4 = 1.341616;
      geom.centralX = -6.342113;
      geom.centralY = 104.290051;
      geom.centralZ = -37.931439;
      geom.pRotationAngleX = -1.908756;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 853 based Row/Col = 85/3
      geom_tower geom;
      geom.id = 853;
      geom.pDz = 6.751923;
      geom.pTheta = 0.036150;
      geom.pPhi = 2.849824;
      geom.pAlp1 = 0.012191;
      geom.pAlp2 = 0.012187;
      geom.pDy1 = 1.127255;
      geom.pDx1 = 1.192500;
      geom.pDx2 = 1.183319;
      geom.pDy2 = 1.267455;
      geom.pDx3 = 1.349093;
      geom.pDx4 = 1.338773;
      geom.centralX = -3.802562;
      geom.centralY = 104.290051;
      geom.centralZ = -37.931439;
      geom.pRotationAngleX = -1.908756;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 854 based Row/Col = 85/4
      geom_tower geom;
      geom.id = 854;
      geom.pDz = 6.751923;
      geom.pTheta = 0.015537;
      geom.pPhi = 2.408038;
      geom.pAlp1 = 0.004060;
      geom.pAlp2 = 0.004058;
      geom.pDy1 = 1.127255;
      geom.pDx1 = 1.191218;
      geom.pDx2 = 1.182065;
      geom.pDy2 = 1.267455;
      geom.pDx3 = 1.347643;
      geom.pDx4 = 1.337355;
      geom.centralX = -1.267070;
      geom.centralY = 104.290051;
      geom.centralZ = -37.931439;
      geom.pRotationAngleX = -1.908756;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 855 based Row/Col = 85/5
      geom_tower geom;
      geom.id = 855;
      geom.pDz = 6.751923;
      geom.pTheta = 0.015537;
      geom.pPhi = 0.733555;
      geom.pAlp1 = -0.004060;
      geom.pAlp2 = -0.004058;
      geom.pDy1 = 1.127255;
      geom.pDx1 = 1.191218;
      geom.pDx2 = 1.182065;
      geom.pDy2 = 1.267455;
      geom.pDx3 = 1.347643;
      geom.pDx4 = 1.337355;
      geom.centralX = 1.267070;
      geom.centralY = 104.290051;
      geom.centralZ = -37.931439;
      geom.pRotationAngleX = -1.908756;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 856 based Row/Col = 85/6
      geom_tower geom;
      geom.id = 856;
      geom.pDz = 6.751923;
      geom.pTheta = 0.036150;
      geom.pPhi = 0.291769;
      geom.pAlp1 = -0.012191;
      geom.pAlp2 = -0.012187;
      geom.pDy1 = 1.127255;
      geom.pDx1 = 1.192500;
      geom.pDx2 = 1.183319;
      geom.pDy2 = 1.267455;
      geom.pDx3 = 1.349093;
      geom.pDx4 = 1.338773;
      geom.centralX = 3.802562;
      geom.centralY = 104.290051;
      geom.centralZ = -37.931439;
      geom.pRotationAngleX = -1.908756;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 857 based Row/Col = 85/7
      geom_tower geom;
      geom.id = 857;
      geom.pDz = 6.751923;
      geom.pTheta = 0.058632;
      geom.pPhi = 0.178166;
      geom.pAlp1 = -0.020358;
      geom.pAlp2 = -0.020353;
      geom.pDy1 = 1.127255;
      geom.pDx1 = 1.195070;
      geom.pDx2 = 1.185832;
      geom.pDy2 = 1.267455;
      geom.pDx3 = 1.351999;
      geom.pDx4 = 1.341616;
      geom.centralX = 6.342113;
      geom.centralY = 104.290051;
      geom.centralZ = -37.931439;
      geom.pRotationAngleX = -1.908756;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 858 based Row/Col = 85/8
      geom_tower geom;
      geom.id = 858;
      geom.pDz = 6.751923;
      geom.pTheta = 0.081450;
      geom.pPhi = 0.127788;
      geom.pAlp1 = -0.028586;
      geom.pAlp2 = -0.028578;
      geom.pDy1 = 1.127255;
      geom.pDx1 = 1.198938;
      geom.pDx2 = 1.189615;
      geom.pDy2 = 1.267455;
      geom.pDx3 = 1.356375;
      geom.pDx4 = 1.345894;
      geom.centralX = 8.888448;
      geom.centralY = 104.290051;
      geom.centralZ = -37.931439;
      geom.pRotationAngleX = -1.908756;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 861 based Row/Col = 86/1
      geom_tower geom;
      geom.id = 861;
      geom.pDz = 6.751919;
      geom.pTheta = 0.082028;
      geom.pPhi = -3.015224;
      geom.pAlp1 = 0.028585;
      geom.pAlp2 = 0.028578;
      geom.pDy1 = 1.127255;
      geom.pDx1 = 1.208263;
      geom.pDx2 = 1.198938;
      geom.pDy2 = 1.267455;
      geom.pDx3 = 1.366857;
      geom.pDx4 = 1.356375;
      geom.centralX = -8.956910;
      geom.centralY = 105.085706;
      geom.centralZ = -35.667474;
      geom.pRotationAngleX = -1.908756;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 862 based Row/Col = 86/2
      geom_tower geom;
      geom.id = 862;
      geom.pDz = 6.751919;
      geom.pTheta = 0.059038;
      geom.pPhi = -2.965382;
      geom.pAlp1 = 0.020357;
      geom.pAlp2 = 0.020351;
      geom.pDy1 = 1.127255;
      geom.pDx1 = 1.204307;
      geom.pDx2 = 1.195070;
      geom.pDy2 = 1.267455;
      geom.pDx3 = 1.362383;
      geom.pDx4 = 1.351999;
      geom.centralX = -6.390862;
      geom.centralY = 105.085706;
      geom.centralZ = -35.667474;
      geom.pRotationAngleX = -1.908756;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 863 based Row/Col = 86/3
      geom_tower geom;
      geom.id = 863;
      geom.pDz = 6.751919;
      geom.pTheta = 0.036380;
      geom.pPhi = -2.852912;
      geom.pAlp1 = 0.012190;
      geom.pAlp2 = 0.012186;
      geom.pDy1 = 1.127255;
      geom.pDx1 = 1.201680;
      geom.pDx2 = 1.192500;
      geom.pDy2 = 1.267455;
      geom.pDx3 = 1.359412;
      geom.pDx4 = 1.349093;
      geom.centralX = -3.831751;
      geom.centralY = 105.085706;
      geom.centralZ = -35.667474;
      geom.pRotationAngleX = -1.908756;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 864 based Row/Col = 86/4
      geom_tower geom;
      geom.id = 864;
      geom.pDz = 6.751919;
      geom.pTheta = 0.015572;
      geom.pPhi = -2.413632;
      geom.pAlp1 = 0.004059;
      geom.pAlp2 = 0.004058;
      geom.pDy1 = 1.127255;
      geom.pDx1 = 1.200369;
      geom.pDx2 = 1.191218;
      geom.pDy2 = 1.267455;
      geom.pDx3 = 1.357930;
      geom.pDx4 = 1.347643;
      geom.centralX = -1.276790;
      geom.centralY = 105.085706;
      geom.centralZ = -35.667474;
      geom.pRotationAngleX = -1.908756;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 865 based Row/Col = 86/5
      geom_tower geom;
      geom.id = 865;
      geom.pDz = 6.751919;
      geom.pTheta = 0.015572;
      geom.pPhi = -0.727960;
      geom.pAlp1 = -0.004059;
      geom.pAlp2 = -0.004058;
      geom.pDy1 = 1.127255;
      geom.pDx1 = 1.200369;
      geom.pDx2 = 1.191218;
      geom.pDy2 = 1.267455;
      geom.pDx3 = 1.357930;
      geom.pDx4 = 1.347643;
      geom.centralX = 1.276790;
      geom.centralY = 105.085706;
      geom.centralZ = -35.667474;
      geom.pRotationAngleX = -1.908756;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 866 based Row/Col = 86/6
      geom_tower geom;
      geom.id = 866;
      geom.pDz = 6.751919;
      geom.pTheta = 0.036380;
      geom.pPhi = -0.288681;
      geom.pAlp1 = -0.012190;
      geom.pAlp2 = -0.012186;
      geom.pDy1 = 1.127255;
      geom.pDx1 = 1.201680;
      geom.pDx2 = 1.192500;
      geom.pDy2 = 1.267455;
      geom.pDx3 = 1.359412;
      geom.pDx4 = 1.349093;
      geom.centralX = 3.831751;
      geom.centralY = 105.085706;
      geom.centralZ = -35.667474;
      geom.pRotationAngleX = -1.908756;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 867 based Row/Col = 86/7
      geom_tower geom;
      geom.id = 867;
      geom.pDz = 6.751919;
      geom.pTheta = 0.059038;
      geom.pPhi = -0.176211;
      geom.pAlp1 = -0.020357;
      geom.pAlp2 = -0.020351;
      geom.pDy1 = 1.127255;
      geom.pDx1 = 1.204307;
      geom.pDx2 = 1.195070;
      geom.pDy2 = 1.267455;
      geom.pDx3 = 1.362383;
      geom.pDx4 = 1.351999;
      geom.centralX = 6.390862;
      geom.centralY = 105.085706;
      geom.centralZ = -35.667474;
      geom.pRotationAngleX = -1.908756;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 868 based Row/Col = 86/8
      geom_tower geom;
      geom.id = 868;
      geom.pDz = 6.751919;
      geom.pTheta = 0.082028;
      geom.pPhi = -0.126369;
      geom.pAlp1 = -0.028585;
      geom.pAlp2 = -0.028578;
      geom.pDy1 = 1.127255;
      geom.pDx1 = 1.208263;
      geom.pDx2 = 1.198938;
      geom.pDy2 = 1.267455;
      geom.pDx3 = 1.366857;
      geom.pDx4 = 1.356375;
      geom.centralX = 8.956910;
      geom.centralY = 105.085706;
      geom.centralZ = -35.667474;
      geom.pRotationAngleX = -1.908756;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1201 based Row/Col = 120/1
      geom_tower geom;
      geom.id = 1201;
      geom.pDz = 6.751937;
      geom.pTheta = 0.078667;
      geom.pPhi = -3.014605;
      geom.pAlp1 = -0.035387;
      geom.pAlp2 = -0.035372;
      geom.pDy1 = 1.128053;
      geom.pDx1 = 1.189943;
      geom.pDx2 = 1.201484;
      geom.pDy2 = 1.262406;
      geom.pDx3 = 1.340707;
      geom.pDx4 = 1.353618;
      geom.centralX = -8.882151;
      geom.centralY = 104.213247;
      geom.centralZ = 48.242347;
      geom.pRotationAngleX = -1.147870;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1202 based Row/Col = 120/2
      geom_tower geom;
      geom.id = 1202;
      geom.pDz = 6.751937;
      geom.pTheta = 0.056623;
      geom.pPhi = -2.964544;
      geom.pAlp1 = -0.025209;
      geom.pAlp2 = -0.025198;
      geom.pDy1 = 1.128053;
      geom.pDx1 = 1.186419;
      geom.pDx2 = 1.197861;
      geom.pDy2 = 1.262406;
      geom.pDx3 = 1.336737;
      geom.pDx4 = 1.349537;
      geom.centralX = -6.338074;
      geom.centralY = 104.213247;
      geom.centralZ = 48.242347;
      geom.pRotationAngleX = -1.147870;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1203 based Row/Col = 120/3
      geom_tower geom;
      geom.id = 1203;
      geom.pDz = 6.751937;
      geom.pTheta = 0.034900;
      geom.pPhi = -2.851602;
      geom.pAlp1 = -0.015099;
      geom.pAlp2 = -0.015092;
      geom.pDy1 = 1.128053;
      geom.pDx1 = 1.184077;
      geom.pDx2 = 1.195454;
      geom.pDy2 = 1.262406;
      geom.pDx3 = 1.334099;
      geom.pDx4 = 1.346826;
      geom.centralX = -3.800322;
      geom.centralY = 104.213247;
      geom.centralZ = 48.242347;
      geom.pRotationAngleX = -1.147870;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1204 based Row/Col = 120/4
      geom_tower geom;
      geom.id = 1204;
      geom.pDz = 6.751937;
      geom.pTheta = 0.014964;
      geom.pPhi = -2.411267;
      geom.pAlp1 = -0.005028;
      geom.pAlp2 = -0.005026;
      geom.pDy1 = 1.128053;
      geom.pDx1 = 1.182908;
      geom.pDx2 = 1.194253;
      geom.pDy2 = 1.262406;
      geom.pDx3 = 1.332782;
      geom.pDx4 = 1.345473;
      geom.centralX = -1.266354;
      geom.centralY = 104.213247;
      geom.centralZ = 48.242347;
      geom.pRotationAngleX = -1.147870;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1205 based Row/Col = 120/5
      geom_tower geom;
      geom.id = 1205;
      geom.pDz = 6.751937;
      geom.pTheta = 0.014964;
      geom.pPhi = -0.730325;
      geom.pAlp1 = 0.005028;
      geom.pAlp2 = 0.005026;
      geom.pDy1 = 1.128053;
      geom.pDx1 = 1.182908;
      geom.pDx2 = 1.194253;
      geom.pDy2 = 1.262406;
      geom.pDx3 = 1.332782;
      geom.pDx4 = 1.345473;
      geom.centralX = 1.266354;
      geom.centralY = 104.213247;
      geom.centralZ = 48.242347;
      geom.pRotationAngleX = -1.147870;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1206 based Row/Col = 120/6
      geom_tower geom;
      geom.id = 1206;
      geom.pDz = 6.751937;
      geom.pTheta = 0.034900;
      geom.pPhi = -0.289990;
      geom.pAlp1 = 0.015099;
      geom.pAlp2 = 0.015092;
      geom.pDy1 = 1.128053;
      geom.pDx1 = 1.184077;
      geom.pDx2 = 1.195454;
      geom.pDy2 = 1.262406;
      geom.pDx3 = 1.334099;
      geom.pDx4 = 1.346826;
      geom.centralX = 3.800322;
      geom.centralY = 104.213247;
      geom.centralZ = 48.242347;
      geom.pRotationAngleX = -1.147870;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1207 based Row/Col = 120/7
      geom_tower geom;
      geom.id = 1207;
      geom.pDz = 6.751937;
      geom.pTheta = 0.056623;
      geom.pPhi = -0.177049;
      geom.pAlp1 = 0.025209;
      geom.pAlp2 = 0.025198;
      geom.pDy1 = 1.128053;
      geom.pDx1 = 1.186419;
      geom.pDx2 = 1.197861;
      geom.pDy2 = 1.262406;
      geom.pDx3 = 1.336737;
      geom.pDx4 = 1.349537;
      geom.centralX = 6.338074;
      geom.centralY = 104.213247;
      geom.centralZ = 48.242347;
      geom.pRotationAngleX = -1.147870;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1208 based Row/Col = 120/8
      geom_tower geom;
      geom.id = 1208;
      geom.pDz = 6.751937;
      geom.pTheta = 0.078667;
      geom.pPhi = -0.126987;
      geom.pAlp1 = 0.035387;
      geom.pAlp2 = 0.035372;
      geom.pDy1 = 1.128053;
      geom.pDx1 = 1.189943;
      geom.pDx2 = 1.201484;
      geom.pDy2 = 1.262406;
      geom.pDx3 = 1.340707;
      geom.pDx4 = 1.353618;
      geom.centralX = 8.882151;
      geom.centralY = 104.213247;
      geom.centralZ = 48.242347;
      geom.pRotationAngleX = -1.147870;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1191 based Row/Col = 119/1
      geom_tower geom;
      geom.id = 1191;
      geom.pDz = 6.751840;
      geom.pTheta = 0.079351;
      geom.pPhi = 3.016577;
      geom.pAlp1 = -0.035386;
      geom.pAlp2 = -0.035371;
      geom.pDy1 = 1.128053;
      geom.pDx1 = 1.201484;
      geom.pDx2 = 1.213027;
      geom.pDy2 = 1.262406;
      geom.pDx3 = 1.353618;
      geom.pDx4 = 1.366531;
      geom.centralX = -8.966758;
      geom.centralY = 105.196418;
      geom.centralZ = 46.057947;
      geom.pRotationAngleX = -1.147870;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1192 based Row/Col = 119/2
      geom_tower geom;
      geom.id = 1192;
      geom.pDz = 6.751840;
      geom.pTheta = 0.057103;
      geom.pPhi = 2.967263;
      geom.pAlp1 = -0.025207;
      geom.pAlp2 = -0.025196;
      geom.pDy1 = 1.128053;
      geom.pDx1 = 1.197861;
      geom.pDx2 = 1.209304;
      geom.pDy2 = 1.262406;
      geom.pDx3 = 1.349537;
      geom.pDx4 = 1.362338;
      geom.centralX = -6.398332;
      geom.centralY = 105.196418;
      geom.centralZ = 46.057947;
      geom.pRotationAngleX = -1.147870;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1193 based Row/Col = 119/3
      geom_tower geom;
      geom.id = 1193;
      geom.pDz = 6.751840;
      geom.pTheta = 0.035168;
      geom.pPhi = 2.855900;
      geom.pAlp1 = -0.015097;
      geom.pAlp2 = -0.015090;
      geom.pDy1 = 1.128053;
      geom.pDx1 = 1.195454;
      geom.pDx2 = 1.206830;
      geom.pDy2 = 1.262406;
      geom.pDx3 = 1.346826;
      geom.pDx4 = 1.359552;
      geom.centralX = -3.836407;
      geom.centralY = 105.196418;
      geom.centralZ = 46.057947;
      geom.pRotationAngleX = -1.147870;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1194 based Row/Col = 119/4
      geom_tower geom;
      geom.id = 1194;
      geom.pDz = 6.751840;
      geom.pTheta = 0.014993;
      geom.pPhi = 2.419104;
      geom.pAlp1 = -0.005028;
      geom.pAlp2 = -0.005026;
      geom.pDy1 = 1.128053;
      geom.pDx1 = 1.194253;
      geom.pDx2 = 1.205596;
      geom.pDy2 = 1.262406;
      geom.pDx3 = 1.345473;
      geom.pDx4 = 1.358161;
      geom.centralX = -1.278371;
      geom.centralY = 105.196418;
      geom.centralZ = 46.057947;
      geom.pRotationAngleX = -1.147870;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1195 based Row/Col = 119/5
      geom_tower geom;
      geom.id = 1195;
      geom.pDz = 6.751840;
      geom.pTheta = 0.014993;
      geom.pPhi = 0.722489;
      geom.pAlp1 = 0.005028;
      geom.pAlp2 = 0.005026;
      geom.pDy1 = 1.128053;
      geom.pDx1 = 1.194253;
      geom.pDx2 = 1.205596;
      geom.pDy2 = 1.262406;
      geom.pDx3 = 1.345473;
      geom.pDx4 = 1.358161;
      geom.centralX = 1.278371;
      geom.centralY = 105.196418;
      geom.centralZ = 46.057947;
      geom.pRotationAngleX = -1.147870;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1196 based Row/Col = 119/6
      geom_tower geom;
      geom.id = 1196;
      geom.pDz = 6.751840;
      geom.pTheta = 0.035168;
      geom.pPhi = 0.285692;
      geom.pAlp1 = 0.015097;
      geom.pAlp2 = 0.015090;
      geom.pDy1 = 1.128053;
      geom.pDx1 = 1.195454;
      geom.pDx2 = 1.206830;
      geom.pDy2 = 1.262406;
      geom.pDx3 = 1.346826;
      geom.pDx4 = 1.359552;
      geom.centralX = 3.836407;
      geom.centralY = 105.196418;
      geom.centralZ = 46.057947;
      geom.pRotationAngleX = -1.147870;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1197 based Row/Col = 119/7
      geom_tower geom;
      geom.id = 1197;
      geom.pDz = 6.751840;
      geom.pTheta = 0.057103;
      geom.pPhi = 0.174330;
      geom.pAlp1 = 0.025207;
      geom.pAlp2 = 0.025196;
      geom.pDy1 = 1.128053;
      geom.pDx1 = 1.197861;
      geom.pDx2 = 1.209304;
      geom.pDy2 = 1.262406;
      geom.pDx3 = 1.349537;
      geom.pDx4 = 1.362338;
      geom.centralX = 6.398332;
      geom.centralY = 105.196418;
      geom.centralZ = 46.057947;
      geom.pRotationAngleX = -1.147870;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1198 based Row/Col = 119/8
      geom_tower geom;
      geom.id = 1198;
      geom.pDz = 6.751840;
      geom.pTheta = 0.079351;
      geom.pPhi = 0.125015;
      geom.pAlp1 = 0.035386;
      geom.pAlp2 = 0.035371;
      geom.pDy1 = 1.128053;
      geom.pDx1 = 1.201484;
      geom.pDx2 = 1.213027;
      geom.pDy2 = 1.262406;
      geom.pDx3 = 1.353618;
      geom.pDx4 = 1.366531;
      geom.centralX = 8.966758;
      geom.centralY = 105.196418;
      geom.centralZ = 46.057947;
      geom.pRotationAngleX = -1.147870;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 811 based Row/Col = 81/1
      geom_tower geom;
      geom.id = 811;
      geom.pDz = 6.751937;
      geom.pTheta = 0.078667;
      geom.pPhi = 3.014605;
      geom.pAlp1 = 0.035387;
      geom.pAlp2 = 0.035372;
      geom.pDy1 = 1.128053;
      geom.pDx1 = 1.201484;
      geom.pDx2 = 1.189943;
      geom.pDy2 = 1.262406;
      geom.pDx3 = 1.353618;
      geom.pDx4 = 1.340707;
      geom.centralX = -8.882151;
      geom.centralY = 104.213247;
      geom.centralZ = -48.242347;
      geom.pRotationAngleX = -1.993723;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 812 based Row/Col = 81/2
      geom_tower geom;
      geom.id = 812;
      geom.pDz = 6.751937;
      geom.pTheta = 0.056623;
      geom.pPhi = 2.964544;
      geom.pAlp1 = 0.025209;
      geom.pAlp2 = 0.025198;
      geom.pDy1 = 1.128053;
      geom.pDx1 = 1.197861;
      geom.pDx2 = 1.186419;
      geom.pDy2 = 1.262406;
      geom.pDx3 = 1.349537;
      geom.pDx4 = 1.336737;
      geom.centralX = -6.338074;
      geom.centralY = 104.213247;
      geom.centralZ = -48.242347;
      geom.pRotationAngleX = -1.993723;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 813 based Row/Col = 81/3
      geom_tower geom;
      geom.id = 813;
      geom.pDz = 6.751937;
      geom.pTheta = 0.034900;
      geom.pPhi = 2.851602;
      geom.pAlp1 = 0.015099;
      geom.pAlp2 = 0.015092;
      geom.pDy1 = 1.128053;
      geom.pDx1 = 1.195454;
      geom.pDx2 = 1.184077;
      geom.pDy2 = 1.262406;
      geom.pDx3 = 1.346826;
      geom.pDx4 = 1.334099;
      geom.centralX = -3.800322;
      geom.centralY = 104.213247;
      geom.centralZ = -48.242347;
      geom.pRotationAngleX = -1.993723;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 814 based Row/Col = 81/4
      geom_tower geom;
      geom.id = 814;
      geom.pDz = 6.751937;
      geom.pTheta = 0.014964;
      geom.pPhi = 2.411267;
      geom.pAlp1 = 0.005028;
      geom.pAlp2 = 0.005026;
      geom.pDy1 = 1.128053;
      geom.pDx1 = 1.194253;
      geom.pDx2 = 1.182908;
      geom.pDy2 = 1.262406;
      geom.pDx3 = 1.345473;
      geom.pDx4 = 1.332782;
      geom.centralX = -1.266354;
      geom.centralY = 104.213247;
      geom.centralZ = -48.242347;
      geom.pRotationAngleX = -1.993723;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 815 based Row/Col = 81/5
      geom_tower geom;
      geom.id = 815;
      geom.pDz = 6.751937;
      geom.pTheta = 0.014964;
      geom.pPhi = 0.730325;
      geom.pAlp1 = -0.005028;
      geom.pAlp2 = -0.005026;
      geom.pDy1 = 1.128053;
      geom.pDx1 = 1.194253;
      geom.pDx2 = 1.182908;
      geom.pDy2 = 1.262406;
      geom.pDx3 = 1.345473;
      geom.pDx4 = 1.332782;
      geom.centralX = 1.266354;
      geom.centralY = 104.213247;
      geom.centralZ = -48.242347;
      geom.pRotationAngleX = -1.993723;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 816 based Row/Col = 81/6
      geom_tower geom;
      geom.id = 816;
      geom.pDz = 6.751937;
      geom.pTheta = 0.034900;
      geom.pPhi = 0.289990;
      geom.pAlp1 = -0.015099;
      geom.pAlp2 = -0.015092;
      geom.pDy1 = 1.128053;
      geom.pDx1 = 1.195454;
      geom.pDx2 = 1.184077;
      geom.pDy2 = 1.262406;
      geom.pDx3 = 1.346826;
      geom.pDx4 = 1.334099;
      geom.centralX = 3.800322;
      geom.centralY = 104.213247;
      geom.centralZ = -48.242347;
      geom.pRotationAngleX = -1.993723;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 817 based Row/Col = 81/7
      geom_tower geom;
      geom.id = 817;
      geom.pDz = 6.751937;
      geom.pTheta = 0.056623;
      geom.pPhi = 0.177049;
      geom.pAlp1 = -0.025209;
      geom.pAlp2 = -0.025198;
      geom.pDy1 = 1.128053;
      geom.pDx1 = 1.197861;
      geom.pDx2 = 1.186419;
      geom.pDy2 = 1.262406;
      geom.pDx3 = 1.349537;
      geom.pDx4 = 1.336737;
      geom.centralX = 6.338074;
      geom.centralY = 104.213247;
      geom.centralZ = -48.242347;
      geom.pRotationAngleX = -1.993723;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 818 based Row/Col = 81/8
      geom_tower geom;
      geom.id = 818;
      geom.pDz = 6.751937;
      geom.pTheta = 0.078667;
      geom.pPhi = 0.126987;
      geom.pAlp1 = -0.035387;
      geom.pAlp2 = -0.035372;
      geom.pDy1 = 1.128053;
      geom.pDx1 = 1.201484;
      geom.pDx2 = 1.189943;
      geom.pDy2 = 1.262406;
      geom.pDx3 = 1.353618;
      geom.pDx4 = 1.340707;
      geom.centralX = 8.882151;
      geom.centralY = 104.213247;
      geom.centralZ = -48.242347;
      geom.pRotationAngleX = -1.993723;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 821 based Row/Col = 82/1
      geom_tower geom;
      geom.id = 821;
      geom.pDz = 6.751840;
      geom.pTheta = 0.079351;
      geom.pPhi = -3.016577;
      geom.pAlp1 = 0.035386;
      geom.pAlp2 = 0.035371;
      geom.pDy1 = 1.128053;
      geom.pDx1 = 1.213027;
      geom.pDx2 = 1.201484;
      geom.pDy2 = 1.262406;
      geom.pDx3 = 1.366531;
      geom.pDx4 = 1.353618;
      geom.centralX = -8.966758;
      geom.centralY = 105.196418;
      geom.centralZ = -46.057947;
      geom.pRotationAngleX = -1.993723;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 822 based Row/Col = 82/2
      geom_tower geom;
      geom.id = 822;
      geom.pDz = 6.751840;
      geom.pTheta = 0.057103;
      geom.pPhi = -2.967263;
      geom.pAlp1 = 0.025207;
      geom.pAlp2 = 0.025196;
      geom.pDy1 = 1.128053;
      geom.pDx1 = 1.209304;
      geom.pDx2 = 1.197861;
      geom.pDy2 = 1.262406;
      geom.pDx3 = 1.362338;
      geom.pDx4 = 1.349537;
      geom.centralX = -6.398332;
      geom.centralY = 105.196418;
      geom.centralZ = -46.057947;
      geom.pRotationAngleX = -1.993723;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 823 based Row/Col = 82/3
      geom_tower geom;
      geom.id = 823;
      geom.pDz = 6.751840;
      geom.pTheta = 0.035168;
      geom.pPhi = -2.855900;
      geom.pAlp1 = 0.015097;
      geom.pAlp2 = 0.015090;
      geom.pDy1 = 1.128053;
      geom.pDx1 = 1.206830;
      geom.pDx2 = 1.195454;
      geom.pDy2 = 1.262406;
      geom.pDx3 = 1.359552;
      geom.pDx4 = 1.346826;
      geom.centralX = -3.836407;
      geom.centralY = 105.196418;
      geom.centralZ = -46.057947;
      geom.pRotationAngleX = -1.993723;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 824 based Row/Col = 82/4
      geom_tower geom;
      geom.id = 824;
      geom.pDz = 6.751840;
      geom.pTheta = 0.014993;
      geom.pPhi = -2.419104;
      geom.pAlp1 = 0.005028;
      geom.pAlp2 = 0.005026;
      geom.pDy1 = 1.128053;
      geom.pDx1 = 1.205596;
      geom.pDx2 = 1.194253;
      geom.pDy2 = 1.262406;
      geom.pDx3 = 1.358161;
      geom.pDx4 = 1.345473;
      geom.centralX = -1.278371;
      geom.centralY = 105.196418;
      geom.centralZ = -46.057947;
      geom.pRotationAngleX = -1.993723;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 825 based Row/Col = 82/5
      geom_tower geom;
      geom.id = 825;
      geom.pDz = 6.751840;
      geom.pTheta = 0.014993;
      geom.pPhi = -0.722489;
      geom.pAlp1 = -0.005028;
      geom.pAlp2 = -0.005026;
      geom.pDy1 = 1.128053;
      geom.pDx1 = 1.205596;
      geom.pDx2 = 1.194253;
      geom.pDy2 = 1.262406;
      geom.pDx3 = 1.358161;
      geom.pDx4 = 1.345473;
      geom.centralX = 1.278371;
      geom.centralY = 105.196418;
      geom.centralZ = -46.057947;
      geom.pRotationAngleX = -1.993723;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 826 based Row/Col = 82/6
      geom_tower geom;
      geom.id = 826;
      geom.pDz = 6.751840;
      geom.pTheta = 0.035168;
      geom.pPhi = -0.285692;
      geom.pAlp1 = -0.015097;
      geom.pAlp2 = -0.015090;
      geom.pDy1 = 1.128053;
      geom.pDx1 = 1.206830;
      geom.pDx2 = 1.195454;
      geom.pDy2 = 1.262406;
      geom.pDx3 = 1.359552;
      geom.pDx4 = 1.346826;
      geom.centralX = 3.836407;
      geom.centralY = 105.196418;
      geom.centralZ = -46.057947;
      geom.pRotationAngleX = -1.993723;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 827 based Row/Col = 82/7
      geom_tower geom;
      geom.id = 827;
      geom.pDz = 6.751840;
      geom.pTheta = 0.057103;
      geom.pPhi = -0.174330;
      geom.pAlp1 = -0.025207;
      geom.pAlp2 = -0.025196;
      geom.pDy1 = 1.128053;
      geom.pDx1 = 1.209304;
      geom.pDx2 = 1.197861;
      geom.pDy2 = 1.262406;
      geom.pDx3 = 1.362338;
      geom.pDx4 = 1.349537;
      geom.centralX = 6.398332;
      geom.centralY = 105.196418;
      geom.centralZ = -46.057947;
      geom.pRotationAngleX = -1.993723;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 828 based Row/Col = 82/8
      geom_tower geom;
      geom.id = 828;
      geom.pDz = 6.751840;
      geom.pTheta = 0.079351;
      geom.pPhi = -0.125015;
      geom.pAlp1 = -0.035386;
      geom.pAlp2 = -0.035371;
      geom.pDy1 = 1.128053;
      geom.pDx1 = 1.213027;
      geom.pDx2 = 1.201484;
      geom.pDy2 = 1.262406;
      geom.pDx3 = 1.366531;
      geom.pDx4 = 1.353618;
      geom.centralX = 8.966758;
      geom.centralY = 105.196418;
      geom.centralZ = -46.057947;
      geom.pRotationAngleX = -1.993723;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1241 based Row/Col = 124/1
      geom_tower geom;
      geom.id = 1241;
      geom.pDz = 6.751745;
      geom.pTheta = 0.075465;
      geom.pPhi = -3.014668;
      geom.pAlp1 = -0.041687;
      geom.pAlp2 = -0.041680;
      geom.pDy1 = 1.130090;
      geom.pDx1 = 1.191410;
      geom.pDx2 = 1.205020;
      geom.pDy2 = 1.258587;
      geom.pDx3 = 1.335848;
      geom.pDx4 = 1.351003;
      geom.centralX = -8.880661;
      geom.centralY = 104.191848;
      geom.centralZ = 58.960290;
      geom.pRotationAngleX = -1.066053;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1242 based Row/Col = 124/2
      geom_tower geom;
      geom.id = 1242;
      geom.pDz = 6.751745;
      geom.pTheta = 0.054318;
      geom.pPhi = -2.964645;
      geom.pAlp1 = -0.029707;
      geom.pAlp2 = -0.029701;
      geom.pDy1 = 1.130090;
      geom.pDx1 = 1.188167;
      geom.pDx2 = 1.201671;
      geom.pDy2 = 1.258587;
      geom.pDx3 = 1.332213;
      geom.pDx4 = 1.347249;
      geom.centralX = -6.337515;
      geom.centralY = 104.191848;
      geom.centralZ = 58.960290;
      geom.pRotationAngleX = -1.066053;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1243 based Row/Col = 124/3
      geom_tower geom;
      geom.id = 1243;
      geom.pDz = 6.751745;
      geom.pTheta = 0.033478;
      geom.pPhi = -2.851777;
      geom.pAlp1 = -0.017796;
      geom.pAlp2 = -0.017792;
      geom.pDy1 = 1.130090;
      geom.pDx1 = 1.186012;
      geom.pDx2 = 1.199445;
      geom.pDy2 = 1.258587;
      geom.pDx3 = 1.329798;
      geom.pDx4 = 1.344754;
      geom.centralX = -3.800188;
      geom.centralY = 104.191848;
      geom.centralZ = 58.960290;
      geom.pRotationAngleX = -1.066053;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1244 based Row/Col = 124/4
      geom_tower geom;
      geom.id = 1244;
      geom.pDz = 6.751745;
      geom.pTheta = 0.014351;
      geom.pPhi = -2.411597;
      geom.pAlp1 = -0.005927;
      geom.pAlp2 = -0.005926;
      geom.pDy1 = 1.130090;
      geom.pDx1 = 1.184937;
      geom.pDx2 = 1.198334;
      geom.pDy2 = 1.258587;
      geom.pDx3 = 1.328592;
      geom.pDx4 = 1.343509;
      geom.centralX = -1.266343;
      geom.centralY = 104.191848;
      geom.centralZ = 58.960290;
      geom.pRotationAngleX = -1.066053;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1245 based Row/Col = 124/5
      geom_tower geom;
      geom.id = 1245;
      geom.pDz = 6.751745;
      geom.pTheta = 0.014351;
      geom.pPhi = -0.729996;
      geom.pAlp1 = 0.005927;
      geom.pAlp2 = 0.005926;
      geom.pDy1 = 1.130090;
      geom.pDx1 = 1.184937;
      geom.pDx2 = 1.198334;
      geom.pDy2 = 1.258587;
      geom.pDx3 = 1.328592;
      geom.pDx4 = 1.343509;
      geom.centralX = 1.266343;
      geom.centralY = 104.191848;
      geom.centralZ = 58.960290;
      geom.pRotationAngleX = -1.066053;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1246 based Row/Col = 124/6
      geom_tower geom;
      geom.id = 1246;
      geom.pDz = 6.751745;
      geom.pTheta = 0.033478;
      geom.pPhi = -0.289816;
      geom.pAlp1 = 0.017796;
      geom.pAlp2 = 0.017792;
      geom.pDy1 = 1.130090;
      geom.pDx1 = 1.186012;
      geom.pDx2 = 1.199445;
      geom.pDy2 = 1.258587;
      geom.pDx3 = 1.329798;
      geom.pDx4 = 1.344754;
      geom.centralX = 3.800188;
      geom.centralY = 104.191848;
      geom.centralZ = 58.960290;
      geom.pRotationAngleX = -1.066053;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1247 based Row/Col = 124/7
      geom_tower geom;
      geom.id = 1247;
      geom.pDz = 6.751745;
      geom.pTheta = 0.054318;
      geom.pPhi = -0.176948;
      geom.pAlp1 = 0.029707;
      geom.pAlp2 = 0.029701;
      geom.pDy1 = 1.130090;
      geom.pDx1 = 1.188167;
      geom.pDx2 = 1.201671;
      geom.pDy2 = 1.258587;
      geom.pDx3 = 1.332213;
      geom.pDx4 = 1.347249;
      geom.centralX = 6.337515;
      geom.centralY = 104.191848;
      geom.centralZ = 58.960290;
      geom.pRotationAngleX = -1.066053;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1248 based Row/Col = 124/8
      geom_tower geom;
      geom.id = 1248;
      geom.pDz = 6.751745;
      geom.pTheta = 0.075465;
      geom.pPhi = -0.126924;
      geom.pAlp1 = 0.041687;
      geom.pAlp2 = 0.041680;
      geom.pDy1 = 1.130090;
      geom.pDx1 = 1.191410;
      geom.pDx2 = 1.205020;
      geom.pDy2 = 1.258587;
      geom.pDx3 = 1.335848;
      geom.pDx4 = 1.351003;
      geom.centralX = 8.880661;
      geom.centralY = 104.191848;
      geom.centralZ = 58.960290;
      geom.pRotationAngleX = -1.066053;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1231 based Row/Col = 123/1
      geom_tower geom;
      geom.id = 1231;
      geom.pDz = 6.751777;
      geom.pTheta = 0.076232;
      geom.pPhi = 3.017407;
      geom.pAlp1 = -0.041686;
      geom.pAlp2 = -0.041678;
      geom.pDy1 = 1.130090;
      geom.pDx1 = 1.205020;
      geom.pDx2 = 1.218632;
      geom.pDy2 = 1.258587;
      geom.pDx3 = 1.351003;
      geom.pDx4 = 1.366162;
      geom.centralX = -8.980285;
      geom.centralY = 105.349389;
      geom.centralZ = 56.865107;
      geom.pRotationAngleX = -1.066053;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1232 based Row/Col = 123/2
      geom_tower geom;
      geom.id = 1232;
      geom.pDz = 6.751777;
      geom.pTheta = 0.054853;
      geom.pPhi = 2.968422;
      geom.pAlp1 = -0.029704;
      geom.pAlp2 = -0.029698;
      geom.pDy1 = 1.130090;
      geom.pDx1 = 1.201671;
      geom.pDx2 = 1.215174;
      geom.pDy2 = 1.258587;
      geom.pDx3 = 1.347249;
      geom.pDx4 = 1.362286;
      geom.centralX = -6.408486;
      geom.centralY = 105.349389;
      geom.centralZ = 56.865107;
      geom.pRotationAngleX = -1.066053;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1233 based Row/Col = 123/3
      geom_tower geom;
      geom.id = 1233;
      geom.pDz = 6.751777;
      geom.pTheta = 0.033771;
      geom.pPhi = 2.857750;
      geom.pAlp1 = -0.017794;
      geom.pAlp2 = -0.017790;
      geom.pDy1 = 1.130090;
      geom.pDx1 = 1.199445;
      geom.pDx2 = 1.212876;
      geom.pDy2 = 1.258587;
      geom.pDx3 = 1.344754;
      geom.pDx4 = 1.359710;
      geom.centralX = -3.842694;
      geom.centralY = 105.349389;
      geom.centralZ = 56.865107;
      geom.pRotationAngleX = -1.066053;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1234 based Row/Col = 123/4
      geom_tower geom;
      geom.id = 1234;
      geom.pDz = 6.751777;
      geom.pTheta = 0.014362;
      geom.pPhi = 2.422519;
      geom.pAlp1 = -0.005927;
      geom.pAlp2 = -0.005925;
      geom.pDy1 = 1.130090;
      geom.pDx1 = 1.198334;
      geom.pDx2 = 1.211729;
      geom.pDy2 = 1.258587;
      geom.pDx3 = 1.343509;
      geom.pDx4 = 1.358424;
      geom.centralX = -1.280499;
      geom.centralY = 105.349389;
      geom.centralZ = 56.865107;
      geom.pRotationAngleX = -1.066053;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1235 based Row/Col = 123/5
      geom_tower geom;
      geom.id = 1235;
      geom.pDz = 6.751777;
      geom.pTheta = 0.014362;
      geom.pPhi = 0.719074;
      geom.pAlp1 = 0.005927;
      geom.pAlp2 = 0.005925;
      geom.pDy1 = 1.130090;
      geom.pDx1 = 1.198334;
      geom.pDx2 = 1.211729;
      geom.pDy2 = 1.258587;
      geom.pDx3 = 1.343509;
      geom.pDx4 = 1.358424;
      geom.centralX = 1.280499;
      geom.centralY = 105.349389;
      geom.centralZ = 56.865107;
      geom.pRotationAngleX = -1.066053;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1236 based Row/Col = 123/6
      geom_tower geom;
      geom.id = 1236;
      geom.pDz = 6.751777;
      geom.pTheta = 0.033771;
      geom.pPhi = 0.283842;
      geom.pAlp1 = 0.017794;
      geom.pAlp2 = 0.017790;
      geom.pDy1 = 1.130090;
      geom.pDx1 = 1.199445;
      geom.pDx2 = 1.212876;
      geom.pDy2 = 1.258587;
      geom.pDx3 = 1.344754;
      geom.pDx4 = 1.359710;
      geom.centralX = 3.842694;
      geom.centralY = 105.349389;
      geom.centralZ = 56.865107;
      geom.pRotationAngleX = -1.066053;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1237 based Row/Col = 123/7
      geom_tower geom;
      geom.id = 1237;
      geom.pDz = 6.751777;
      geom.pTheta = 0.054853;
      geom.pPhi = 0.173170;
      geom.pAlp1 = 0.029704;
      geom.pAlp2 = 0.029698;
      geom.pDy1 = 1.130090;
      geom.pDx1 = 1.201671;
      geom.pDx2 = 1.215174;
      geom.pDy2 = 1.258587;
      geom.pDx3 = 1.347249;
      geom.pDx4 = 1.362286;
      geom.centralX = 6.408486;
      geom.centralY = 105.349389;
      geom.centralZ = 56.865107;
      geom.pRotationAngleX = -1.066053;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1238 based Row/Col = 123/8
      geom_tower geom;
      geom.id = 1238;
      geom.pDz = 6.751777;
      geom.pTheta = 0.076232;
      geom.pPhi = 0.124185;
      geom.pAlp1 = 0.041686;
      geom.pAlp2 = 0.041678;
      geom.pDy1 = 1.130090;
      geom.pDx1 = 1.205020;
      geom.pDx2 = 1.218632;
      geom.pDy2 = 1.258587;
      geom.pDx3 = 1.351003;
      geom.pDx4 = 1.366162;
      geom.centralX = 8.980285;
      geom.centralY = 105.349389;
      geom.centralZ = 56.865107;
      geom.pRotationAngleX = -1.066053;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 771 based Row/Col = 77/1
      geom_tower geom;
      geom.id = 771;
      geom.pDz = 6.751745;
      geom.pTheta = 0.075465;
      geom.pPhi = 3.014668;
      geom.pAlp1 = 0.041687;
      geom.pAlp2 = 0.041680;
      geom.pDy1 = 1.130090;
      geom.pDx1 = 1.205020;
      geom.pDx2 = 1.191410;
      geom.pDy2 = 1.258587;
      geom.pDx3 = 1.351003;
      geom.pDx4 = 1.335848;
      geom.centralX = -8.880661;
      geom.centralY = 104.191848;
      geom.centralZ = -58.960290;
      geom.pRotationAngleX = -2.075540;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 772 based Row/Col = 77/2
      geom_tower geom;
      geom.id = 772;
      geom.pDz = 6.751745;
      geom.pTheta = 0.054318;
      geom.pPhi = 2.964645;
      geom.pAlp1 = 0.029707;
      geom.pAlp2 = 0.029701;
      geom.pDy1 = 1.130090;
      geom.pDx1 = 1.201671;
      geom.pDx2 = 1.188167;
      geom.pDy2 = 1.258587;
      geom.pDx3 = 1.347249;
      geom.pDx4 = 1.332213;
      geom.centralX = -6.337515;
      geom.centralY = 104.191848;
      geom.centralZ = -58.960290;
      geom.pRotationAngleX = -2.075540;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 773 based Row/Col = 77/3
      geom_tower geom;
      geom.id = 773;
      geom.pDz = 6.751745;
      geom.pTheta = 0.033478;
      geom.pPhi = 2.851777;
      geom.pAlp1 = 0.017796;
      geom.pAlp2 = 0.017792;
      geom.pDy1 = 1.130090;
      geom.pDx1 = 1.199445;
      geom.pDx2 = 1.186012;
      geom.pDy2 = 1.258587;
      geom.pDx3 = 1.344754;
      geom.pDx4 = 1.329798;
      geom.centralX = -3.800188;
      geom.centralY = 104.191848;
      geom.centralZ = -58.960290;
      geom.pRotationAngleX = -2.075540;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 774 based Row/Col = 77/4
      geom_tower geom;
      geom.id = 774;
      geom.pDz = 6.751745;
      geom.pTheta = 0.014351;
      geom.pPhi = 2.411597;
      geom.pAlp1 = 0.005927;
      geom.pAlp2 = 0.005926;
      geom.pDy1 = 1.130090;
      geom.pDx1 = 1.198334;
      geom.pDx2 = 1.184937;
      geom.pDy2 = 1.258587;
      geom.pDx3 = 1.343509;
      geom.pDx4 = 1.328592;
      geom.centralX = -1.266343;
      geom.centralY = 104.191848;
      geom.centralZ = -58.960290;
      geom.pRotationAngleX = -2.075540;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 775 based Row/Col = 77/5
      geom_tower geom;
      geom.id = 775;
      geom.pDz = 6.751745;
      geom.pTheta = 0.014351;
      geom.pPhi = 0.729996;
      geom.pAlp1 = -0.005927;
      geom.pAlp2 = -0.005926;
      geom.pDy1 = 1.130090;
      geom.pDx1 = 1.198334;
      geom.pDx2 = 1.184937;
      geom.pDy2 = 1.258587;
      geom.pDx3 = 1.343509;
      geom.pDx4 = 1.328592;
      geom.centralX = 1.266343;
      geom.centralY = 104.191848;
      geom.centralZ = -58.960290;
      geom.pRotationAngleX = -2.075540;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 776 based Row/Col = 77/6
      geom_tower geom;
      geom.id = 776;
      geom.pDz = 6.751745;
      geom.pTheta = 0.033478;
      geom.pPhi = 0.289816;
      geom.pAlp1 = -0.017796;
      geom.pAlp2 = -0.017792;
      geom.pDy1 = 1.130090;
      geom.pDx1 = 1.199445;
      geom.pDx2 = 1.186012;
      geom.pDy2 = 1.258587;
      geom.pDx3 = 1.344754;
      geom.pDx4 = 1.329798;
      geom.centralX = 3.800188;
      geom.centralY = 104.191848;
      geom.centralZ = -58.960290;
      geom.pRotationAngleX = -2.075540;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 777 based Row/Col = 77/7
      geom_tower geom;
      geom.id = 777;
      geom.pDz = 6.751745;
      geom.pTheta = 0.054318;
      geom.pPhi = 0.176948;
      geom.pAlp1 = -0.029707;
      geom.pAlp2 = -0.029701;
      geom.pDy1 = 1.130090;
      geom.pDx1 = 1.201671;
      geom.pDx2 = 1.188167;
      geom.pDy2 = 1.258587;
      geom.pDx3 = 1.347249;
      geom.pDx4 = 1.332213;
      geom.centralX = 6.337515;
      geom.centralY = 104.191848;
      geom.centralZ = -58.960290;
      geom.pRotationAngleX = -2.075540;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 778 based Row/Col = 77/8
      geom_tower geom;
      geom.id = 778;
      geom.pDz = 6.751745;
      geom.pTheta = 0.075465;
      geom.pPhi = 0.126924;
      geom.pAlp1 = -0.041687;
      geom.pAlp2 = -0.041680;
      geom.pDy1 = 1.130090;
      geom.pDx1 = 1.205020;
      geom.pDx2 = 1.191410;
      geom.pDy2 = 1.258587;
      geom.pDx3 = 1.351003;
      geom.pDx4 = 1.335848;
      geom.centralX = 8.880661;
      geom.centralY = 104.191848;
      geom.centralZ = -58.960290;
      geom.pRotationAngleX = -2.075540;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 781 based Row/Col = 78/1
      geom_tower geom;
      geom.id = 781;
      geom.pDz = 6.751777;
      geom.pTheta = 0.076232;
      geom.pPhi = -3.017407;
      geom.pAlp1 = 0.041686;
      geom.pAlp2 = 0.041678;
      geom.pDy1 = 1.130090;
      geom.pDx1 = 1.218632;
      geom.pDx2 = 1.205020;
      geom.pDy2 = 1.258587;
      geom.pDx3 = 1.366162;
      geom.pDx4 = 1.351003;
      geom.centralX = -8.980285;
      geom.centralY = 105.349389;
      geom.centralZ = -56.865107;
      geom.pRotationAngleX = -2.075540;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 782 based Row/Col = 78/2
      geom_tower geom;
      geom.id = 782;
      geom.pDz = 6.751777;
      geom.pTheta = 0.054853;
      geom.pPhi = -2.968422;
      geom.pAlp1 = 0.029704;
      geom.pAlp2 = 0.029698;
      geom.pDy1 = 1.130090;
      geom.pDx1 = 1.215174;
      geom.pDx2 = 1.201671;
      geom.pDy2 = 1.258587;
      geom.pDx3 = 1.362286;
      geom.pDx4 = 1.347249;
      geom.centralX = -6.408486;
      geom.centralY = 105.349389;
      geom.centralZ = -56.865107;
      geom.pRotationAngleX = -2.075540;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 783 based Row/Col = 78/3
      geom_tower geom;
      geom.id = 783;
      geom.pDz = 6.751777;
      geom.pTheta = 0.033771;
      geom.pPhi = -2.857750;
      geom.pAlp1 = 0.017794;
      geom.pAlp2 = 0.017790;
      geom.pDy1 = 1.130090;
      geom.pDx1 = 1.212876;
      geom.pDx2 = 1.199445;
      geom.pDy2 = 1.258587;
      geom.pDx3 = 1.359710;
      geom.pDx4 = 1.344754;
      geom.centralX = -3.842694;
      geom.centralY = 105.349389;
      geom.centralZ = -56.865107;
      geom.pRotationAngleX = -2.075540;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 784 based Row/Col = 78/4
      geom_tower geom;
      geom.id = 784;
      geom.pDz = 6.751777;
      geom.pTheta = 0.014362;
      geom.pPhi = -2.422519;
      geom.pAlp1 = 0.005927;
      geom.pAlp2 = 0.005925;
      geom.pDy1 = 1.130090;
      geom.pDx1 = 1.211729;
      geom.pDx2 = 1.198334;
      geom.pDy2 = 1.258587;
      geom.pDx3 = 1.358424;
      geom.pDx4 = 1.343509;
      geom.centralX = -1.280499;
      geom.centralY = 105.349389;
      geom.centralZ = -56.865107;
      geom.pRotationAngleX = -2.075540;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 785 based Row/Col = 78/5
      geom_tower geom;
      geom.id = 785;
      geom.pDz = 6.751777;
      geom.pTheta = 0.014362;
      geom.pPhi = -0.719074;
      geom.pAlp1 = -0.005927;
      geom.pAlp2 = -0.005925;
      geom.pDy1 = 1.130090;
      geom.pDx1 = 1.211729;
      geom.pDx2 = 1.198334;
      geom.pDy2 = 1.258587;
      geom.pDx3 = 1.358424;
      geom.pDx4 = 1.343509;
      geom.centralX = 1.280499;
      geom.centralY = 105.349389;
      geom.centralZ = -56.865107;
      geom.pRotationAngleX = -2.075540;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 786 based Row/Col = 78/6
      geom_tower geom;
      geom.id = 786;
      geom.pDz = 6.751777;
      geom.pTheta = 0.033771;
      geom.pPhi = -0.283842;
      geom.pAlp1 = -0.017794;
      geom.pAlp2 = -0.017790;
      geom.pDy1 = 1.130090;
      geom.pDx1 = 1.212876;
      geom.pDx2 = 1.199445;
      geom.pDy2 = 1.258587;
      geom.pDx3 = 1.359710;
      geom.pDx4 = 1.344754;
      geom.centralX = 3.842694;
      geom.centralY = 105.349389;
      geom.centralZ = -56.865107;
      geom.pRotationAngleX = -2.075540;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 787 based Row/Col = 78/7
      geom_tower geom;
      geom.id = 787;
      geom.pDz = 6.751777;
      geom.pTheta = 0.054853;
      geom.pPhi = -0.173170;
      geom.pAlp1 = -0.029704;
      geom.pAlp2 = -0.029698;
      geom.pDy1 = 1.130090;
      geom.pDx1 = 1.215174;
      geom.pDx2 = 1.201671;
      geom.pDy2 = 1.258587;
      geom.pDx3 = 1.362286;
      geom.pDx4 = 1.347249;
      geom.centralX = 6.408486;
      geom.centralY = 105.349389;
      geom.centralZ = -56.865107;
      geom.pRotationAngleX = -2.075540;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 788 based Row/Col = 78/8
      geom_tower geom;
      geom.id = 788;
      geom.pDz = 6.751777;
      geom.pTheta = 0.076232;
      geom.pPhi = -0.124185;
      geom.pAlp1 = -0.041686;
      geom.pAlp2 = -0.041678;
      geom.pDy1 = 1.130090;
      geom.pDx1 = 1.218632;
      geom.pDx2 = 1.205020;
      geom.pDy2 = 1.258587;
      geom.pDx3 = 1.366162;
      geom.pDx4 = 1.351003;
      geom.centralX = 8.980285;
      geom.centralY = 105.349389;
      geom.centralZ = -56.865107;
      geom.pRotationAngleX = -2.075540;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1281 based Row/Col = 128/1
      geom_tower geom;
      geom.id = 1281;
      geom.pDz = 6.751723;
      geom.pTheta = 0.071943;
      geom.pPhi = -3.014883;
      geom.pAlp1 = -0.047448;
      geom.pAlp2 = -0.047437;
      geom.pDy1 = 1.132133;
      geom.pDx1 = 1.193895;
      geom.pDx2 = 1.209401;
      geom.pDy2 = 1.253515;
      geom.pDx3 = 1.331433;
      geom.pDx4 = 1.348597;
      geom.centralX = -8.883705;
      geom.centralY = 104.222868;
      geom.centralZ = 70.180073;
      geom.pRotationAngleX = -0.987936;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1282 based Row/Col = 128/2
      geom_tower geom;
      geom.id = 1282;
      geom.pDz = 6.751723;
      geom.pTheta = 0.051781;
      geom.pPhi = -2.964956;
      geom.pAlp1 = -0.033823;
      geom.pAlp2 = -0.033815;
      geom.pDy1 = 1.132133;
      geom.pDx1 = 1.190946;
      geom.pDx2 = 1.206342;
      geom.pDy2 = 1.253515;
      geom.pDx3 = 1.328144;
      geom.pDx4 = 1.345187;
      geom.centralX = -6.340218;
      geom.centralY = 104.222868;
      geom.centralZ = 70.180073;
      geom.pRotationAngleX = -0.987936;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1283 based Row/Col = 128/3
      geom_tower geom;
      geom.id = 1283;
      geom.pDz = 6.751723;
      geom.pTheta = 0.031912;
      geom.pPhi = -2.852284;
      geom.pAlp1 = -0.020267;
      geom.pAlp2 = -0.020261;
      geom.pDy1 = 1.132133;
      geom.pDx1 = 1.188986;
      geom.pDx2 = 1.204308;
      geom.pDy2 = 1.253515;
      geom.pDx3 = 1.325958;
      geom.pDx4 = 1.342920;
      geom.centralX = -3.802020;
      geom.centralY = 104.222868;
      geom.centralZ = 70.180073;
      geom.pRotationAngleX = -0.987936;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1284 based Row/Col = 128/4
      geom_tower geom;
      geom.id = 1284;
      geom.pDz = 6.751723;
      geom.pTheta = 0.013670;
      geom.pPhi = -2.412531;
      geom.pAlp1 = -0.006751;
      geom.pAlp2 = -0.006749;
      geom.pDy1 = 1.132133;
      geom.pDx1 = 1.188007;
      geom.pDx2 = 1.203293;
      geom.pDy2 = 1.253515;
      geom.pDx3 = 1.324867;
      geom.pDx4 = 1.341788;
      geom.centralX = -1.266989;
      geom.centralY = 104.222868;
      geom.centralZ = 70.180073;
      geom.pRotationAngleX = -0.987936;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1285 based Row/Col = 128/5
      geom_tower geom;
      geom.id = 1285;
      geom.pDz = 6.751723;
      geom.pTheta = 0.013670;
      geom.pPhi = -0.729061;
      geom.pAlp1 = 0.006751;
      geom.pAlp2 = 0.006749;
      geom.pDy1 = 1.132133;
      geom.pDx1 = 1.188007;
      geom.pDx2 = 1.203293;
      geom.pDy2 = 1.253515;
      geom.pDx3 = 1.324867;
      geom.pDx4 = 1.341788;
      geom.centralX = 1.266989;
      geom.centralY = 104.222868;
      geom.centralZ = 70.180073;
      geom.pRotationAngleX = -0.987936;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1286 based Row/Col = 128/6
      geom_tower geom;
      geom.id = 1286;
      geom.pDz = 6.751723;
      geom.pTheta = 0.031912;
      geom.pPhi = -0.289309;
      geom.pAlp1 = 0.020267;
      geom.pAlp2 = 0.020261;
      geom.pDy1 = 1.132133;
      geom.pDx1 = 1.188986;
      geom.pDx2 = 1.204308;
      geom.pDy2 = 1.253515;
      geom.pDx3 = 1.325958;
      geom.pDx4 = 1.342920;
      geom.centralX = 3.802020;
      geom.centralY = 104.222868;
      geom.centralZ = 70.180073;
      geom.pRotationAngleX = -0.987936;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1287 based Row/Col = 128/7
      geom_tower geom;
      geom.id = 1287;
      geom.pDz = 6.751723;
      geom.pTheta = 0.051781;
      geom.pPhi = -0.176637;
      geom.pAlp1 = 0.033823;
      geom.pAlp2 = 0.033815;
      geom.pDy1 = 1.132133;
      geom.pDx1 = 1.190946;
      geom.pDx2 = 1.206342;
      geom.pDy2 = 1.253515;
      geom.pDx3 = 1.328144;
      geom.pDx4 = 1.345187;
      geom.centralX = 6.340218;
      geom.centralY = 104.222868;
      geom.centralZ = 70.180073;
      geom.pRotationAngleX = -0.987936;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1288 based Row/Col = 128/8
      geom_tower geom;
      geom.id = 1288;
      geom.pDz = 6.751723;
      geom.pTheta = 0.071943;
      geom.pPhi = -0.126709;
      geom.pAlp1 = 0.047448;
      geom.pAlp2 = 0.047437;
      geom.pDy1 = 1.132133;
      geom.pDx1 = 1.193895;
      geom.pDx2 = 1.209401;
      geom.pDy2 = 1.253515;
      geom.pDx3 = 1.331433;
      geom.pDx4 = 1.348597;
      geom.centralX = 8.883705;
      geom.centralY = 104.222868;
      geom.centralZ = 70.180073;
      geom.pRotationAngleX = -0.987936;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1271 based Row/Col = 127/1
      geom_tower geom;
      geom.id = 1271;
      geom.pDz = 6.751714;
      geom.pTheta = 0.072754;
      geom.pPhi = 3.019577;
      geom.pAlp1 = -0.047446;
      geom.pAlp2 = -0.047435;
      geom.pDy1 = 1.132133;
      geom.pDx1 = 1.209401;
      geom.pDx2 = 1.224910;
      geom.pDy2 = 1.253515;
      geom.pDx3 = 1.348597;
      geom.pDx4 = 1.365766;
      geom.centralX = -8.996967;
      geom.centralY = 105.538715;
      geom.centralZ = 68.184141;
      geom.pRotationAngleX = -0.987936;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1272 based Row/Col = 127/2
      geom_tower geom;
      geom.id = 1272;
      geom.pDz = 6.751714;
      geom.pTheta = 0.052337;
      geom.pPhi = 2.971433;
      geom.pAlp1 = -0.033820;
      geom.pAlp2 = -0.033812;
      geom.pDy1 = 1.132133;
      geom.pDx1 = 1.206342;
      geom.pDx2 = 1.221738;
      geom.pDy2 = 1.253515;
      geom.pDx3 = 1.345187;
      geom.pDx4 = 1.362230;
      geom.centralX = -6.420925;
      geom.centralY = 105.538715;
      geom.centralZ = 68.184141;
      geom.pRotationAngleX = -0.987936;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1273 based Row/Col = 127/3
      geom_tower geom;
      geom.id = 1273;
      geom.pDz = 6.751714;
      geom.pTheta = 0.032194;
      geom.pPhi = 2.862538;
      geom.pAlp1 = -0.020264;
      geom.pAlp2 = -0.020259;
      geom.pDy1 = 1.132133;
      geom.pDx1 = 1.204308;
      geom.pDx2 = 1.219629;
      geom.pDy2 = 1.253515;
      geom.pDx3 = 1.342920;
      geom.pDx4 = 1.359880;
      geom.centralX = -3.850366;
      geom.centralY = 105.538715;
      geom.centralZ = 68.184141;
      geom.pRotationAngleX = -0.987936;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1274 based Row/Col = 127/4
      geom_tower geom;
      geom.id = 1274;
      geom.pDz = 6.751714;
      geom.pTheta = 0.013605;
      geom.pPhi = 2.431415;
      geom.pAlp1 = -0.006750;
      geom.pAlp2 = -0.006748;
      geom.pDy1 = 1.132133;
      geom.pDx1 = 1.203293;
      geom.pDx2 = 1.218577;
      geom.pDy2 = 1.253515;
      geom.pDx3 = 1.341788;
      geom.pDx4 = 1.358706;
      geom.centralX = -1.283091;
      geom.centralY = 105.538715;
      geom.centralZ = 68.184141;
      geom.pRotationAngleX = -0.987936;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1275 based Row/Col = 127/5
      geom_tower geom;
      geom.id = 1275;
      geom.pDz = 6.751714;
      geom.pTheta = 0.013605;
      geom.pPhi = 0.710178;
      geom.pAlp1 = 0.006750;
      geom.pAlp2 = 0.006748;
      geom.pDy1 = 1.132133;
      geom.pDx1 = 1.203293;
      geom.pDx2 = 1.218577;
      geom.pDy2 = 1.253515;
      geom.pDx3 = 1.341788;
      geom.pDx4 = 1.358706;
      geom.centralX = 1.283091;
      geom.centralY = 105.538715;
      geom.centralZ = 68.184141;
      geom.pRotationAngleX = -0.987936;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1276 based Row/Col = 127/6
      geom_tower geom;
      geom.id = 1276;
      geom.pDz = 6.751714;
      geom.pTheta = 0.032194;
      geom.pPhi = 0.279055;
      geom.pAlp1 = 0.020264;
      geom.pAlp2 = 0.020259;
      geom.pDy1 = 1.132133;
      geom.pDx1 = 1.204308;
      geom.pDx2 = 1.219629;
      geom.pDy2 = 1.253515;
      geom.pDx3 = 1.342920;
      geom.pDx4 = 1.359880;
      geom.centralX = 3.850366;
      geom.centralY = 105.538715;
      geom.centralZ = 68.184141;
      geom.pRotationAngleX = -0.987936;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1277 based Row/Col = 127/7
      geom_tower geom;
      geom.id = 1277;
      geom.pDz = 6.751714;
      geom.pTheta = 0.052337;
      geom.pPhi = 0.170160;
      geom.pAlp1 = 0.033820;
      geom.pAlp2 = 0.033812;
      geom.pDy1 = 1.132133;
      geom.pDx1 = 1.206342;
      geom.pDx2 = 1.221738;
      geom.pDy2 = 1.253515;
      geom.pDx3 = 1.345187;
      geom.pDx4 = 1.362230;
      geom.centralX = 6.420925;
      geom.centralY = 105.538715;
      geom.centralZ = 68.184141;
      geom.pRotationAngleX = -0.987936;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1278 based Row/Col = 127/8
      geom_tower geom;
      geom.id = 1278;
      geom.pDz = 6.751714;
      geom.pTheta = 0.072754;
      geom.pPhi = 0.122016;
      geom.pAlp1 = 0.047446;
      geom.pAlp2 = 0.047435;
      geom.pDy1 = 1.132133;
      geom.pDx1 = 1.209401;
      geom.pDx2 = 1.224910;
      geom.pDy2 = 1.253515;
      geom.pDx3 = 1.348597;
      geom.pDx4 = 1.365766;
      geom.centralX = 8.996967;
      geom.centralY = 105.538715;
      geom.centralZ = 68.184141;
      geom.pRotationAngleX = -0.987936;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 731 based Row/Col = 73/1
      geom_tower geom;
      geom.id = 731;
      geom.pDz = 6.751723;
      geom.pTheta = 0.071943;
      geom.pPhi = 3.014883;
      geom.pAlp1 = 0.047448;
      geom.pAlp2 = 0.047437;
      geom.pDy1 = 1.132133;
      geom.pDx1 = 1.209401;
      geom.pDx2 = 1.193895;
      geom.pDy2 = 1.253515;
      geom.pDx3 = 1.348597;
      geom.pDx4 = 1.331433;
      geom.centralX = -8.883705;
      geom.centralY = 104.222868;
      geom.centralZ = -70.180073;
      geom.pRotationAngleX = -2.153657;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 732 based Row/Col = 73/2
      geom_tower geom;
      geom.id = 732;
      geom.pDz = 6.751723;
      geom.pTheta = 0.051781;
      geom.pPhi = 2.964956;
      geom.pAlp1 = 0.033823;
      geom.pAlp2 = 0.033815;
      geom.pDy1 = 1.132133;
      geom.pDx1 = 1.206342;
      geom.pDx2 = 1.190946;
      geom.pDy2 = 1.253515;
      geom.pDx3 = 1.345187;
      geom.pDx4 = 1.328144;
      geom.centralX = -6.340218;
      geom.centralY = 104.222868;
      geom.centralZ = -70.180073;
      geom.pRotationAngleX = -2.153657;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 733 based Row/Col = 73/3
      geom_tower geom;
      geom.id = 733;
      geom.pDz = 6.751723;
      geom.pTheta = 0.031912;
      geom.pPhi = 2.852284;
      geom.pAlp1 = 0.020267;
      geom.pAlp2 = 0.020261;
      geom.pDy1 = 1.132133;
      geom.pDx1 = 1.204308;
      geom.pDx2 = 1.188986;
      geom.pDy2 = 1.253515;
      geom.pDx3 = 1.342920;
      geom.pDx4 = 1.325958;
      geom.centralX = -3.802020;
      geom.centralY = 104.222868;
      geom.centralZ = -70.180073;
      geom.pRotationAngleX = -2.153657;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 734 based Row/Col = 73/4
      geom_tower geom;
      geom.id = 734;
      geom.pDz = 6.751723;
      geom.pTheta = 0.013670;
      geom.pPhi = 2.412531;
      geom.pAlp1 = 0.006751;
      geom.pAlp2 = 0.006749;
      geom.pDy1 = 1.132133;
      geom.pDx1 = 1.203293;
      geom.pDx2 = 1.188007;
      geom.pDy2 = 1.253515;
      geom.pDx3 = 1.341788;
      geom.pDx4 = 1.324867;
      geom.centralX = -1.266989;
      geom.centralY = 104.222868;
      geom.centralZ = -70.180073;
      geom.pRotationAngleX = -2.153657;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 735 based Row/Col = 73/5
      geom_tower geom;
      geom.id = 735;
      geom.pDz = 6.751723;
      geom.pTheta = 0.013670;
      geom.pPhi = 0.729061;
      geom.pAlp1 = -0.006751;
      geom.pAlp2 = -0.006749;
      geom.pDy1 = 1.132133;
      geom.pDx1 = 1.203293;
      geom.pDx2 = 1.188007;
      geom.pDy2 = 1.253515;
      geom.pDx3 = 1.341788;
      geom.pDx4 = 1.324867;
      geom.centralX = 1.266989;
      geom.centralY = 104.222868;
      geom.centralZ = -70.180073;
      geom.pRotationAngleX = -2.153657;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 736 based Row/Col = 73/6
      geom_tower geom;
      geom.id = 736;
      geom.pDz = 6.751723;
      geom.pTheta = 0.031912;
      geom.pPhi = 0.289309;
      geom.pAlp1 = -0.020267;
      geom.pAlp2 = -0.020261;
      geom.pDy1 = 1.132133;
      geom.pDx1 = 1.204308;
      geom.pDx2 = 1.188986;
      geom.pDy2 = 1.253515;
      geom.pDx3 = 1.342920;
      geom.pDx4 = 1.325958;
      geom.centralX = 3.802020;
      geom.centralY = 104.222868;
      geom.centralZ = -70.180073;
      geom.pRotationAngleX = -2.153657;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 737 based Row/Col = 73/7
      geom_tower geom;
      geom.id = 737;
      geom.pDz = 6.751723;
      geom.pTheta = 0.051781;
      geom.pPhi = 0.176637;
      geom.pAlp1 = -0.033823;
      geom.pAlp2 = -0.033815;
      geom.pDy1 = 1.132133;
      geom.pDx1 = 1.206342;
      geom.pDx2 = 1.190946;
      geom.pDy2 = 1.253515;
      geom.pDx3 = 1.345187;
      geom.pDx4 = 1.328144;
      geom.centralX = 6.340218;
      geom.centralY = 104.222868;
      geom.centralZ = -70.180073;
      geom.pRotationAngleX = -2.153657;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 738 based Row/Col = 73/8
      geom_tower geom;
      geom.id = 738;
      geom.pDz = 6.751723;
      geom.pTheta = 0.071943;
      geom.pPhi = 0.126709;
      geom.pAlp1 = -0.047448;
      geom.pAlp2 = -0.047437;
      geom.pDy1 = 1.132133;
      geom.pDx1 = 1.209401;
      geom.pDx2 = 1.193895;
      geom.pDy2 = 1.253515;
      geom.pDx3 = 1.348597;
      geom.pDx4 = 1.331433;
      geom.centralX = 8.883705;
      geom.centralY = 104.222868;
      geom.centralZ = -70.180073;
      geom.pRotationAngleX = -2.153657;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 741 based Row/Col = 74/1
      geom_tower geom;
      geom.id = 741;
      geom.pDz = 6.751714;
      geom.pTheta = 0.072754;
      geom.pPhi = -3.019577;
      geom.pAlp1 = 0.047446;
      geom.pAlp2 = 0.047435;
      geom.pDy1 = 1.132133;
      geom.pDx1 = 1.224910;
      geom.pDx2 = 1.209401;
      geom.pDy2 = 1.253515;
      geom.pDx3 = 1.365766;
      geom.pDx4 = 1.348597;
      geom.centralX = -8.996967;
      geom.centralY = 105.538715;
      geom.centralZ = -68.184141;
      geom.pRotationAngleX = -2.153657;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 742 based Row/Col = 74/2
      geom_tower geom;
      geom.id = 742;
      geom.pDz = 6.751714;
      geom.pTheta = 0.052337;
      geom.pPhi = -2.971433;
      geom.pAlp1 = 0.033820;
      geom.pAlp2 = 0.033812;
      geom.pDy1 = 1.132133;
      geom.pDx1 = 1.221738;
      geom.pDx2 = 1.206342;
      geom.pDy2 = 1.253515;
      geom.pDx3 = 1.362230;
      geom.pDx4 = 1.345187;
      geom.centralX = -6.420925;
      geom.centralY = 105.538715;
      geom.centralZ = -68.184141;
      geom.pRotationAngleX = -2.153657;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 743 based Row/Col = 74/3
      geom_tower geom;
      geom.id = 743;
      geom.pDz = 6.751714;
      geom.pTheta = 0.032194;
      geom.pPhi = -2.862538;
      geom.pAlp1 = 0.020264;
      geom.pAlp2 = 0.020259;
      geom.pDy1 = 1.132133;
      geom.pDx1 = 1.219629;
      geom.pDx2 = 1.204308;
      geom.pDy2 = 1.253515;
      geom.pDx3 = 1.359880;
      geom.pDx4 = 1.342920;
      geom.centralX = -3.850366;
      geom.centralY = 105.538715;
      geom.centralZ = -68.184141;
      geom.pRotationAngleX = -2.153657;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 744 based Row/Col = 74/4
      geom_tower geom;
      geom.id = 744;
      geom.pDz = 6.751714;
      geom.pTheta = 0.013605;
      geom.pPhi = -2.431415;
      geom.pAlp1 = 0.006750;
      geom.pAlp2 = 0.006748;
      geom.pDy1 = 1.132133;
      geom.pDx1 = 1.218577;
      geom.pDx2 = 1.203293;
      geom.pDy2 = 1.253515;
      geom.pDx3 = 1.358706;
      geom.pDx4 = 1.341788;
      geom.centralX = -1.283091;
      geom.centralY = 105.538715;
      geom.centralZ = -68.184141;
      geom.pRotationAngleX = -2.153657;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 745 based Row/Col = 74/5
      geom_tower geom;
      geom.id = 745;
      geom.pDz = 6.751714;
      geom.pTheta = 0.013605;
      geom.pPhi = -0.710178;
      geom.pAlp1 = -0.006750;
      geom.pAlp2 = -0.006748;
      geom.pDy1 = 1.132133;
      geom.pDx1 = 1.218577;
      geom.pDx2 = 1.203293;
      geom.pDy2 = 1.253515;
      geom.pDx3 = 1.358706;
      geom.pDx4 = 1.341788;
      geom.centralX = 1.283091;
      geom.centralY = 105.538715;
      geom.centralZ = -68.184141;
      geom.pRotationAngleX = -2.153657;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 746 based Row/Col = 74/6
      geom_tower geom;
      geom.id = 746;
      geom.pDz = 6.751714;
      geom.pTheta = 0.032194;
      geom.pPhi = -0.279055;
      geom.pAlp1 = -0.020264;
      geom.pAlp2 = -0.020259;
      geom.pDy1 = 1.132133;
      geom.pDx1 = 1.219629;
      geom.pDx2 = 1.204308;
      geom.pDy2 = 1.253515;
      geom.pDx3 = 1.359880;
      geom.pDx4 = 1.342920;
      geom.centralX = 3.850366;
      geom.centralY = 105.538715;
      geom.centralZ = -68.184141;
      geom.pRotationAngleX = -2.153657;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 747 based Row/Col = 74/7
      geom_tower geom;
      geom.id = 747;
      geom.pDz = 6.751714;
      geom.pTheta = 0.052337;
      geom.pPhi = -0.170160;
      geom.pAlp1 = -0.033820;
      geom.pAlp2 = -0.033812;
      geom.pDy1 = 1.132133;
      geom.pDx1 = 1.221738;
      geom.pDx2 = 1.206342;
      geom.pDy2 = 1.253515;
      geom.pDx3 = 1.362230;
      geom.pDx4 = 1.345187;
      geom.centralX = 6.420925;
      geom.centralY = 105.538715;
      geom.centralZ = -68.184141;
      geom.pRotationAngleX = -2.153657;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 748 based Row/Col = 74/8
      geom_tower geom;
      geom.id = 748;
      geom.pDz = 6.751714;
      geom.pTheta = 0.072754;
      geom.pPhi = -0.122016;
      geom.pAlp1 = -0.047446;
      geom.pAlp2 = -0.047435;
      geom.pDy1 = 1.132133;
      geom.pDx1 = 1.224910;
      geom.pDx2 = 1.209401;
      geom.pDy2 = 1.253515;
      geom.pDx3 = 1.365766;
      geom.pDx4 = 1.348597;
      geom.centralX = 8.996967;
      geom.centralY = 105.538715;
      geom.centralZ = -68.184141;
      geom.pRotationAngleX = -2.153657;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1321 based Row/Col = 132/1
      geom_tower geom;
      geom.id = 1321;
      geom.pDz = 6.923155;
      geom.pTheta = 0.068202;
      geom.pPhi = -3.014437;
      geom.pAlp1 = -0.052649;
      geom.pAlp2 = -0.052634;
      geom.pDy1 = 1.129974;
      geom.pDx1 = 1.193922;
      geom.pDx2 = 1.211081;
      geom.pDy2 = 1.248374;
      geom.pDx3 = 1.327452;
      geom.pDx4 = 1.346405;
      geom.centralX = -8.878911;
      geom.centralY = 104.162860;
      geom.centralZ = 81.887600;
      geom.pRotationAngleX = -0.913765;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1322 based Row/Col = 132/2
      geom_tower geom;
      geom.id = 1322;
      geom.pDz = 6.923155;
      geom.pTheta = 0.049092;
      geom.pPhi = -2.964354;
      geom.pAlp1 = -0.037543;
      geom.pAlp2 = -0.037532;
      geom.pDy1 = 1.129974;
      geom.pDx1 = 1.191275;
      geom.pDx2 = 1.208325;
      geom.pDy2 = 1.248374;
      geom.pDx3 = 1.324510;
      geom.pDx4 = 1.343341;
      geom.centralX = -6.337333;
      geom.centralY = 104.162860;
      geom.centralZ = 81.887600;
      geom.pRotationAngleX = -0.913765;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1323 based Row/Col = 132/3
      geom_tower geom;
      geom.id = 1323;
      geom.pDz = 6.923155;
      geom.pTheta = 0.030259;
      geom.pPhi = -2.851349;
      geom.pAlp1 = -0.022500;
      geom.pAlp2 = -0.022493;
      geom.pDy1 = 1.129974;
      geom.pDx1 = 1.189515;
      geom.pDx2 = 1.206492;
      geom.pDy2 = 1.248374;
      geom.pDx3 = 1.322554;
      geom.pDx4 = 1.341304;
      geom.centralX = -3.800504;
      geom.centralY = 104.162860;
      geom.centralZ = 81.887600;
      geom.pRotationAngleX = -0.913765;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1324 based Row/Col = 132/4
      geom_tower geom;
      geom.id = 1324;
      geom.pDz = 6.923155;
      geom.pTheta = 0.012978;
      geom.pPhi = -2.410848;
      geom.pAlp1 = -0.007496;
      geom.pAlp2 = -0.007494;
      geom.pDy1 = 1.129974;
      geom.pDx1 = 1.188636;
      geom.pDx2 = 1.205577;
      geom.pDy2 = 1.248374;
      geom.pDx3 = 1.321577;
      geom.pDx4 = 1.340287;
      geom.centralX = -1.266519;
      geom.centralY = 104.162860;
      geom.centralZ = 81.887600;
      geom.pRotationAngleX = -0.913765;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1325 based Row/Col = 132/5
      geom_tower geom;
      geom.id = 1325;
      geom.pDz = 6.923155;
      geom.pTheta = 0.012978;
      geom.pPhi = -0.730745;
      geom.pAlp1 = 0.007496;
      geom.pAlp2 = 0.007494;
      geom.pDy1 = 1.129974;
      geom.pDx1 = 1.188636;
      geom.pDx2 = 1.205577;
      geom.pDy2 = 1.248374;
      geom.pDx3 = 1.321577;
      geom.pDx4 = 1.340287;
      geom.centralX = 1.266519;
      geom.centralY = 104.162860;
      geom.centralZ = 81.887600;
      geom.pRotationAngleX = -0.913765;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1326 based Row/Col = 132/6
      geom_tower geom;
      geom.id = 1326;
      geom.pDz = 6.923155;
      geom.pTheta = 0.030259;
      geom.pPhi = -0.290244;
      geom.pAlp1 = 0.022500;
      geom.pAlp2 = 0.022493;
      geom.pDy1 = 1.129974;
      geom.pDx1 = 1.189515;
      geom.pDx2 = 1.206492;
      geom.pDy2 = 1.248374;
      geom.pDx3 = 1.322554;
      geom.pDx4 = 1.341304;
      geom.centralX = 3.800504;
      geom.centralY = 104.162860;
      geom.centralZ = 81.887600;
      geom.pRotationAngleX = -0.913765;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1327 based Row/Col = 132/7
      geom_tower geom;
      geom.id = 1327;
      geom.pDz = 6.923155;
      geom.pTheta = 0.049092;
      geom.pPhi = -0.177238;
      geom.pAlp1 = 0.037543;
      geom.pAlp2 = 0.037532;
      geom.pDy1 = 1.129974;
      geom.pDx1 = 1.191275;
      geom.pDx2 = 1.208325;
      geom.pDy2 = 1.248374;
      geom.pDx3 = 1.324510;
      geom.pDx4 = 1.343341;
      geom.centralX = 6.337333;
      geom.centralY = 104.162860;
      geom.centralZ = 81.887600;
      geom.pRotationAngleX = -0.913765;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1328 based Row/Col = 132/8
      geom_tower geom;
      geom.id = 1328;
      geom.pDz = 6.923155;
      geom.pTheta = 0.068202;
      geom.pPhi = -0.127156;
      geom.pAlp1 = 0.052649;
      geom.pAlp2 = 0.052634;
      geom.pDy1 = 1.129974;
      geom.pDx1 = 1.193922;
      geom.pDx2 = 1.211081;
      geom.pDy2 = 1.248374;
      geom.pDx3 = 1.327452;
      geom.pDx4 = 1.346405;
      geom.centralX = 8.878911;
      geom.centralY = 104.162860;
      geom.centralZ = 81.887600;
      geom.pRotationAngleX = -0.913765;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1311 based Row/Col = 131/1
      geom_tower geom;
      geom.id = 1311;
      geom.pDz = 6.923089;
      geom.pTheta = 0.069063;
      geom.pPhi = 3.019280;
      geom.pAlp1 = -0.052648;
      geom.pAlp2 = -0.052633;
      geom.pDy1 = 1.129974;
      geom.pDx1 = 1.211081;
      geom.pDx2 = 1.228245;
      geom.pDy2 = 1.248374;
      geom.pDx3 = 1.346405;
      geom.pDx4 = 1.365362;
      geom.centralX = -9.004224;
      geom.centralY = 105.618534;
      geom.centralZ = 80.000443;
      geom.pRotationAngleX = -0.913765;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1312 based Row/Col = 131/2
      geom_tower geom;
      geom.id = 1312;
      geom.pDz = 6.923089;
      geom.pTheta = 0.049683;
      geom.pPhi = 2.971037;
      geom.pAlp1 = -0.037540;
      geom.pAlp2 = -0.037528;
      geom.pDy1 = 1.129974;
      geom.pDx1 = 1.208325;
      geom.pDx2 = 1.225375;
      geom.pDy2 = 1.248374;
      geom.pDx3 = 1.343341;
      geom.pDx4 = 1.362173;
      geom.centralX = -6.426647;
      geom.centralY = 105.618534;
      geom.centralZ = 80.000443;
      geom.pRotationAngleX = -0.913765;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1313 based Row/Col = 131/3
      geom_tower geom;
      geom.id = 1313;
      geom.pDz = 6.923089;
      geom.pTheta = 0.030565;
      geom.pPhi = 2.861924;
      geom.pAlp1 = -0.022497;
      geom.pAlp2 = -0.022491;
      geom.pDy1 = 1.129974;
      geom.pDx1 = 1.206492;
      geom.pDx2 = 1.223467;
      geom.pDy2 = 1.248374;
      geom.pDx3 = 1.341304;
      geom.pDx4 = 1.360053;
      geom.centralX = -3.854015;
      geom.centralY = 105.618534;
      geom.centralZ = 80.000443;
      geom.pRotationAngleX = -0.913765;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1314 based Row/Col = 131/4
      geom_tower geom;
      geom.id = 1314;
      geom.pDz = 6.923089;
      geom.pTheta = 0.012927;
      geom.pPhi = 2.430284;
      geom.pAlp1 = -0.007495;
      geom.pAlp2 = -0.007492;
      geom.pDy1 = 1.129974;
      geom.pDx1 = 1.205577;
      geom.pDx2 = 1.222515;
      geom.pDy2 = 1.248374;
      geom.pDx3 = 1.340287;
      geom.pDx4 = 1.358994;
      geom.centralX = -1.284343;
      geom.centralY = 105.618534;
      geom.centralZ = 80.000443;
      geom.pRotationAngleX = -0.913765;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1315 based Row/Col = 131/5
      geom_tower geom;
      geom.id = 1315;
      geom.pDz = 6.923089;
      geom.pTheta = 0.012927;
      geom.pPhi = 0.711309;
      geom.pAlp1 = 0.007495;
      geom.pAlp2 = 0.007492;
      geom.pDy1 = 1.129974;
      geom.pDx1 = 1.205577;
      geom.pDx2 = 1.222515;
      geom.pDy2 = 1.248374;
      geom.pDx3 = 1.340287;
      geom.pDx4 = 1.358994;
      geom.centralX = 1.284343;
      geom.centralY = 105.618534;
      geom.centralZ = 80.000443;
      geom.pRotationAngleX = -0.913765;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1316 based Row/Col = 131/6
      geom_tower geom;
      geom.id = 1316;
      geom.pDz = 6.923089;
      geom.pTheta = 0.030565;
      geom.pPhi = 0.279668;
      geom.pAlp1 = 0.022497;
      geom.pAlp2 = 0.022491;
      geom.pDy1 = 1.129974;
      geom.pDx1 = 1.206492;
      geom.pDx2 = 1.223467;
      geom.pDy2 = 1.248374;
      geom.pDx3 = 1.341304;
      geom.pDx4 = 1.360053;
      geom.centralX = 3.854015;
      geom.centralY = 105.618534;
      geom.centralZ = 80.000443;
      geom.pRotationAngleX = -0.913765;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1317 based Row/Col = 131/7
      geom_tower geom;
      geom.id = 1317;
      geom.pDz = 6.923089;
      geom.pTheta = 0.049683;
      geom.pPhi = 0.170556;
      geom.pAlp1 = 0.037540;
      geom.pAlp2 = 0.037528;
      geom.pDy1 = 1.129974;
      geom.pDx1 = 1.208325;
      geom.pDx2 = 1.225375;
      geom.pDy2 = 1.248374;
      geom.pDx3 = 1.343341;
      geom.pDx4 = 1.362173;
      geom.centralX = 6.426647;
      geom.centralY = 105.618534;
      geom.centralZ = 80.000443;
      geom.pRotationAngleX = -0.913765;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1318 based Row/Col = 131/8
      geom_tower geom;
      geom.id = 1318;
      geom.pDz = 6.923089;
      geom.pTheta = 0.069063;
      geom.pPhi = 0.122313;
      geom.pAlp1 = 0.052648;
      geom.pAlp2 = 0.052633;
      geom.pDy1 = 1.129974;
      geom.pDx1 = 1.211081;
      geom.pDx2 = 1.228245;
      geom.pDy2 = 1.248374;
      geom.pDx3 = 1.346405;
      geom.pDx4 = 1.365362;
      geom.centralX = 9.004224;
      geom.centralY = 105.618534;
      geom.centralZ = 80.000443;
      geom.pRotationAngleX = -0.913765;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 691 based Row/Col = 69/1
      geom_tower geom;
      geom.id = 691;
      geom.pDz = 6.923155;
      geom.pTheta = 0.068202;
      geom.pPhi = 3.014437;
      geom.pAlp1 = 0.052649;
      geom.pAlp2 = 0.052634;
      geom.pDy1 = 1.129974;
      geom.pDx1 = 1.211081;
      geom.pDx2 = 1.193922;
      geom.pDy2 = 1.248374;
      geom.pDx3 = 1.346405;
      geom.pDx4 = 1.327452;
      geom.centralX = -8.878911;
      geom.centralY = 104.162860;
      geom.centralZ = -81.887600;
      geom.pRotationAngleX = -2.227827;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 692 based Row/Col = 69/2
      geom_tower geom;
      geom.id = 692;
      geom.pDz = 6.923155;
      geom.pTheta = 0.049092;
      geom.pPhi = 2.964354;
      geom.pAlp1 = 0.037543;
      geom.pAlp2 = 0.037532;
      geom.pDy1 = 1.129974;
      geom.pDx1 = 1.208325;
      geom.pDx2 = 1.191275;
      geom.pDy2 = 1.248374;
      geom.pDx3 = 1.343341;
      geom.pDx4 = 1.324510;
      geom.centralX = -6.337333;
      geom.centralY = 104.162860;
      geom.centralZ = -81.887600;
      geom.pRotationAngleX = -2.227827;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 693 based Row/Col = 69/3
      geom_tower geom;
      geom.id = 693;
      geom.pDz = 6.923155;
      geom.pTheta = 0.030259;
      geom.pPhi = 2.851349;
      geom.pAlp1 = 0.022500;
      geom.pAlp2 = 0.022493;
      geom.pDy1 = 1.129974;
      geom.pDx1 = 1.206492;
      geom.pDx2 = 1.189515;
      geom.pDy2 = 1.248374;
      geom.pDx3 = 1.341304;
      geom.pDx4 = 1.322554;
      geom.centralX = -3.800504;
      geom.centralY = 104.162860;
      geom.centralZ = -81.887600;
      geom.pRotationAngleX = -2.227827;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 694 based Row/Col = 69/4
      geom_tower geom;
      geom.id = 694;
      geom.pDz = 6.923155;
      geom.pTheta = 0.012978;
      geom.pPhi = 2.410848;
      geom.pAlp1 = 0.007496;
      geom.pAlp2 = 0.007494;
      geom.pDy1 = 1.129974;
      geom.pDx1 = 1.205577;
      geom.pDx2 = 1.188636;
      geom.pDy2 = 1.248374;
      geom.pDx3 = 1.340287;
      geom.pDx4 = 1.321577;
      geom.centralX = -1.266519;
      geom.centralY = 104.162860;
      geom.centralZ = -81.887600;
      geom.pRotationAngleX = -2.227827;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 695 based Row/Col = 69/5
      geom_tower geom;
      geom.id = 695;
      geom.pDz = 6.923155;
      geom.pTheta = 0.012978;
      geom.pPhi = 0.730745;
      geom.pAlp1 = -0.007496;
      geom.pAlp2 = -0.007494;
      geom.pDy1 = 1.129974;
      geom.pDx1 = 1.205577;
      geom.pDx2 = 1.188636;
      geom.pDy2 = 1.248374;
      geom.pDx3 = 1.340287;
      geom.pDx4 = 1.321577;
      geom.centralX = 1.266519;
      geom.centralY = 104.162860;
      geom.centralZ = -81.887600;
      geom.pRotationAngleX = -2.227827;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 696 based Row/Col = 69/6
      geom_tower geom;
      geom.id = 696;
      geom.pDz = 6.923155;
      geom.pTheta = 0.030259;
      geom.pPhi = 0.290244;
      geom.pAlp1 = -0.022500;
      geom.pAlp2 = -0.022493;
      geom.pDy1 = 1.129974;
      geom.pDx1 = 1.206492;
      geom.pDx2 = 1.189515;
      geom.pDy2 = 1.248374;
      geom.pDx3 = 1.341304;
      geom.pDx4 = 1.322554;
      geom.centralX = 3.800504;
      geom.centralY = 104.162860;
      geom.centralZ = -81.887600;
      geom.pRotationAngleX = -2.227827;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 697 based Row/Col = 69/7
      geom_tower geom;
      geom.id = 697;
      geom.pDz = 6.923155;
      geom.pTheta = 0.049092;
      geom.pPhi = 0.177238;
      geom.pAlp1 = -0.037543;
      geom.pAlp2 = -0.037532;
      geom.pDy1 = 1.129974;
      geom.pDx1 = 1.208325;
      geom.pDx2 = 1.191275;
      geom.pDy2 = 1.248374;
      geom.pDx3 = 1.343341;
      geom.pDx4 = 1.324510;
      geom.centralX = 6.337333;
      geom.centralY = 104.162860;
      geom.centralZ = -81.887600;
      geom.pRotationAngleX = -2.227827;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 698 based Row/Col = 69/8
      geom_tower geom;
      geom.id = 698;
      geom.pDz = 6.923155;
      geom.pTheta = 0.068202;
      geom.pPhi = 0.127156;
      geom.pAlp1 = -0.052649;
      geom.pAlp2 = -0.052634;
      geom.pDy1 = 1.129974;
      geom.pDx1 = 1.211081;
      geom.pDx2 = 1.193922;
      geom.pDy2 = 1.248374;
      geom.pDx3 = 1.346405;
      geom.pDx4 = 1.327452;
      geom.centralX = 8.878911;
      geom.centralY = 104.162860;
      geom.centralZ = -81.887600;
      geom.pRotationAngleX = -2.227827;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 701 based Row/Col = 70/1
      geom_tower geom;
      geom.id = 701;
      geom.pDz = 6.923089;
      geom.pTheta = 0.069063;
      geom.pPhi = -3.019280;
      geom.pAlp1 = 0.052648;
      geom.pAlp2 = 0.052633;
      geom.pDy1 = 1.129974;
      geom.pDx1 = 1.228245;
      geom.pDx2 = 1.211081;
      geom.pDy2 = 1.248374;
      geom.pDx3 = 1.365362;
      geom.pDx4 = 1.346405;
      geom.centralX = -9.004224;
      geom.centralY = 105.618534;
      geom.centralZ = -80.000443;
      geom.pRotationAngleX = -2.227827;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 702 based Row/Col = 70/2
      geom_tower geom;
      geom.id = 702;
      geom.pDz = 6.923089;
      geom.pTheta = 0.049683;
      geom.pPhi = -2.971037;
      geom.pAlp1 = 0.037540;
      geom.pAlp2 = 0.037528;
      geom.pDy1 = 1.129974;
      geom.pDx1 = 1.225375;
      geom.pDx2 = 1.208325;
      geom.pDy2 = 1.248374;
      geom.pDx3 = 1.362173;
      geom.pDx4 = 1.343341;
      geom.centralX = -6.426647;
      geom.centralY = 105.618534;
      geom.centralZ = -80.000443;
      geom.pRotationAngleX = -2.227827;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 703 based Row/Col = 70/3
      geom_tower geom;
      geom.id = 703;
      geom.pDz = 6.923089;
      geom.pTheta = 0.030565;
      geom.pPhi = -2.861924;
      geom.pAlp1 = 0.022497;
      geom.pAlp2 = 0.022491;
      geom.pDy1 = 1.129974;
      geom.pDx1 = 1.223467;
      geom.pDx2 = 1.206492;
      geom.pDy2 = 1.248374;
      geom.pDx3 = 1.360053;
      geom.pDx4 = 1.341304;
      geom.centralX = -3.854015;
      geom.centralY = 105.618534;
      geom.centralZ = -80.000443;
      geom.pRotationAngleX = -2.227827;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 704 based Row/Col = 70/4
      geom_tower geom;
      geom.id = 704;
      geom.pDz = 6.923089;
      geom.pTheta = 0.012927;
      geom.pPhi = -2.430284;
      geom.pAlp1 = 0.007495;
      geom.pAlp2 = 0.007492;
      geom.pDy1 = 1.129974;
      geom.pDx1 = 1.222515;
      geom.pDx2 = 1.205577;
      geom.pDy2 = 1.248374;
      geom.pDx3 = 1.358994;
      geom.pDx4 = 1.340287;
      geom.centralX = -1.284343;
      geom.centralY = 105.618534;
      geom.centralZ = -80.000443;
      geom.pRotationAngleX = -2.227827;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 705 based Row/Col = 70/5
      geom_tower geom;
      geom.id = 705;
      geom.pDz = 6.923089;
      geom.pTheta = 0.012927;
      geom.pPhi = -0.711309;
      geom.pAlp1 = -0.007495;
      geom.pAlp2 = -0.007492;
      geom.pDy1 = 1.129974;
      geom.pDx1 = 1.222515;
      geom.pDx2 = 1.205577;
      geom.pDy2 = 1.248374;
      geom.pDx3 = 1.358994;
      geom.pDx4 = 1.340287;
      geom.centralX = 1.284343;
      geom.centralY = 105.618534;
      geom.centralZ = -80.000443;
      geom.pRotationAngleX = -2.227827;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 706 based Row/Col = 70/6
      geom_tower geom;
      geom.id = 706;
      geom.pDz = 6.923089;
      geom.pTheta = 0.030565;
      geom.pPhi = -0.279668;
      geom.pAlp1 = -0.022497;
      geom.pAlp2 = -0.022491;
      geom.pDy1 = 1.129974;
      geom.pDx1 = 1.223467;
      geom.pDx2 = 1.206492;
      geom.pDy2 = 1.248374;
      geom.pDx3 = 1.360053;
      geom.pDx4 = 1.341304;
      geom.centralX = 3.854015;
      geom.centralY = 105.618534;
      geom.centralZ = -80.000443;
      geom.pRotationAngleX = -2.227827;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 707 based Row/Col = 70/7
      geom_tower geom;
      geom.id = 707;
      geom.pDz = 6.923089;
      geom.pTheta = 0.049683;
      geom.pPhi = -0.170556;
      geom.pAlp1 = -0.037540;
      geom.pAlp2 = -0.037528;
      geom.pDy1 = 1.129974;
      geom.pDx1 = 1.225375;
      geom.pDx2 = 1.208325;
      geom.pDy2 = 1.248374;
      geom.pDx3 = 1.362173;
      geom.pDx4 = 1.343341;
      geom.centralX = 6.426647;
      geom.centralY = 105.618534;
      geom.centralZ = -80.000443;
      geom.pRotationAngleX = -2.227827;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 708 based Row/Col = 70/8
      geom_tower geom;
      geom.id = 708;
      geom.pDz = 6.923089;
      geom.pTheta = 0.069063;
      geom.pPhi = -0.122313;
      geom.pAlp1 = -0.052648;
      geom.pAlp2 = -0.052633;
      geom.pDy1 = 1.129974;
      geom.pDx1 = 1.228245;
      geom.pDx2 = 1.211081;
      geom.pDy2 = 1.248374;
      geom.pDx3 = 1.365362;
      geom.pDx4 = 1.346405;
      geom.centralX = 9.004224;
      geom.centralY = 105.618534;
      geom.centralZ = -80.000443;
      geom.pRotationAngleX = -2.227827;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1361 based Row/Col = 136/1
      geom_tower geom;
      geom.id = 1361;
      geom.pDz = 7.177098;
      geom.pTheta = 0.064336;
      geom.pPhi = -3.013671;
      geom.pAlp1 = -0.057289;
      geom.pAlp2 = -0.057275;
      geom.pDy1 = 1.128109;
      geom.pDx1 = 1.193449;
      geom.pDx2 = 1.212076;
      geom.pDy2 = 1.244514;
      geom.pDx3 = 1.323869;
      geom.pDx4 = 1.344414;
      geom.centralX = -8.873025;
      geom.centralY = 104.090279;
      geom.centralZ = 94.212410;
      geom.pRotationAngleX = -0.843786;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1362 based Row/Col = 136/2
      geom_tower geom;
      geom.id = 1362;
      geom.pDz = 7.177098;
      geom.pTheta = 0.046313;
      geom.pPhi = -2.963313;
      geom.pAlp1 = -0.040864;
      geom.pAlp2 = -0.040854;
      geom.pDy1 = 1.128109;
      geom.pDx1 = 1.191098;
      geom.pDx2 = 1.209619;
      geom.pDy2 = 1.244514;
      geom.pDx3 = 1.321262;
      geom.pDx4 = 1.341689;
      geom.centralX = -6.333656;
      geom.centralY = 104.090279;
      geom.centralZ = 94.212410;
      geom.pRotationAngleX = -0.843786;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1363 based Row/Col = 136/3
      geom_tower geom;
      geom.id = 1363;
      geom.pDz = 7.177098;
      geom.pTheta = 0.028555;
      geom.pPhi = -2.849718;
      geom.pAlp1 = -0.024496;
      geom.pAlp2 = -0.024489;
      geom.pDy1 = 1.128109;
      geom.pDx1 = 1.189534;
      geom.pDx2 = 1.207984;
      geom.pDy2 = 1.244514;
      geom.pDx3 = 1.319528;
      geom.pDx4 = 1.339877;
      geom.centralX = -3.798508;
      geom.centralY = 104.090279;
      geom.centralZ = 94.212410;
      geom.pRotationAngleX = -0.843786;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1364 based Row/Col = 136/4
      geom_tower geom;
      geom.id = 1364;
      geom.pDz = 7.177098;
      geom.pTheta = 0.012274;
      geom.pPhi = -2.407913;
      geom.pAlp1 = -0.008161;
      geom.pAlp2 = -0.008159;
      geom.pDy1 = 1.128109;
      geom.pDx1 = 1.188753;
      geom.pDx2 = 1.207168;
      geom.pDy2 = 1.244514;
      geom.pDx3 = 1.318662;
      geom.pDx4 = 1.338971;
      geom.centralX = -1.265889;
      geom.centralY = 104.090279;
      geom.centralZ = 94.212410;
      geom.pRotationAngleX = -0.843786;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1365 based Row/Col = 136/5
      geom_tower geom;
      geom.id = 1365;
      geom.pDz = 7.177098;
      geom.pTheta = 0.012274;
      geom.pPhi = -0.733679;
      geom.pAlp1 = 0.008161;
      geom.pAlp2 = 0.008159;
      geom.pDy1 = 1.128109;
      geom.pDx1 = 1.188753;
      geom.pDx2 = 1.207168;
      geom.pDy2 = 1.244514;
      geom.pDx3 = 1.318662;
      geom.pDx4 = 1.338971;
      geom.centralX = 1.265889;
      geom.centralY = 104.090279;
      geom.centralZ = 94.212410;
      geom.pRotationAngleX = -0.843786;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1366 based Row/Col = 136/6
      geom_tower geom;
      geom.id = 1366;
      geom.pDz = 7.177098;
      geom.pTheta = 0.028555;
      geom.pPhi = -0.291874;
      geom.pAlp1 = 0.024496;
      geom.pAlp2 = 0.024489;
      geom.pDy1 = 1.128109;
      geom.pDx1 = 1.189534;
      geom.pDx2 = 1.207984;
      geom.pDy2 = 1.244514;
      geom.pDx3 = 1.319528;
      geom.pDx4 = 1.339877;
      geom.centralX = 3.798508;
      geom.centralY = 104.090279;
      geom.centralZ = 94.212410;
      geom.pRotationAngleX = -0.843786;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1367 based Row/Col = 136/7
      geom_tower geom;
      geom.id = 1367;
      geom.pDz = 7.177098;
      geom.pTheta = 0.046313;
      geom.pPhi = -0.178280;
      geom.pAlp1 = 0.040864;
      geom.pAlp2 = 0.040854;
      geom.pDy1 = 1.128109;
      geom.pDx1 = 1.191098;
      geom.pDx2 = 1.209619;
      geom.pDy2 = 1.244514;
      geom.pDx3 = 1.321262;
      geom.pDx4 = 1.341689;
      geom.centralX = 6.333656;
      geom.centralY = 104.090279;
      geom.centralZ = 94.212410;
      geom.pRotationAngleX = -0.843786;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1368 based Row/Col = 136/8
      geom_tower geom;
      geom.id = 1368;
      geom.pDz = 7.177098;
      geom.pTheta = 0.064336;
      geom.pPhi = -0.127921;
      geom.pAlp1 = 0.057289;
      geom.pAlp2 = 0.057275;
      geom.pDy1 = 1.128109;
      geom.pDx1 = 1.193449;
      geom.pDx2 = 1.212076;
      geom.pDy2 = 1.244514;
      geom.pDx3 = 1.323869;
      geom.pDx4 = 1.344414;
      geom.centralX = 8.873025;
      geom.centralY = 104.090279;
      geom.centralZ = 94.212410;
      geom.pRotationAngleX = -0.843786;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1351 based Row/Col = 135/1
      geom_tower geom;
      geom.id = 1351;
      geom.pDz = 7.177066;
      geom.pTheta = 0.065226;
      geom.pPhi = 3.018806;
      geom.pAlp1 = -0.057287;
      geom.pAlp2 = -0.057273;
      geom.pDy1 = 1.128109;
      geom.pDx1 = 1.212076;
      geom.pDx2 = 1.230707;
      geom.pDy2 = 1.244514;
      geom.pDx3 = 1.344415;
      geom.pDx4 = 1.364964;
      geom.centralX = -9.009079;
      geom.centralY = 105.670541;
      geom.centralZ = 92.435934;
      geom.pRotationAngleX = -0.843786;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1352 based Row/Col = 135/2
      geom_tower geom;
      geom.id = 1352;
      geom.pDz = 7.177066;
      geom.pTheta = 0.046925;
      geom.pPhi = 2.970396;
      geom.pAlp1 = -0.040860;
      geom.pAlp2 = -0.040850;
      geom.pDy1 = 1.128109;
      geom.pDx1 = 1.209619;
      geom.pDx2 = 1.228140;
      geom.pDy2 = 1.244514;
      geom.pDx3 = 1.341689;
      geom.pDx4 = 1.362117;
      geom.centralX = -6.430648;
      geom.centralY = 105.670541;
      geom.centralZ = 92.435934;
      geom.pRotationAngleX = -0.843786;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1353 based Row/Col = 135/3
      geom_tower geom;
      geom.id = 1353;
      geom.pDz = 7.177066;
      geom.pTheta = 0.028873;
      geom.pPhi = 2.860923;
      geom.pAlp1 = -0.024493;
      geom.pAlp2 = -0.024486;
      geom.pDy1 = 1.128109;
      geom.pDx1 = 1.207984;
      geom.pDx2 = 1.226432;
      geom.pDy2 = 1.244514;
      geom.pDx3 = 1.339877;
      geom.pDx4 = 1.360223;
      geom.centralX = -3.856627;
      geom.centralY = 105.670541;
      geom.centralZ = 92.435934;
      geom.pRotationAngleX = -0.843786;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1354 based Row/Col = 135/4
      geom_tower geom;
      geom.id = 1354;
      geom.pDz = 7.177066;
      geom.pTheta = 0.012228;
      geom.pPhi = 2.428434;
      geom.pAlp1 = -0.008160;
      geom.pAlp2 = -0.008158;
      geom.pDy1 = 1.128109;
      geom.pDx1 = 1.207168;
      geom.pDx2 = 1.225580;
      geom.pDy2 = 1.244514;
      geom.pDx3 = 1.338971;
      geom.pDx4 = 1.359278;
      geom.centralX = -1.285249;
      geom.centralY = 105.670541;
      geom.centralZ = 92.435934;
      geom.pRotationAngleX = -0.843786;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1355 based Row/Col = 135/5
      geom_tower geom;
      geom.id = 1355;
      geom.pDz = 7.177066;
      geom.pTheta = 0.012228;
      geom.pPhi = 0.713158;
      geom.pAlp1 = 0.008160;
      geom.pAlp2 = 0.008158;
      geom.pDy1 = 1.128109;
      geom.pDx1 = 1.207168;
      geom.pDx2 = 1.225580;
      geom.pDy2 = 1.244514;
      geom.pDx3 = 1.338971;
      geom.pDx4 = 1.359278;
      geom.centralX = 1.285249;
      geom.centralY = 105.670541;
      geom.centralZ = 92.435934;
      geom.pRotationAngleX = -0.843786;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1356 based Row/Col = 135/6
      geom_tower geom;
      geom.id = 1356;
      geom.pDz = 7.177066;
      geom.pTheta = 0.028873;
      geom.pPhi = 0.280670;
      geom.pAlp1 = 0.024493;
      geom.pAlp2 = 0.024486;
      geom.pDy1 = 1.128109;
      geom.pDx1 = 1.207984;
      geom.pDx2 = 1.226432;
      geom.pDy2 = 1.244514;
      geom.pDx3 = 1.339877;
      geom.pDx4 = 1.360223;
      geom.centralX = 3.856627;
      geom.centralY = 105.670541;
      geom.centralZ = 92.435934;
      geom.pRotationAngleX = -0.843786;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1357 based Row/Col = 135/7
      geom_tower geom;
      geom.id = 1357;
      geom.pDz = 7.177066;
      geom.pTheta = 0.046925;
      geom.pPhi = 0.171196;
      geom.pAlp1 = 0.040860;
      geom.pAlp2 = 0.040850;
      geom.pDy1 = 1.128109;
      geom.pDx1 = 1.209619;
      geom.pDx2 = 1.228140;
      geom.pDy2 = 1.244514;
      geom.pDx3 = 1.341689;
      geom.pDx4 = 1.362117;
      geom.centralX = 6.430648;
      geom.centralY = 105.670541;
      geom.centralZ = 92.435934;
      geom.pRotationAngleX = -0.843786;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1358 based Row/Col = 135/8
      geom_tower geom;
      geom.id = 1358;
      geom.pDz = 7.177066;
      geom.pTheta = 0.065226;
      geom.pPhi = 0.122786;
      geom.pAlp1 = 0.057287;
      geom.pAlp2 = 0.057273;
      geom.pDy1 = 1.128109;
      geom.pDx1 = 1.212076;
      geom.pDx2 = 1.230707;
      geom.pDy2 = 1.244514;
      geom.pDx3 = 1.344415;
      geom.pDx4 = 1.364964;
      geom.centralX = 9.009079;
      geom.centralY = 105.670541;
      geom.centralZ = 92.435934;
      geom.pRotationAngleX = -0.843786;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 651 based Row/Col = 65/1
      geom_tower geom;
      geom.id = 651;
      geom.pDz = 7.177098;
      geom.pTheta = 0.064336;
      geom.pPhi = 3.013671;
      geom.pAlp1 = 0.057289;
      geom.pAlp2 = 0.057275;
      geom.pDy1 = 1.128109;
      geom.pDx1 = 1.212076;
      geom.pDx2 = 1.193449;
      geom.pDy2 = 1.244514;
      geom.pDx3 = 1.344414;
      geom.pDx4 = 1.323869;
      geom.centralX = -8.873025;
      geom.centralY = 104.090279;
      geom.centralZ = -94.212410;
      geom.pRotationAngleX = -2.297807;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 652 based Row/Col = 65/2
      geom_tower geom;
      geom.id = 652;
      geom.pDz = 7.177098;
      geom.pTheta = 0.046313;
      geom.pPhi = 2.963313;
      geom.pAlp1 = 0.040864;
      geom.pAlp2 = 0.040854;
      geom.pDy1 = 1.128109;
      geom.pDx1 = 1.209619;
      geom.pDx2 = 1.191098;
      geom.pDy2 = 1.244514;
      geom.pDx3 = 1.341689;
      geom.pDx4 = 1.321262;
      geom.centralX = -6.333656;
      geom.centralY = 104.090279;
      geom.centralZ = -94.212410;
      geom.pRotationAngleX = -2.297807;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 653 based Row/Col = 65/3
      geom_tower geom;
      geom.id = 653;
      geom.pDz = 7.177098;
      geom.pTheta = 0.028555;
      geom.pPhi = 2.849718;
      geom.pAlp1 = 0.024496;
      geom.pAlp2 = 0.024489;
      geom.pDy1 = 1.128109;
      geom.pDx1 = 1.207984;
      geom.pDx2 = 1.189534;
      geom.pDy2 = 1.244514;
      geom.pDx3 = 1.339877;
      geom.pDx4 = 1.319528;
      geom.centralX = -3.798508;
      geom.centralY = 104.090279;
      geom.centralZ = -94.212410;
      geom.pRotationAngleX = -2.297807;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 654 based Row/Col = 65/4
      geom_tower geom;
      geom.id = 654;
      geom.pDz = 7.177098;
      geom.pTheta = 0.012274;
      geom.pPhi = 2.407913;
      geom.pAlp1 = 0.008161;
      geom.pAlp2 = 0.008159;
      geom.pDy1 = 1.128109;
      geom.pDx1 = 1.207168;
      geom.pDx2 = 1.188753;
      geom.pDy2 = 1.244514;
      geom.pDx3 = 1.338971;
      geom.pDx4 = 1.318662;
      geom.centralX = -1.265889;
      geom.centralY = 104.090279;
      geom.centralZ = -94.212410;
      geom.pRotationAngleX = -2.297807;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 655 based Row/Col = 65/5
      geom_tower geom;
      geom.id = 655;
      geom.pDz = 7.177098;
      geom.pTheta = 0.012274;
      geom.pPhi = 0.733679;
      geom.pAlp1 = -0.008161;
      geom.pAlp2 = -0.008159;
      geom.pDy1 = 1.128109;
      geom.pDx1 = 1.207168;
      geom.pDx2 = 1.188753;
      geom.pDy2 = 1.244514;
      geom.pDx3 = 1.338971;
      geom.pDx4 = 1.318662;
      geom.centralX = 1.265889;
      geom.centralY = 104.090279;
      geom.centralZ = -94.212410;
      geom.pRotationAngleX = -2.297807;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 656 based Row/Col = 65/6
      geom_tower geom;
      geom.id = 656;
      geom.pDz = 7.177098;
      geom.pTheta = 0.028555;
      geom.pPhi = 0.291874;
      geom.pAlp1 = -0.024496;
      geom.pAlp2 = -0.024489;
      geom.pDy1 = 1.128109;
      geom.pDx1 = 1.207984;
      geom.pDx2 = 1.189534;
      geom.pDy2 = 1.244514;
      geom.pDx3 = 1.339877;
      geom.pDx4 = 1.319528;
      geom.centralX = 3.798508;
      geom.centralY = 104.090279;
      geom.centralZ = -94.212410;
      geom.pRotationAngleX = -2.297807;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 657 based Row/Col = 65/7
      geom_tower geom;
      geom.id = 657;
      geom.pDz = 7.177098;
      geom.pTheta = 0.046313;
      geom.pPhi = 0.178280;
      geom.pAlp1 = -0.040864;
      geom.pAlp2 = -0.040854;
      geom.pDy1 = 1.128109;
      geom.pDx1 = 1.209619;
      geom.pDx2 = 1.191098;
      geom.pDy2 = 1.244514;
      geom.pDx3 = 1.341689;
      geom.pDx4 = 1.321262;
      geom.centralX = 6.333656;
      geom.centralY = 104.090279;
      geom.centralZ = -94.212410;
      geom.pRotationAngleX = -2.297807;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 658 based Row/Col = 65/8
      geom_tower geom;
      geom.id = 658;
      geom.pDz = 7.177098;
      geom.pTheta = 0.064336;
      geom.pPhi = 0.127921;
      geom.pAlp1 = -0.057289;
      geom.pAlp2 = -0.057275;
      geom.pDy1 = 1.128109;
      geom.pDx1 = 1.212076;
      geom.pDx2 = 1.193449;
      geom.pDy2 = 1.244514;
      geom.pDx3 = 1.344414;
      geom.pDx4 = 1.323869;
      geom.centralX = 8.873025;
      geom.centralY = 104.090279;
      geom.centralZ = -94.212410;
      geom.pRotationAngleX = -2.297807;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 661 based Row/Col = 66/1
      geom_tower geom;
      geom.id = 661;
      geom.pDz = 7.177066;
      geom.pTheta = 0.065226;
      geom.pPhi = -3.018806;
      geom.pAlp1 = 0.057287;
      geom.pAlp2 = 0.057273;
      geom.pDy1 = 1.128109;
      geom.pDx1 = 1.230707;
      geom.pDx2 = 1.212076;
      geom.pDy2 = 1.244514;
      geom.pDx3 = 1.364964;
      geom.pDx4 = 1.344415;
      geom.centralX = -9.009079;
      geom.centralY = 105.670541;
      geom.centralZ = -92.435934;
      geom.pRotationAngleX = -2.297807;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 662 based Row/Col = 66/2
      geom_tower geom;
      geom.id = 662;
      geom.pDz = 7.177066;
      geom.pTheta = 0.046925;
      geom.pPhi = -2.970396;
      geom.pAlp1 = 0.040860;
      geom.pAlp2 = 0.040850;
      geom.pDy1 = 1.128109;
      geom.pDx1 = 1.228140;
      geom.pDx2 = 1.209619;
      geom.pDy2 = 1.244514;
      geom.pDx3 = 1.362117;
      geom.pDx4 = 1.341689;
      geom.centralX = -6.430648;
      geom.centralY = 105.670541;
      geom.centralZ = -92.435934;
      geom.pRotationAngleX = -2.297807;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 663 based Row/Col = 66/3
      geom_tower geom;
      geom.id = 663;
      geom.pDz = 7.177066;
      geom.pTheta = 0.028873;
      geom.pPhi = -2.860923;
      geom.pAlp1 = 0.024493;
      geom.pAlp2 = 0.024486;
      geom.pDy1 = 1.128109;
      geom.pDx1 = 1.226432;
      geom.pDx2 = 1.207984;
      geom.pDy2 = 1.244514;
      geom.pDx3 = 1.360223;
      geom.pDx4 = 1.339877;
      geom.centralX = -3.856627;
      geom.centralY = 105.670541;
      geom.centralZ = -92.435934;
      geom.pRotationAngleX = -2.297807;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 664 based Row/Col = 66/4
      geom_tower geom;
      geom.id = 664;
      geom.pDz = 7.177066;
      geom.pTheta = 0.012228;
      geom.pPhi = -2.428434;
      geom.pAlp1 = 0.008160;
      geom.pAlp2 = 0.008158;
      geom.pDy1 = 1.128109;
      geom.pDx1 = 1.225580;
      geom.pDx2 = 1.207168;
      geom.pDy2 = 1.244514;
      geom.pDx3 = 1.359278;
      geom.pDx4 = 1.338971;
      geom.centralX = -1.285249;
      geom.centralY = 105.670541;
      geom.centralZ = -92.435934;
      geom.pRotationAngleX = -2.297807;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 665 based Row/Col = 66/5
      geom_tower geom;
      geom.id = 665;
      geom.pDz = 7.177066;
      geom.pTheta = 0.012228;
      geom.pPhi = -0.713158;
      geom.pAlp1 = -0.008160;
      geom.pAlp2 = -0.008158;
      geom.pDy1 = 1.128109;
      geom.pDx1 = 1.225580;
      geom.pDx2 = 1.207168;
      geom.pDy2 = 1.244514;
      geom.pDx3 = 1.359278;
      geom.pDx4 = 1.338971;
      geom.centralX = 1.285249;
      geom.centralY = 105.670541;
      geom.centralZ = -92.435934;
      geom.pRotationAngleX = -2.297807;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 666 based Row/Col = 66/6
      geom_tower geom;
      geom.id = 666;
      geom.pDz = 7.177066;
      geom.pTheta = 0.028873;
      geom.pPhi = -0.280670;
      geom.pAlp1 = -0.024493;
      geom.pAlp2 = -0.024486;
      geom.pDy1 = 1.128109;
      geom.pDx1 = 1.226432;
      geom.pDx2 = 1.207984;
      geom.pDy2 = 1.244514;
      geom.pDx3 = 1.360223;
      geom.pDx4 = 1.339877;
      geom.centralX = 3.856627;
      geom.centralY = 105.670541;
      geom.centralZ = -92.435934;
      geom.pRotationAngleX = -2.297807;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 667 based Row/Col = 66/7
      geom_tower geom;
      geom.id = 667;
      geom.pDz = 7.177066;
      geom.pTheta = 0.046925;
      geom.pPhi = -0.171196;
      geom.pAlp1 = -0.040860;
      geom.pAlp2 = -0.040850;
      geom.pDy1 = 1.128109;
      geom.pDx1 = 1.228140;
      geom.pDx2 = 1.209619;
      geom.pDy2 = 1.244514;
      geom.pDx3 = 1.362117;
      geom.pDx4 = 1.341689;
      geom.centralX = 6.430648;
      geom.centralY = 105.670541;
      geom.centralZ = -92.435934;
      geom.pRotationAngleX = -2.297807;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 668 based Row/Col = 66/8
      geom_tower geom;
      geom.id = 668;
      geom.pDz = 7.177066;
      geom.pTheta = 0.065226;
      geom.pPhi = -0.122786;
      geom.pAlp1 = -0.057287;
      geom.pAlp2 = -0.057273;
      geom.pDy1 = 1.128109;
      geom.pDx1 = 1.230707;
      geom.pDx2 = 1.212076;
      geom.pDy2 = 1.244514;
      geom.pDx3 = 1.364964;
      geom.pDx4 = 1.344415;
      geom.centralX = 9.009079;
      geom.centralY = 105.670541;
      geom.centralZ = -92.435934;
      geom.pRotationAngleX = -2.297807;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1401 based Row/Col = 140/1
      geom_tower geom;
      geom.id = 1401;
      geom.pDz = 7.507343;
      geom.pTheta = 0.060429;
      geom.pPhi = -3.013392;
      geom.pAlp1 = -0.061383;
      geom.pAlp2 = -0.061370;
      geom.pDy1 = 1.125494;
      geom.pDx1 = 1.192753;
      geom.pDx2 = 1.212651;
      geom.pDy2 = 1.239457;
      geom.pDx3 = 1.320760;
      geom.pDx4 = 1.342668;
      geom.centralX = -8.867121;
      geom.centralY = 104.017699;
      geom.centralZ = 107.254167;
      geom.pRotationAngleX = -0.778199;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1402 based Row/Col = 140/2
      geom_tower geom;
      geom.id = 1402;
      geom.pDz = 7.507343;
      geom.pTheta = 0.043501;
      geom.pPhi = -2.962941;
      geom.pAlp1 = -0.043798;
      geom.pAlp2 = -0.043788;
      geom.pDy1 = 1.125494;
      geom.pDx1 = 1.190683;
      geom.pDx2 = 1.210480;
      geom.pDy2 = 1.239457;
      geom.pDx3 = 1.318468;
      geom.pDx4 = 1.340265;
      geom.centralX = -6.329939;
      geom.centralY = 104.017699;
      geom.centralZ = 107.254167;
      geom.pRotationAngleX = -0.778199;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1403 based Row/Col = 140/3
      geom_tower geom;
      geom.id = 1403;
      geom.pDz = 7.507343;
      geom.pTheta = 0.026824;
      geom.pPhi = -2.849146;
      geom.pAlp1 = -0.026260;
      geom.pAlp2 = -0.026254;
      geom.pDy1 = 1.125494;
      geom.pDx1 = 1.189305;
      geom.pDx2 = 1.209035;
      geom.pDy2 = 1.239457;
      geom.pDx3 = 1.316943;
      geom.pDx4 = 1.338666;
      geom.centralX = -3.796477;
      geom.centralY = 104.017699;
      geom.centralZ = 107.254167;
      geom.pRotationAngleX = -0.778199;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1404 based Row/Col = 140/4
      geom_tower geom;
      geom.id = 1404;
      geom.pDz = 7.507343;
      geom.pTheta = 0.011538;
      geom.pPhi = -2.406895;
      geom.pAlp1 = -0.008750;
      geom.pAlp2 = -0.008748;
      geom.pDy1 = 1.125494;
      geom.pDx1 = 1.188617;
      geom.pDx2 = 1.208314;
      geom.pDy2 = 1.239457;
      geom.pDx3 = 1.316181;
      geom.pDx4 = 1.337867;
      geom.centralX = -1.265245;
      geom.centralY = 104.017699;
      geom.centralZ = 107.254167;
      geom.pRotationAngleX = -0.778199;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1405 based Row/Col = 140/5
      geom_tower geom;
      geom.id = 1405;
      geom.pDz = 7.507343;
      geom.pTheta = 0.011538;
      geom.pPhi = -0.734698;
      geom.pAlp1 = 0.008750;
      geom.pAlp2 = 0.008748;
      geom.pDy1 = 1.125494;
      geom.pDx1 = 1.188617;
      geom.pDx2 = 1.208314;
      geom.pDy2 = 1.239457;
      geom.pDx3 = 1.316181;
      geom.pDx4 = 1.337867;
      geom.centralX = 1.265245;
      geom.centralY = 104.017699;
      geom.centralZ = 107.254167;
      geom.pRotationAngleX = -0.778199;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1406 based Row/Col = 140/6
      geom_tower geom;
      geom.id = 1406;
      geom.pDz = 7.507343;
      geom.pTheta = 0.026824;
      geom.pPhi = -0.292446;
      geom.pAlp1 = 0.026260;
      geom.pAlp2 = 0.026254;
      geom.pDy1 = 1.125494;
      geom.pDx1 = 1.189305;
      geom.pDx2 = 1.209035;
      geom.pDy2 = 1.239457;
      geom.pDx3 = 1.316943;
      geom.pDx4 = 1.338666;
      geom.centralX = 3.796477;
      geom.centralY = 104.017699;
      geom.centralZ = 107.254167;
      geom.pRotationAngleX = -0.778199;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1407 based Row/Col = 140/7
      geom_tower geom;
      geom.id = 1407;
      geom.pDz = 7.507343;
      geom.pTheta = 0.043501;
      geom.pPhi = -0.178651;
      geom.pAlp1 = 0.043798;
      geom.pAlp2 = 0.043788;
      geom.pDy1 = 1.125494;
      geom.pDx1 = 1.190683;
      geom.pDx2 = 1.210480;
      geom.pDy2 = 1.239457;
      geom.pDx3 = 1.318468;
      geom.pDx4 = 1.340265;
      geom.centralX = 6.329939;
      geom.centralY = 104.017699;
      geom.centralZ = 107.254167;
      geom.pRotationAngleX = -0.778199;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1408 based Row/Col = 140/8
      geom_tower geom;
      geom.id = 1408;
      geom.pDz = 7.507343;
      geom.pTheta = 0.060429;
      geom.pPhi = -0.128201;
      geom.pAlp1 = 0.061383;
      geom.pAlp2 = 0.061370;
      geom.pDy1 = 1.125494;
      geom.pDx1 = 1.192753;
      geom.pDx2 = 1.212651;
      geom.pDy2 = 1.239457;
      geom.pDx3 = 1.320760;
      geom.pDx4 = 1.342668;
      geom.centralX = 8.867121;
      geom.centralY = 104.017699;
      geom.centralZ = 107.254167;
      geom.pRotationAngleX = -0.778199;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1391 based Row/Col = 139/1
      geom_tower geom;
      geom.id = 1391;
      geom.pDz = 7.507331;
      geom.pTheta = 0.061313;
      geom.pPhi = 3.020020;
      geom.pAlp1 = -0.061381;
      geom.pAlp2 = -0.061368;
      geom.pDy1 = 1.125494;
      geom.pDx1 = 1.212651;
      geom.pDx2 = 1.232552;
      geom.pDy2 = 1.239457;
      geom.pDx3 = 1.342668;
      geom.pDx4 = 1.364581;
      geom.centralX = -9.012453;
      geom.centralY = 105.705529;
      geom.centralZ = 105.590467;
      geom.pRotationAngleX = -0.778199;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1392 based Row/Col = 139/2
      geom_tower geom;
      geom.id = 1392;
      geom.pDz = 7.507331;
      geom.pTheta = 0.044104;
      geom.pPhi = 2.972087;
      geom.pAlp1 = -0.043794;
      geom.pAlp2 = -0.043784;
      geom.pDy1 = 1.125494;
      geom.pDx1 = 1.210480;
      geom.pDx2 = 1.230278;
      geom.pDy2 = 1.239457;
      geom.pDx3 = 1.340265;
      geom.pDx4 = 1.362063;
      geom.centralX = -6.433568;
      geom.centralY = 105.705529;
      geom.centralZ = 105.590467;
      geom.pRotationAngleX = -0.778199;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1393 based Row/Col = 139/3
      geom_tower geom;
      geom.id = 1393;
      geom.pDz = 7.507331;
      geom.pTheta = 0.027124;
      geom.pPhi = 2.863618;
      geom.pAlp1 = -0.026257;
      geom.pAlp2 = -0.026251;
      geom.pDy1 = 1.125494;
      geom.pDx1 = 1.209035;
      geom.pDx2 = 1.228764;
      geom.pDy2 = 1.239457;
      geom.pDx3 = 1.338666;
      geom.pDx4 = 1.360387;
      geom.centralX = -3.858584;
      geom.centralY = 105.705529;
      geom.centralZ = 105.590467;
      geom.pRotationAngleX = -0.778199;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1394 based Row/Col = 139/4
      geom_tower geom;
      geom.id = 1394;
      geom.pDz = 7.507331;
      geom.pTheta = 0.011446;
      geom.pPhi = 2.433475;
      geom.pAlp1 = -0.008749;
      geom.pAlp2 = -0.008747;
      geom.pDy1 = 1.125494;
      geom.pDx1 = 1.208314;
      geom.pDx2 = 1.228008;
      geom.pDy2 = 1.239457;
      geom.pDx3 = 1.337867;
      geom.pDx4 = 1.359551;
      geom.centralX = -1.285935;
      geom.centralY = 105.705529;
      geom.centralZ = 105.590467;
      geom.pRotationAngleX = -0.778199;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1395 based Row/Col = 139/5
      geom_tower geom;
      geom.id = 1395;
      geom.pDz = 7.507331;
      geom.pTheta = 0.011446;
      geom.pPhi = 0.708117;
      geom.pAlp1 = 0.008749;
      geom.pAlp2 = 0.008747;
      geom.pDy1 = 1.125494;
      geom.pDx1 = 1.208314;
      geom.pDx2 = 1.228008;
      geom.pDy2 = 1.239457;
      geom.pDx3 = 1.337867;
      geom.pDx4 = 1.359551;
      geom.centralX = 1.285935;
      geom.centralY = 105.705529;
      geom.centralZ = 105.590467;
      geom.pRotationAngleX = -0.778199;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1396 based Row/Col = 139/6
      geom_tower geom;
      geom.id = 1396;
      geom.pDz = 7.507331;
      geom.pTheta = 0.027124;
      geom.pPhi = 0.277974;
      geom.pAlp1 = 0.026257;
      geom.pAlp2 = 0.026251;
      geom.pDy1 = 1.125494;
      geom.pDx1 = 1.209035;
      geom.pDx2 = 1.228764;
      geom.pDy2 = 1.239457;
      geom.pDx3 = 1.338666;
      geom.pDx4 = 1.360387;
      geom.centralX = 3.858584;
      geom.centralY = 105.705529;
      geom.centralZ = 105.590467;
      geom.pRotationAngleX = -0.778199;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1397 based Row/Col = 139/7
      geom_tower geom;
      geom.id = 1397;
      geom.pDz = 7.507331;
      geom.pTheta = 0.044104;
      geom.pPhi = 0.169506;
      geom.pAlp1 = 0.043794;
      geom.pAlp2 = 0.043784;
      geom.pDy1 = 1.125494;
      geom.pDx1 = 1.210480;
      geom.pDx2 = 1.230278;
      geom.pDy2 = 1.239457;
      geom.pDx3 = 1.340265;
      geom.pDx4 = 1.362063;
      geom.centralX = 6.433568;
      geom.centralY = 105.705529;
      geom.centralZ = 105.590467;
      geom.pRotationAngleX = -0.778199;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1398 based Row/Col = 139/8
      geom_tower geom;
      geom.id = 1398;
      geom.pDz = 7.507331;
      geom.pTheta = 0.061313;
      geom.pPhi = 0.121572;
      geom.pAlp1 = 0.061381;
      geom.pAlp2 = 0.061368;
      geom.pDy1 = 1.125494;
      geom.pDx1 = 1.212651;
      geom.pDx2 = 1.232552;
      geom.pDy2 = 1.239457;
      geom.pDx3 = 1.342668;
      geom.pDx4 = 1.364581;
      geom.centralX = 9.012453;
      geom.centralY = 105.705529;
      geom.centralZ = 105.590467;
      geom.pRotationAngleX = -0.778199;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 611 based Row/Col = 61/1
      geom_tower geom;
      geom.id = 611;
      geom.pDz = 7.507343;
      geom.pTheta = 0.060429;
      geom.pPhi = 3.013392;
      geom.pAlp1 = 0.061383;
      geom.pAlp2 = 0.061370;
      geom.pDy1 = 1.125494;
      geom.pDx1 = 1.212651;
      geom.pDx2 = 1.192753;
      geom.pDy2 = 1.239457;
      geom.pDx3 = 1.342668;
      geom.pDx4 = 1.320760;
      geom.centralX = -8.867121;
      geom.centralY = 104.017699;
      geom.centralZ = -107.254167;
      geom.pRotationAngleX = -2.363394;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 612 based Row/Col = 61/2
      geom_tower geom;
      geom.id = 612;
      geom.pDz = 7.507343;
      geom.pTheta = 0.043501;
      geom.pPhi = 2.962941;
      geom.pAlp1 = 0.043798;
      geom.pAlp2 = 0.043788;
      geom.pDy1 = 1.125494;
      geom.pDx1 = 1.210480;
      geom.pDx2 = 1.190683;
      geom.pDy2 = 1.239457;
      geom.pDx3 = 1.340265;
      geom.pDx4 = 1.318468;
      geom.centralX = -6.329939;
      geom.centralY = 104.017699;
      geom.centralZ = -107.254167;
      geom.pRotationAngleX = -2.363394;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 613 based Row/Col = 61/3
      geom_tower geom;
      geom.id = 613;
      geom.pDz = 7.507343;
      geom.pTheta = 0.026824;
      geom.pPhi = 2.849146;
      geom.pAlp1 = 0.026260;
      geom.pAlp2 = 0.026254;
      geom.pDy1 = 1.125494;
      geom.pDx1 = 1.209035;
      geom.pDx2 = 1.189305;
      geom.pDy2 = 1.239457;
      geom.pDx3 = 1.338666;
      geom.pDx4 = 1.316943;
      geom.centralX = -3.796477;
      geom.centralY = 104.017699;
      geom.centralZ = -107.254167;
      geom.pRotationAngleX = -2.363394;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 614 based Row/Col = 61/4
      geom_tower geom;
      geom.id = 614;
      geom.pDz = 7.507343;
      geom.pTheta = 0.011538;
      geom.pPhi = 2.406895;
      geom.pAlp1 = 0.008750;
      geom.pAlp2 = 0.008748;
      geom.pDy1 = 1.125494;
      geom.pDx1 = 1.208314;
      geom.pDx2 = 1.188617;
      geom.pDy2 = 1.239457;
      geom.pDx3 = 1.337867;
      geom.pDx4 = 1.316181;
      geom.centralX = -1.265245;
      geom.centralY = 104.017699;
      geom.centralZ = -107.254167;
      geom.pRotationAngleX = -2.363394;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 615 based Row/Col = 61/5
      geom_tower geom;
      geom.id = 615;
      geom.pDz = 7.507343;
      geom.pTheta = 0.011538;
      geom.pPhi = 0.734698;
      geom.pAlp1 = -0.008750;
      geom.pAlp2 = -0.008748;
      geom.pDy1 = 1.125494;
      geom.pDx1 = 1.208314;
      geom.pDx2 = 1.188617;
      geom.pDy2 = 1.239457;
      geom.pDx3 = 1.337867;
      geom.pDx4 = 1.316181;
      geom.centralX = 1.265245;
      geom.centralY = 104.017699;
      geom.centralZ = -107.254167;
      geom.pRotationAngleX = -2.363394;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 616 based Row/Col = 61/6
      geom_tower geom;
      geom.id = 616;
      geom.pDz = 7.507343;
      geom.pTheta = 0.026824;
      geom.pPhi = 0.292446;
      geom.pAlp1 = -0.026260;
      geom.pAlp2 = -0.026254;
      geom.pDy1 = 1.125494;
      geom.pDx1 = 1.209035;
      geom.pDx2 = 1.189305;
      geom.pDy2 = 1.239457;
      geom.pDx3 = 1.338666;
      geom.pDx4 = 1.316943;
      geom.centralX = 3.796477;
      geom.centralY = 104.017699;
      geom.centralZ = -107.254167;
      geom.pRotationAngleX = -2.363394;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 617 based Row/Col = 61/7
      geom_tower geom;
      geom.id = 617;
      geom.pDz = 7.507343;
      geom.pTheta = 0.043501;
      geom.pPhi = 0.178651;
      geom.pAlp1 = -0.043798;
      geom.pAlp2 = -0.043788;
      geom.pDy1 = 1.125494;
      geom.pDx1 = 1.210480;
      geom.pDx2 = 1.190683;
      geom.pDy2 = 1.239457;
      geom.pDx3 = 1.340265;
      geom.pDx4 = 1.318468;
      geom.centralX = 6.329939;
      geom.centralY = 104.017699;
      geom.centralZ = -107.254167;
      geom.pRotationAngleX = -2.363394;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 618 based Row/Col = 61/8
      geom_tower geom;
      geom.id = 618;
      geom.pDz = 7.507343;
      geom.pTheta = 0.060429;
      geom.pPhi = 0.128201;
      geom.pAlp1 = -0.061383;
      geom.pAlp2 = -0.061370;
      geom.pDy1 = 1.125494;
      geom.pDx1 = 1.212651;
      geom.pDx2 = 1.192753;
      geom.pDy2 = 1.239457;
      geom.pDx3 = 1.342668;
      geom.pDx4 = 1.320760;
      geom.centralX = 8.867121;
      geom.centralY = 104.017699;
      geom.centralZ = -107.254167;
      geom.pRotationAngleX = -2.363394;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 621 based Row/Col = 62/1
      geom_tower geom;
      geom.id = 621;
      geom.pDz = 7.507331;
      geom.pTheta = 0.061313;
      geom.pPhi = -3.020020;
      geom.pAlp1 = 0.061381;
      geom.pAlp2 = 0.061368;
      geom.pDy1 = 1.125494;
      geom.pDx1 = 1.232552;
      geom.pDx2 = 1.212651;
      geom.pDy2 = 1.239457;
      geom.pDx3 = 1.364581;
      geom.pDx4 = 1.342668;
      geom.centralX = -9.012453;
      geom.centralY = 105.705529;
      geom.centralZ = -105.590467;
      geom.pRotationAngleX = -2.363394;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 622 based Row/Col = 62/2
      geom_tower geom;
      geom.id = 622;
      geom.pDz = 7.507331;
      geom.pTheta = 0.044104;
      geom.pPhi = -2.972087;
      geom.pAlp1 = 0.043794;
      geom.pAlp2 = 0.043784;
      geom.pDy1 = 1.125494;
      geom.pDx1 = 1.230278;
      geom.pDx2 = 1.210480;
      geom.pDy2 = 1.239457;
      geom.pDx3 = 1.362063;
      geom.pDx4 = 1.340265;
      geom.centralX = -6.433568;
      geom.centralY = 105.705529;
      geom.centralZ = -105.590467;
      geom.pRotationAngleX = -2.363394;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 623 based Row/Col = 62/3
      geom_tower geom;
      geom.id = 623;
      geom.pDz = 7.507331;
      geom.pTheta = 0.027124;
      geom.pPhi = -2.863618;
      geom.pAlp1 = 0.026257;
      geom.pAlp2 = 0.026251;
      geom.pDy1 = 1.125494;
      geom.pDx1 = 1.228764;
      geom.pDx2 = 1.209035;
      geom.pDy2 = 1.239457;
      geom.pDx3 = 1.360387;
      geom.pDx4 = 1.338666;
      geom.centralX = -3.858584;
      geom.centralY = 105.705529;
      geom.centralZ = -105.590467;
      geom.pRotationAngleX = -2.363394;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 624 based Row/Col = 62/4
      geom_tower geom;
      geom.id = 624;
      geom.pDz = 7.507331;
      geom.pTheta = 0.011446;
      geom.pPhi = -2.433475;
      geom.pAlp1 = 0.008749;
      geom.pAlp2 = 0.008747;
      geom.pDy1 = 1.125494;
      geom.pDx1 = 1.228008;
      geom.pDx2 = 1.208314;
      geom.pDy2 = 1.239457;
      geom.pDx3 = 1.359551;
      geom.pDx4 = 1.337867;
      geom.centralX = -1.285935;
      geom.centralY = 105.705529;
      geom.centralZ = -105.590467;
      geom.pRotationAngleX = -2.363394;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 625 based Row/Col = 62/5
      geom_tower geom;
      geom.id = 625;
      geom.pDz = 7.507331;
      geom.pTheta = 0.011446;
      geom.pPhi = -0.708117;
      geom.pAlp1 = -0.008749;
      geom.pAlp2 = -0.008747;
      geom.pDy1 = 1.125494;
      geom.pDx1 = 1.228008;
      geom.pDx2 = 1.208314;
      geom.pDy2 = 1.239457;
      geom.pDx3 = 1.359551;
      geom.pDx4 = 1.337867;
      geom.centralX = 1.285935;
      geom.centralY = 105.705529;
      geom.centralZ = -105.590467;
      geom.pRotationAngleX = -2.363394;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 626 based Row/Col = 62/6
      geom_tower geom;
      geom.id = 626;
      geom.pDz = 7.507331;
      geom.pTheta = 0.027124;
      geom.pPhi = -0.277974;
      geom.pAlp1 = -0.026257;
      geom.pAlp2 = -0.026251;
      geom.pDy1 = 1.125494;
      geom.pDx1 = 1.228764;
      geom.pDx2 = 1.209035;
      geom.pDy2 = 1.239457;
      geom.pDx3 = 1.360387;
      geom.pDx4 = 1.338666;
      geom.centralX = 3.858584;
      geom.centralY = 105.705529;
      geom.centralZ = -105.590467;
      geom.pRotationAngleX = -2.363394;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 627 based Row/Col = 62/7
      geom_tower geom;
      geom.id = 627;
      geom.pDz = 7.507331;
      geom.pTheta = 0.044104;
      geom.pPhi = -0.169506;
      geom.pAlp1 = -0.043794;
      geom.pAlp2 = -0.043784;
      geom.pDy1 = 1.125494;
      geom.pDx1 = 1.230278;
      geom.pDx2 = 1.210480;
      geom.pDy2 = 1.239457;
      geom.pDx3 = 1.362063;
      geom.pDx4 = 1.340265;
      geom.centralX = 6.433568;
      geom.centralY = 105.705529;
      geom.centralZ = -105.590467;
      geom.pRotationAngleX = -2.363394;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 628 based Row/Col = 62/8
      geom_tower geom;
      geom.id = 628;
      geom.pDz = 7.507331;
      geom.pTheta = 0.061313;
      geom.pPhi = -0.121572;
      geom.pAlp1 = -0.061381;
      geom.pAlp2 = -0.061368;
      geom.pDy1 = 1.125494;
      geom.pDx1 = 1.232552;
      geom.pDx2 = 1.212651;
      geom.pDy2 = 1.239457;
      geom.pDx3 = 1.364581;
      geom.pDx4 = 1.342668;
      geom.centralX = 9.012453;
      geom.centralY = 105.705529;
      geom.centralZ = -105.590467;
      geom.pRotationAngleX = -2.363394;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1441 based Row/Col = 144/1
      geom_tower geom;
      geom.id = 1441;
      geom.pDz = 7.875766;
      geom.pTheta = 0.056555;
      geom.pPhi = -3.011788;
      geom.pAlp1 = -0.064967;
      geom.pAlp2 = -0.064953;
      geom.pDy1 = 1.121675;
      geom.pDx1 = 1.192532;
      geom.pDx2 = 1.213506;
      geom.pDy2 = 1.234325;
      geom.pDx3 = 1.318066;
      geom.pDx4 = 1.341142;
      geom.centralX = -8.863449;
      geom.centralY = 103.971312;
      geom.centralZ = 121.139966;
      geom.pRotationAngleX = -0.716975;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1442 based Row/Col = 144/2
      geom_tower geom;
      geom.id = 1442;
      geom.pDz = 7.875766;
      geom.pTheta = 0.040721;
      geom.pPhi = -2.960742;
      geom.pAlp1 = -0.046368;
      geom.pAlp2 = -0.046358;
      geom.pDy1 = 1.121675;
      geom.pDx1 = 1.190721;
      geom.pDx2 = 1.211602;
      geom.pDy2 = 1.234325;
      geom.pDx3 = 1.316065;
      geom.pDx4 = 1.339038;
      geom.centralX = -6.327780;
      geom.centralY = 103.971312;
      geom.centralZ = 121.139966;
      geom.pRotationAngleX = -0.716975;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1443 based Row/Col = 144/3
      geom_tower geom;
      geom.id = 1443;
      geom.pDz = 7.875766;
      geom.pTheta = 0.025126;
      geom.pPhi = -2.845691;
      geom.pAlp1 = -0.027806;
      geom.pAlp2 = -0.027799;
      geom.pDy1 = 1.121675;
      geom.pDx1 = 1.189516;
      geom.pDx2 = 1.210335;
      geom.pDy2 = 1.234325;
      geom.pDx3 = 1.314734;
      geom.pDx4 = 1.337638;
      geom.centralX = -3.795368;
      geom.centralY = 103.971312;
      geom.centralZ = 121.139966;
      geom.pRotationAngleX = -0.716975;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1444 based Row/Col = 144/4
      geom_tower geom;
      geom.id = 1444;
      geom.pDz = 7.875766;
      geom.pTheta = 0.010857;
      geom.pPhi = -2.400710;
      geom.pAlp1 = -0.009266;
      geom.pAlp2 = -0.009264;
      geom.pDy1 = 1.121675;
      geom.pDx1 = 1.188914;
      geom.pDx2 = 1.209702;
      geom.pDy2 = 1.234325;
      geom.pDx3 = 1.314069;
      geom.pDx4 = 1.336939;
      geom.centralX = -1.264906;
      geom.centralY = 103.971312;
      geom.centralZ = 121.139966;
      geom.pRotationAngleX = -0.716975;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1445 based Row/Col = 144/5
      geom_tower geom;
      geom.id = 1445;
      geom.pDz = 7.875766;
      geom.pTheta = 0.010857;
      geom.pPhi = -0.740883;
      geom.pAlp1 = 0.009266;
      geom.pAlp2 = 0.009264;
      geom.pDy1 = 1.121675;
      geom.pDx1 = 1.188914;
      geom.pDx2 = 1.209702;
      geom.pDy2 = 1.234325;
      geom.pDx3 = 1.314069;
      geom.pDx4 = 1.336939;
      geom.centralX = 1.264906;
      geom.centralY = 103.971312;
      geom.centralZ = 121.139966;
      geom.pRotationAngleX = -0.716975;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1446 based Row/Col = 144/6
      geom_tower geom;
      geom.id = 1446;
      geom.pDz = 7.875766;
      geom.pTheta = 0.025126;
      geom.pPhi = -0.295901;
      geom.pAlp1 = 0.027806;
      geom.pAlp2 = 0.027799;
      geom.pDy1 = 1.121675;
      geom.pDx1 = 1.189516;
      geom.pDx2 = 1.210335;
      geom.pDy2 = 1.234325;
      geom.pDx3 = 1.314734;
      geom.pDx4 = 1.337638;
      geom.centralX = 3.795368;
      geom.centralY = 103.971312;
      geom.centralZ = 121.139966;
      geom.pRotationAngleX = -0.716975;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1447 based Row/Col = 144/7
      geom_tower geom;
      geom.id = 1447;
      geom.pDz = 7.875766;
      geom.pTheta = 0.040721;
      geom.pPhi = -0.180850;
      geom.pAlp1 = 0.046368;
      geom.pAlp2 = 0.046358;
      geom.pDy1 = 1.121675;
      geom.pDx1 = 1.190721;
      geom.pDx2 = 1.211602;
      geom.pDy2 = 1.234325;
      geom.pDx3 = 1.316065;
      geom.pDx4 = 1.339038;
      geom.centralX = 6.327780;
      geom.centralY = 103.971312;
      geom.centralZ = 121.139966;
      geom.pRotationAngleX = -0.716975;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1448 based Row/Col = 144/8
      geom_tower geom;
      geom.id = 1448;
      geom.pDz = 7.875766;
      geom.pTheta = 0.056555;
      geom.pPhi = -0.129805;
      geom.pAlp1 = 0.064967;
      geom.pAlp2 = 0.064953;
      geom.pDy1 = 1.121675;
      geom.pDx1 = 1.192532;
      geom.pDx2 = 1.213506;
      geom.pDy2 = 1.234325;
      geom.pDx3 = 1.318066;
      geom.pDx4 = 1.341142;
      geom.centralX = 8.863449;
      geom.centralY = 103.971312;
      geom.centralZ = 121.139966;
      geom.pRotationAngleX = -0.716975;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1431 based Row/Col = 143/1
      geom_tower geom;
      geom.id = 1431;
      geom.pDz = 7.875742;
      geom.pTheta = 0.057429;
      geom.pPhi = 3.019971;
      geom.pAlp1 = -0.064966;
      geom.pAlp2 = -0.064951;
      geom.pDy1 = 1.121675;
      geom.pDx1 = 1.213506;
      geom.pDx2 = 1.234485;
      geom.pDy2 = 1.234325;
      geom.pDx3 = 1.341142;
      geom.pDx4 = 1.364223;
      geom.centralX = -9.016708;
      geom.centralY = 105.751027;
      geom.centralZ = 119.588534;
      geom.pRotationAngleX = -0.716975;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1432 based Row/Col = 143/2
      geom_tower geom;
      geom.id = 1432;
      geom.pDz = 7.875742;
      geom.pTheta = 0.041310;
      geom.pPhi = 2.972031;
      geom.pAlp1 = -0.046365;
      geom.pAlp2 = -0.046354;
      geom.pDy1 = 1.121675;
      geom.pDx1 = 1.211602;
      geom.pDx2 = 1.232484;
      geom.pDy2 = 1.234325;
      geom.pDx3 = 1.339038;
      geom.pDx4 = 1.362012;
      geom.centralX = -6.437085;
      geom.centralY = 105.751027;
      geom.centralZ = 119.588534;
      geom.pRotationAngleX = -0.716975;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1433 based Row/Col = 143/3
      geom_tower geom;
      geom.id = 1433;
      geom.pDz = 7.875742;
      geom.pTheta = 0.025406;
      geom.pPhi = 2.863543;
      geom.pAlp1 = -0.027803;
      geom.pAlp2 = -0.027796;
      geom.pDy1 = 1.121675;
      geom.pDx1 = 1.210335;
      geom.pDx2 = 1.231152;
      geom.pDy2 = 1.234325;
      geom.pDx3 = 1.337638;
      geom.pDx4 = 1.360541;
      geom.centralX = -3.860884;
      geom.centralY = 105.751027;
      geom.centralZ = 119.588534;
      geom.pRotationAngleX = -0.716975;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1434 based Row/Col = 143/4
      geom_tower geom;
      geom.id = 1434;
      geom.pDz = 7.875742;
      geom.pTheta = 0.010722;
      geom.pPhi = 2.433347;
      geom.pAlp1 = -0.009265;
      geom.pAlp2 = -0.009263;
      geom.pDy1 = 1.121675;
      geom.pDx1 = 1.209702;
      geom.pDx2 = 1.230487;
      geom.pDy2 = 1.234325;
      geom.pDx3 = 1.336939;
      geom.pDx4 = 1.359806;
      geom.centralX = -1.286734;
      geom.centralY = 105.751027;
      geom.centralZ = 119.588534;
      geom.pRotationAngleX = -0.716975;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1435 based Row/Col = 143/5
      geom_tower geom;
      geom.id = 1435;
      geom.pDz = 7.875742;
      geom.pTheta = 0.010722;
      geom.pPhi = 0.708246;
      geom.pAlp1 = 0.009265;
      geom.pAlp2 = 0.009263;
      geom.pDy1 = 1.121675;
      geom.pDx1 = 1.209702;
      geom.pDx2 = 1.230487;
      geom.pDy2 = 1.234325;
      geom.pDx3 = 1.336939;
      geom.pDx4 = 1.359806;
      geom.centralX = 1.286734;
      geom.centralY = 105.751027;
      geom.centralZ = 119.588534;
      geom.pRotationAngleX = -0.716975;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1436 based Row/Col = 143/6
      geom_tower geom;
      geom.id = 1436;
      geom.pDz = 7.875742;
      geom.pTheta = 0.025406;
      geom.pPhi = 0.278050;
      geom.pAlp1 = 0.027803;
      geom.pAlp2 = 0.027796;
      geom.pDy1 = 1.121675;
      geom.pDx1 = 1.210335;
      geom.pDx2 = 1.231152;
      geom.pDy2 = 1.234325;
      geom.pDx3 = 1.337638;
      geom.pDx4 = 1.360541;
      geom.centralX = 3.860884;
      geom.centralY = 105.751027;
      geom.centralZ = 119.588534;
      geom.pRotationAngleX = -0.716975;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1437 based Row/Col = 143/7
      geom_tower geom;
      geom.id = 1437;
      geom.pDz = 7.875742;
      geom.pTheta = 0.041310;
      geom.pPhi = 0.169562;
      geom.pAlp1 = 0.046365;
      geom.pAlp2 = 0.046354;
      geom.pDy1 = 1.121675;
      geom.pDx1 = 1.211602;
      geom.pDx2 = 1.232484;
      geom.pDy2 = 1.234325;
      geom.pDx3 = 1.339038;
      geom.pDx4 = 1.362012;
      geom.centralX = 6.437085;
      geom.centralY = 105.751027;
      geom.centralZ = 119.588534;
      geom.pRotationAngleX = -0.716975;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1438 based Row/Col = 143/8
      geom_tower geom;
      geom.id = 1438;
      geom.pDz = 7.875742;
      geom.pTheta = 0.057429;
      geom.pPhi = 0.121622;
      geom.pAlp1 = 0.064966;
      geom.pAlp2 = 0.064951;
      geom.pDy1 = 1.121675;
      geom.pDx1 = 1.213506;
      geom.pDx2 = 1.234485;
      geom.pDy2 = 1.234325;
      geom.pDx3 = 1.341142;
      geom.pDx4 = 1.364223;
      geom.centralX = 9.016708;
      geom.centralY = 105.751027;
      geom.centralZ = 119.588534;
      geom.pRotationAngleX = -0.716975;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 571 based Row/Col = 57/1
      geom_tower geom;
      geom.id = 571;
      geom.pDz = 7.875766;
      geom.pTheta = 0.056555;
      geom.pPhi = 3.011788;
      geom.pAlp1 = 0.064967;
      geom.pAlp2 = 0.064953;
      geom.pDy1 = 1.121675;
      geom.pDx1 = 1.213506;
      geom.pDx2 = 1.192532;
      geom.pDy2 = 1.234325;
      geom.pDx3 = 1.341142;
      geom.pDx4 = 1.318066;
      geom.centralX = -8.863449;
      geom.centralY = 103.971312;
      geom.centralZ = -121.139966;
      geom.pRotationAngleX = -2.424617;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 572 based Row/Col = 57/2
      geom_tower geom;
      geom.id = 572;
      geom.pDz = 7.875766;
      geom.pTheta = 0.040721;
      geom.pPhi = 2.960742;
      geom.pAlp1 = 0.046368;
      geom.pAlp2 = 0.046358;
      geom.pDy1 = 1.121675;
      geom.pDx1 = 1.211602;
      geom.pDx2 = 1.190721;
      geom.pDy2 = 1.234325;
      geom.pDx3 = 1.339038;
      geom.pDx4 = 1.316065;
      geom.centralX = -6.327780;
      geom.centralY = 103.971312;
      geom.centralZ = -121.139966;
      geom.pRotationAngleX = -2.424617;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 573 based Row/Col = 57/3
      geom_tower geom;
      geom.id = 573;
      geom.pDz = 7.875766;
      geom.pTheta = 0.025126;
      geom.pPhi = 2.845691;
      geom.pAlp1 = 0.027806;
      geom.pAlp2 = 0.027799;
      geom.pDy1 = 1.121675;
      geom.pDx1 = 1.210335;
      geom.pDx2 = 1.189516;
      geom.pDy2 = 1.234325;
      geom.pDx3 = 1.337638;
      geom.pDx4 = 1.314734;
      geom.centralX = -3.795368;
      geom.centralY = 103.971312;
      geom.centralZ = -121.139966;
      geom.pRotationAngleX = -2.424617;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 574 based Row/Col = 57/4
      geom_tower geom;
      geom.id = 574;
      geom.pDz = 7.875766;
      geom.pTheta = 0.010857;
      geom.pPhi = 2.400710;
      geom.pAlp1 = 0.009266;
      geom.pAlp2 = 0.009264;
      geom.pDy1 = 1.121675;
      geom.pDx1 = 1.209702;
      geom.pDx2 = 1.188914;
      geom.pDy2 = 1.234325;
      geom.pDx3 = 1.336939;
      geom.pDx4 = 1.314069;
      geom.centralX = -1.264906;
      geom.centralY = 103.971312;
      geom.centralZ = -121.139966;
      geom.pRotationAngleX = -2.424617;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 575 based Row/Col = 57/5
      geom_tower geom;
      geom.id = 575;
      geom.pDz = 7.875766;
      geom.pTheta = 0.010857;
      geom.pPhi = 0.740883;
      geom.pAlp1 = -0.009266;
      geom.pAlp2 = -0.009264;
      geom.pDy1 = 1.121675;
      geom.pDx1 = 1.209702;
      geom.pDx2 = 1.188914;
      geom.pDy2 = 1.234325;
      geom.pDx3 = 1.336939;
      geom.pDx4 = 1.314069;
      geom.centralX = 1.264906;
      geom.centralY = 103.971312;
      geom.centralZ = -121.139966;
      geom.pRotationAngleX = -2.424617;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 576 based Row/Col = 57/6
      geom_tower geom;
      geom.id = 576;
      geom.pDz = 7.875766;
      geom.pTheta = 0.025126;
      geom.pPhi = 0.295901;
      geom.pAlp1 = -0.027806;
      geom.pAlp2 = -0.027799;
      geom.pDy1 = 1.121675;
      geom.pDx1 = 1.210335;
      geom.pDx2 = 1.189516;
      geom.pDy2 = 1.234325;
      geom.pDx3 = 1.337638;
      geom.pDx4 = 1.314734;
      geom.centralX = 3.795368;
      geom.centralY = 103.971312;
      geom.centralZ = -121.139966;
      geom.pRotationAngleX = -2.424617;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 577 based Row/Col = 57/7
      geom_tower geom;
      geom.id = 577;
      geom.pDz = 7.875766;
      geom.pTheta = 0.040721;
      geom.pPhi = 0.180850;
      geom.pAlp1 = -0.046368;
      geom.pAlp2 = -0.046358;
      geom.pDy1 = 1.121675;
      geom.pDx1 = 1.211602;
      geom.pDx2 = 1.190721;
      geom.pDy2 = 1.234325;
      geom.pDx3 = 1.339038;
      geom.pDx4 = 1.316065;
      geom.centralX = 6.327780;
      geom.centralY = 103.971312;
      geom.centralZ = -121.139966;
      geom.pRotationAngleX = -2.424617;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 578 based Row/Col = 57/8
      geom_tower geom;
      geom.id = 578;
      geom.pDz = 7.875766;
      geom.pTheta = 0.056555;
      geom.pPhi = 0.129805;
      geom.pAlp1 = -0.064967;
      geom.pAlp2 = -0.064953;
      geom.pDy1 = 1.121675;
      geom.pDx1 = 1.213506;
      geom.pDx2 = 1.192532;
      geom.pDy2 = 1.234325;
      geom.pDx3 = 1.341142;
      geom.pDx4 = 1.318066;
      geom.centralX = 8.863449;
      geom.centralY = 103.971312;
      geom.centralZ = -121.139966;
      geom.pRotationAngleX = -2.424617;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 581 based Row/Col = 58/1
      geom_tower geom;
      geom.id = 581;
      geom.pDz = 7.875742;
      geom.pTheta = 0.057429;
      geom.pPhi = -3.019971;
      geom.pAlp1 = 0.064966;
      geom.pAlp2 = 0.064951;
      geom.pDy1 = 1.121675;
      geom.pDx1 = 1.234485;
      geom.pDx2 = 1.213506;
      geom.pDy2 = 1.234325;
      geom.pDx3 = 1.364223;
      geom.pDx4 = 1.341142;
      geom.centralX = -9.016708;
      geom.centralY = 105.751027;
      geom.centralZ = -119.588534;
      geom.pRotationAngleX = -2.424617;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 582 based Row/Col = 58/2
      geom_tower geom;
      geom.id = 582;
      geom.pDz = 7.875742;
      geom.pTheta = 0.041310;
      geom.pPhi = -2.972031;
      geom.pAlp1 = 0.046365;
      geom.pAlp2 = 0.046354;
      geom.pDy1 = 1.121675;
      geom.pDx1 = 1.232484;
      geom.pDx2 = 1.211602;
      geom.pDy2 = 1.234325;
      geom.pDx3 = 1.362012;
      geom.pDx4 = 1.339038;
      geom.centralX = -6.437085;
      geom.centralY = 105.751027;
      geom.centralZ = -119.588534;
      geom.pRotationAngleX = -2.424617;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 583 based Row/Col = 58/3
      geom_tower geom;
      geom.id = 583;
      geom.pDz = 7.875742;
      geom.pTheta = 0.025406;
      geom.pPhi = -2.863543;
      geom.pAlp1 = 0.027803;
      geom.pAlp2 = 0.027796;
      geom.pDy1 = 1.121675;
      geom.pDx1 = 1.231152;
      geom.pDx2 = 1.210335;
      geom.pDy2 = 1.234325;
      geom.pDx3 = 1.360541;
      geom.pDx4 = 1.337638;
      geom.centralX = -3.860884;
      geom.centralY = 105.751027;
      geom.centralZ = -119.588534;
      geom.pRotationAngleX = -2.424617;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 584 based Row/Col = 58/4
      geom_tower geom;
      geom.id = 584;
      geom.pDz = 7.875742;
      geom.pTheta = 0.010722;
      geom.pPhi = -2.433347;
      geom.pAlp1 = 0.009265;
      geom.pAlp2 = 0.009263;
      geom.pDy1 = 1.121675;
      geom.pDx1 = 1.230487;
      geom.pDx2 = 1.209702;
      geom.pDy2 = 1.234325;
      geom.pDx3 = 1.359806;
      geom.pDx4 = 1.336939;
      geom.centralX = -1.286734;
      geom.centralY = 105.751027;
      geom.centralZ = -119.588534;
      geom.pRotationAngleX = -2.424617;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 585 based Row/Col = 58/5
      geom_tower geom;
      geom.id = 585;
      geom.pDz = 7.875742;
      geom.pTheta = 0.010722;
      geom.pPhi = -0.708246;
      geom.pAlp1 = -0.009265;
      geom.pAlp2 = -0.009263;
      geom.pDy1 = 1.121675;
      geom.pDx1 = 1.230487;
      geom.pDx2 = 1.209702;
      geom.pDy2 = 1.234325;
      geom.pDx3 = 1.359806;
      geom.pDx4 = 1.336939;
      geom.centralX = 1.286734;
      geom.centralY = 105.751027;
      geom.centralZ = -119.588534;
      geom.pRotationAngleX = -2.424617;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 586 based Row/Col = 58/6
      geom_tower geom;
      geom.id = 586;
      geom.pDz = 7.875742;
      geom.pTheta = 0.025406;
      geom.pPhi = -0.278050;
      geom.pAlp1 = -0.027803;
      geom.pAlp2 = -0.027796;
      geom.pDy1 = 1.121675;
      geom.pDx1 = 1.231152;
      geom.pDx2 = 1.210335;
      geom.pDy2 = 1.234325;
      geom.pDx3 = 1.360541;
      geom.pDx4 = 1.337638;
      geom.centralX = 3.860884;
      geom.centralY = 105.751027;
      geom.centralZ = -119.588534;
      geom.pRotationAngleX = -2.424617;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 587 based Row/Col = 58/7
      geom_tower geom;
      geom.id = 587;
      geom.pDz = 7.875742;
      geom.pTheta = 0.041310;
      geom.pPhi = -0.169562;
      geom.pAlp1 = -0.046365;
      geom.pAlp2 = -0.046354;
      geom.pDy1 = 1.121675;
      geom.pDx1 = 1.232484;
      geom.pDx2 = 1.211602;
      geom.pDy2 = 1.234325;
      geom.pDx3 = 1.362012;
      geom.pDx4 = 1.339038;
      geom.centralX = 6.437085;
      geom.centralY = 105.751027;
      geom.centralZ = -119.588534;
      geom.pRotationAngleX = -2.424617;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 588 based Row/Col = 58/8
      geom_tower geom;
      geom.id = 588;
      geom.pDz = 7.875742;
      geom.pTheta = 0.057429;
      geom.pPhi = -0.121622;
      geom.pAlp1 = -0.064966;
      geom.pAlp2 = -0.064951;
      geom.pDy1 = 1.121675;
      geom.pDx1 = 1.234485;
      geom.pDx2 = 1.213506;
      geom.pDy2 = 1.234325;
      geom.pDx3 = 1.364223;
      geom.pDx4 = 1.341142;
      geom.centralX = 9.016708;
      geom.centralY = 105.751027;
      geom.centralZ = -119.588534;
      geom.pRotationAngleX = -2.424617;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1481 based Row/Col = 148/1
      geom_tower geom;
      geom.id = 1481;
      geom.pDz = 8.326725;
      geom.pTheta = 0.052693;
      geom.pPhi = -3.007283;
      geom.pAlp1 = -0.068144;
      geom.pAlp2 = -0.068129;
      geom.pDy1 = 1.116997;
      geom.pDx1 = 1.192154;
      geom.pDx2 = 1.214052;
      geom.pDy2 = 1.231781;
      geom.pDx3 = 1.315601;
      geom.pDx4 = 1.339743;
      geom.centralX = -8.859413;
      geom.centralY = 103.920957;
      geom.centralZ = 135.945213;
      geom.pRotationAngleX = -0.658852;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1482 based Row/Col = 148/2
      geom_tower geom;
      geom.id = 1482;
      geom.pDz = 8.326725;
      geom.pTheta = 0.037961;
      geom.pPhi = -2.954548;
      geom.pAlp1 = -0.048647;
      geom.pAlp2 = -0.048636;
      geom.pDy1 = 1.116997;
      geom.pDx1 = 1.190587;
      geom.pDx2 = 1.212398;
      geom.pDy2 = 1.231781;
      geom.pDx3 = 1.313872;
      geom.pDx4 = 1.337918;
      geom.centralX = -6.325332;
      geom.centralY = 103.920957;
      geom.centralZ = 135.945213;
      geom.pRotationAngleX = -0.658852;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1483 based Row/Col = 148/3
      geom_tower geom;
      geom.id = 1483;
      geom.pDz = 8.326725;
      geom.pTheta = 0.023467;
      geom.pPhi = -2.835958;
      geom.pAlp1 = -0.029177;
      geom.pAlp2 = -0.029170;
      geom.pDy1 = 1.116997;
      geom.pDx1 = 1.189544;
      geom.pDx2 = 1.211297;
      geom.pDy2 = 1.231781;
      geom.pDx3 = 1.312721;
      geom.pDx4 = 1.336703;
      geom.centralX = -3.794072;
      geom.centralY = 103.920957;
      geom.centralZ = 135.945213;
      geom.pRotationAngleX = -0.658852;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1484 based Row/Col = 148/4
      geom_tower geom;
      geom.id = 1484;
      geom.pDz = 8.326725;
      geom.pTheta = 0.010273;
      geom.pPhi = -2.383565;
      geom.pAlp1 = -0.009724;
      geom.pAlp2 = -0.009722;
      geom.pDy1 = 1.116997;
      geom.pDx1 = 1.189023;
      geom.pDx2 = 1.210746;
      geom.pDy2 = 1.231781;
      geom.pDx3 = 1.312146;
      geom.pDx4 = 1.336096;
      geom.centralX = -1.264503;
      geom.centralY = 103.920957;
      geom.centralZ = 135.945213;
      geom.pRotationAngleX = -0.658852;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1485 based Row/Col = 148/5
      geom_tower geom;
      geom.id = 1485;
      geom.pDz = 8.326725;
      geom.pTheta = 0.010273;
      geom.pPhi = -0.758027;
      geom.pAlp1 = 0.009724;
      geom.pAlp2 = 0.009722;
      geom.pDy1 = 1.116997;
      geom.pDx1 = 1.189023;
      geom.pDx2 = 1.210746;
      geom.pDy2 = 1.231781;
      geom.pDx3 = 1.312146;
      geom.pDx4 = 1.336096;
      geom.centralX = 1.264503;
      geom.centralY = 103.920957;
      geom.centralZ = 135.945213;
      geom.pRotationAngleX = -0.658852;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1486 based Row/Col = 148/6
      geom_tower geom;
      geom.id = 1486;
      geom.pDz = 8.326725;
      geom.pTheta = 0.023467;
      geom.pPhi = -0.305635;
      geom.pAlp1 = 0.029177;
      geom.pAlp2 = 0.029170;
      geom.pDy1 = 1.116997;
      geom.pDx1 = 1.189544;
      geom.pDx2 = 1.211297;
      geom.pDy2 = 1.231781;
      geom.pDx3 = 1.312721;
      geom.pDx4 = 1.336703;
      geom.centralX = 3.794072;
      geom.centralY = 103.920957;
      geom.centralZ = 135.945213;
      geom.pRotationAngleX = -0.658852;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1487 based Row/Col = 148/7
      geom_tower geom;
      geom.id = 1487;
      geom.pDz = 8.326725;
      geom.pTheta = 0.037961;
      geom.pPhi = -0.187044;
      geom.pAlp1 = 0.048647;
      geom.pAlp2 = 0.048636;
      geom.pDy1 = 1.116997;
      geom.pDx1 = 1.190587;
      geom.pDx2 = 1.212398;
      geom.pDy2 = 1.231781;
      geom.pDx3 = 1.313872;
      geom.pDx4 = 1.337918;
      geom.centralX = 6.325332;
      geom.centralY = 103.920957;
      geom.centralZ = 135.945213;
      geom.pRotationAngleX = -0.658852;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1488 based Row/Col = 148/8
      geom_tower geom;
      geom.id = 1488;
      geom.pDz = 8.326725;
      geom.pTheta = 0.052693;
      geom.pPhi = -0.134309;
      geom.pAlp1 = 0.068144;
      geom.pAlp2 = 0.068129;
      geom.pDy1 = 1.116997;
      geom.pDx1 = 1.192154;
      geom.pDx2 = 1.214052;
      geom.pDy2 = 1.231781;
      geom.pDx3 = 1.315601;
      geom.pDx4 = 1.339743;
      geom.centralX = 8.859413;
      geom.centralY = 103.920957;
      geom.centralZ = 135.945213;
      geom.pRotationAngleX = -0.658852;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1471 based Row/Col = 147/1
      geom_tower geom;
      geom.id = 1471;
      geom.pDz = 8.326697;
      geom.pTheta = 0.053578;
      geom.pPhi = 3.015910;
      geom.pAlp1 = -0.068143;
      geom.pAlp2 = -0.068127;
      geom.pDy1 = 1.116997;
      geom.pDx1 = 1.214052;
      geom.pDx2 = 1.235953;
      geom.pDy2 = 1.231781;
      geom.pDx3 = 1.339743;
      geom.pDx4 = 1.363889;
      geom.centralX = -9.019696;
      geom.centralY = 105.782078;
      geom.centralZ = 134.504208;
      geom.pRotationAngleX = -0.658852;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1472 based Row/Col = 147/2
      geom_tower geom;
      geom.id = 1472;
      geom.pDz = 8.326697;
      geom.pTheta = 0.038558;
      geom.pPhi = 2.966436;
      geom.pAlp1 = -0.048643;
      geom.pAlp2 = -0.048632;
      geom.pDy1 = 1.116997;
      geom.pDx1 = 1.212398;
      geom.pDx2 = 1.234208;
      geom.pDy2 = 1.231781;
      geom.pDx3 = 1.337918;
      geom.pDx4 = 1.361965;
      geom.centralX = -6.439665;
      geom.centralY = 105.782078;
      geom.centralZ = 134.504208;
      geom.pRotationAngleX = -0.658852;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1473 based Row/Col = 147/3
      geom_tower geom;
      geom.id = 1473;
      geom.pDz = 8.326697;
      geom.pTheta = 0.023752;
      geom.pPhi = 2.854692;
      geom.pAlp1 = -0.029174;
      geom.pAlp2 = -0.029167;
      geom.pDy1 = 1.116997;
      geom.pDx1 = 1.211297;
      geom.pDx2 = 1.233047;
      geom.pDy2 = 1.231781;
      geom.pDx3 = 1.336703;
      geom.pDx4 = 1.360684;
      geom.centralX = -3.862610;
      geom.centralY = 105.782078;
      geom.centralZ = 134.504208;
      geom.pRotationAngleX = -0.658852;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1474 based Row/Col = 147/4
      geom_tower geom;
      geom.id = 1474;
      geom.pDz = 8.326697;
      geom.pTheta = 0.010142;
      geom.pPhi = 2.416982;
      geom.pAlp1 = -0.009723;
      geom.pAlp2 = -0.009720;
      geom.pDy1 = 1.116997;
      geom.pDx1 = 1.210746;
      geom.pDx2 = 1.232467;
      geom.pDy2 = 1.231781;
      geom.pDx3 = 1.336096;
      geom.pDx4 = 1.360044;
      geom.centralX = -1.287339;
      geom.centralY = 105.782078;
      geom.centralZ = 134.504208;
      geom.pRotationAngleX = -0.658852;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1475 based Row/Col = 147/5
      geom_tower geom;
      geom.id = 1475;
      geom.pDz = 8.326697;
      geom.pTheta = 0.010142;
      geom.pPhi = 0.724610;
      geom.pAlp1 = 0.009723;
      geom.pAlp2 = 0.009720;
      geom.pDy1 = 1.116997;
      geom.pDx1 = 1.210746;
      geom.pDx2 = 1.232467;
      geom.pDy2 = 1.231781;
      geom.pDx3 = 1.336096;
      geom.pDx4 = 1.360044;
      geom.centralX = 1.287339;
      geom.centralY = 105.782078;
      geom.centralZ = 134.504208;
      geom.pRotationAngleX = -0.658852;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1476 based Row/Col = 147/6
      geom_tower geom;
      geom.id = 1476;
      geom.pDz = 8.326697;
      geom.pTheta = 0.023752;
      geom.pPhi = 0.286901;
      geom.pAlp1 = 0.029174;
      geom.pAlp2 = 0.029167;
      geom.pDy1 = 1.116997;
      geom.pDx1 = 1.211297;
      geom.pDx2 = 1.233047;
      geom.pDy2 = 1.231781;
      geom.pDx3 = 1.336703;
      geom.pDx4 = 1.360684;
      geom.centralX = 3.862610;
      geom.centralY = 105.782078;
      geom.centralZ = 134.504208;
      geom.pRotationAngleX = -0.658852;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1477 based Row/Col = 147/7
      geom_tower geom;
      geom.id = 1477;
      geom.pDz = 8.326697;
      geom.pTheta = 0.038558;
      geom.pPhi = 0.175156;
      geom.pAlp1 = 0.048643;
      geom.pAlp2 = 0.048632;
      geom.pDy1 = 1.116997;
      geom.pDx1 = 1.212398;
      geom.pDx2 = 1.234208;
      geom.pDy2 = 1.231781;
      geom.pDx3 = 1.337918;
      geom.pDx4 = 1.361965;
      geom.centralX = 6.439665;
      geom.centralY = 105.782078;
      geom.centralZ = 134.504208;
      geom.pRotationAngleX = -0.658852;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 1478 based Row/Col = 147/8
      geom_tower geom;
      geom.id = 1478;
      geom.pDz = 8.326697;
      geom.pTheta = 0.053578;
      geom.pPhi = 0.125683;
      geom.pAlp1 = 0.068143;
      geom.pAlp2 = 0.068127;
      geom.pDy1 = 1.116997;
      geom.pDx1 = 1.214052;
      geom.pDx2 = 1.235953;
      geom.pDy2 = 1.231781;
      geom.pDx3 = 1.339743;
      geom.pDx4 = 1.363889;
      geom.centralX = 9.019696;
      geom.centralY = 105.782078;
      geom.centralZ = 134.504208;
      geom.pRotationAngleX = -0.658852;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 531 based Row/Col = 53/1
      geom_tower geom;
      geom.id = 531;
      geom.pDz = 8.326725;
      geom.pTheta = 0.052693;
      geom.pPhi = 3.007283;
      geom.pAlp1 = 0.068144;
      geom.pAlp2 = 0.068129;
      geom.pDy1 = 1.116997;
      geom.pDx1 = 1.214052;
      geom.pDx2 = 1.192154;
      geom.pDy2 = 1.231781;
      geom.pDx3 = 1.339743;
      geom.pDx4 = 1.315601;
      geom.centralX = -8.859413;
      geom.centralY = 103.920957;
      geom.centralZ = -135.945213;
      geom.pRotationAngleX = -2.482741;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 532 based Row/Col = 53/2
      geom_tower geom;
      geom.id = 532;
      geom.pDz = 8.326725;
      geom.pTheta = 0.037961;
      geom.pPhi = 2.954548;
      geom.pAlp1 = 0.048647;
      geom.pAlp2 = 0.048636;
      geom.pDy1 = 1.116997;
      geom.pDx1 = 1.212398;
      geom.pDx2 = 1.190587;
      geom.pDy2 = 1.231781;
      geom.pDx3 = 1.337918;
      geom.pDx4 = 1.313872;
      geom.centralX = -6.325332;
      geom.centralY = 103.920957;
      geom.centralZ = -135.945213;
      geom.pRotationAngleX = -2.482741;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 533 based Row/Col = 53/3
      geom_tower geom;
      geom.id = 533;
      geom.pDz = 8.326725;
      geom.pTheta = 0.023467;
      geom.pPhi = 2.835958;
      geom.pAlp1 = 0.029177;
      geom.pAlp2 = 0.029170;
      geom.pDy1 = 1.116997;
      geom.pDx1 = 1.211297;
      geom.pDx2 = 1.189544;
      geom.pDy2 = 1.231781;
      geom.pDx3 = 1.336703;
      geom.pDx4 = 1.312721;
      geom.centralX = -3.794072;
      geom.centralY = 103.920957;
      geom.centralZ = -135.945213;
      geom.pRotationAngleX = -2.482741;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 534 based Row/Col = 53/4
      geom_tower geom;
      geom.id = 534;
      geom.pDz = 8.326725;
      geom.pTheta = 0.010273;
      geom.pPhi = 2.383565;
      geom.pAlp1 = 0.009724;
      geom.pAlp2 = 0.009722;
      geom.pDy1 = 1.116997;
      geom.pDx1 = 1.210746;
      geom.pDx2 = 1.189023;
      geom.pDy2 = 1.231781;
      geom.pDx3 = 1.336096;
      geom.pDx4 = 1.312146;
      geom.centralX = -1.264503;
      geom.centralY = 103.920957;
      geom.centralZ = -135.945213;
      geom.pRotationAngleX = -2.482741;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 535 based Row/Col = 53/5
      geom_tower geom;
      geom.id = 535;
      geom.pDz = 8.326725;
      geom.pTheta = 0.010273;
      geom.pPhi = 0.758027;
      geom.pAlp1 = -0.009724;
      geom.pAlp2 = -0.009722;
      geom.pDy1 = 1.116997;
      geom.pDx1 = 1.210746;
      geom.pDx2 = 1.189023;
      geom.pDy2 = 1.231781;
      geom.pDx3 = 1.336096;
      geom.pDx4 = 1.312146;
      geom.centralX = 1.264503;
      geom.centralY = 103.920957;
      geom.centralZ = -135.945213;
      geom.pRotationAngleX = -2.482741;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 536 based Row/Col = 53/6
      geom_tower geom;
      geom.id = 536;
      geom.pDz = 8.326725;
      geom.pTheta = 0.023467;
      geom.pPhi = 0.305635;
      geom.pAlp1 = -0.029177;
      geom.pAlp2 = -0.029170;
      geom.pDy1 = 1.116997;
      geom.pDx1 = 1.211297;
      geom.pDx2 = 1.189544;
      geom.pDy2 = 1.231781;
      geom.pDx3 = 1.336703;
      geom.pDx4 = 1.312721;
      geom.centralX = 3.794072;
      geom.centralY = 103.920957;
      geom.centralZ = -135.945213;
      geom.pRotationAngleX = -2.482741;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 537 based Row/Col = 53/7
      geom_tower geom;
      geom.id = 537;
      geom.pDz = 8.326725;
      geom.pTheta = 0.037961;
      geom.pPhi = 0.187044;
      geom.pAlp1 = -0.048647;
      geom.pAlp2 = -0.048636;
      geom.pDy1 = 1.116997;
      geom.pDx1 = 1.212398;
      geom.pDx2 = 1.190587;
      geom.pDy2 = 1.231781;
      geom.pDx3 = 1.337918;
      geom.pDx4 = 1.313872;
      geom.centralX = 6.325332;
      geom.centralY = 103.920957;
      geom.centralZ = -135.945213;
      geom.pRotationAngleX = -2.482741;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 538 based Row/Col = 53/8
      geom_tower geom;
      geom.id = 538;
      geom.pDz = 8.326725;
      geom.pTheta = 0.052693;
      geom.pPhi = 0.134309;
      geom.pAlp1 = -0.068144;
      geom.pAlp2 = -0.068129;
      geom.pDy1 = 1.116997;
      geom.pDx1 = 1.214052;
      geom.pDx2 = 1.192154;
      geom.pDy2 = 1.231781;
      geom.pDx3 = 1.339743;
      geom.pDx4 = 1.315601;
      geom.centralX = 8.859413;
      geom.centralY = 103.920957;
      geom.centralZ = -135.945213;
      geom.pRotationAngleX = -2.482741;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 541 based Row/Col = 54/1
      geom_tower geom;
      geom.id = 541;
      geom.pDz = 8.326697;
      geom.pTheta = 0.053578;
      geom.pPhi = -3.015910;
      geom.pAlp1 = 0.068143;
      geom.pAlp2 = 0.068127;
      geom.pDy1 = 1.116997;
      geom.pDx1 = 1.235953;
      geom.pDx2 = 1.214052;
      geom.pDy2 = 1.231781;
      geom.pDx3 = 1.363889;
      geom.pDx4 = 1.339743;
      geom.centralX = -9.019696;
      geom.centralY = 105.782078;
      geom.centralZ = -134.504208;
      geom.pRotationAngleX = -2.482741;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 542 based Row/Col = 54/2
      geom_tower geom;
      geom.id = 542;
      geom.pDz = 8.326697;
      geom.pTheta = 0.038558;
      geom.pPhi = -2.966436;
      geom.pAlp1 = 0.048643;
      geom.pAlp2 = 0.048632;
      geom.pDy1 = 1.116997;
      geom.pDx1 = 1.234208;
      geom.pDx2 = 1.212398;
      geom.pDy2 = 1.231781;
      geom.pDx3 = 1.361965;
      geom.pDx4 = 1.337918;
      geom.centralX = -6.439665;
      geom.centralY = 105.782078;
      geom.centralZ = -134.504208;
      geom.pRotationAngleX = -2.482741;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 543 based Row/Col = 54/3
      geom_tower geom;
      geom.id = 543;
      geom.pDz = 8.326697;
      geom.pTheta = 0.023752;
      geom.pPhi = -2.854692;
      geom.pAlp1 = 0.029174;
      geom.pAlp2 = 0.029167;
      geom.pDy1 = 1.116997;
      geom.pDx1 = 1.233047;
      geom.pDx2 = 1.211297;
      geom.pDy2 = 1.231781;
      geom.pDx3 = 1.360684;
      geom.pDx4 = 1.336703;
      geom.centralX = -3.862610;
      geom.centralY = 105.782078;
      geom.centralZ = -134.504208;
      geom.pRotationAngleX = -2.482741;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 544 based Row/Col = 54/4
      geom_tower geom;
      geom.id = 544;
      geom.pDz = 8.326697;
      geom.pTheta = 0.010142;
      geom.pPhi = -2.416982;
      geom.pAlp1 = 0.009723;
      geom.pAlp2 = 0.009720;
      geom.pDy1 = 1.116997;
      geom.pDx1 = 1.232467;
      geom.pDx2 = 1.210746;
      geom.pDy2 = 1.231781;
      geom.pDx3 = 1.360044;
      geom.pDx4 = 1.336096;
      geom.centralX = -1.287339;
      geom.centralY = 105.782078;
      geom.centralZ = -134.504208;
      geom.pRotationAngleX = -2.482741;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 545 based Row/Col = 54/5
      geom_tower geom;
      geom.id = 545;
      geom.pDz = 8.326697;
      geom.pTheta = 0.010142;
      geom.pPhi = -0.724610;
      geom.pAlp1 = -0.009723;
      geom.pAlp2 = -0.009720;
      geom.pDy1 = 1.116997;
      geom.pDx1 = 1.232467;
      geom.pDx2 = 1.210746;
      geom.pDy2 = 1.231781;
      geom.pDx3 = 1.360044;
      geom.pDx4 = 1.336096;
      geom.centralX = 1.287339;
      geom.centralY = 105.782078;
      geom.centralZ = -134.504208;
      geom.pRotationAngleX = -2.482741;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 546 based Row/Col = 54/6
      geom_tower geom;
      geom.id = 546;
      geom.pDz = 8.326697;
      geom.pTheta = 0.023752;
      geom.pPhi = -0.286901;
      geom.pAlp1 = -0.029174;
      geom.pAlp2 = -0.029167;
      geom.pDy1 = 1.116997;
      geom.pDx1 = 1.233047;
      geom.pDx2 = 1.211297;
      geom.pDy2 = 1.231781;
      geom.pDx3 = 1.360684;
      geom.pDx4 = 1.336703;
      geom.centralX = 3.862610;
      geom.centralY = 105.782078;
      geom.centralZ = -134.504208;
      geom.pRotationAngleX = -2.482741;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 547 based Row/Col = 54/7
      geom_tower geom;
      geom.id = 547;
      geom.pDz = 8.326697;
      geom.pTheta = 0.038558;
      geom.pPhi = -0.175156;
      geom.pAlp1 = -0.048643;
      geom.pAlp2 = -0.048632;
      geom.pDy1 = 1.116997;
      geom.pDx1 = 1.234208;
      geom.pDx2 = 1.212398;
      geom.pDy2 = 1.231781;
      geom.pDx3 = 1.361965;
      geom.pDx4 = 1.337918;
      geom.centralX = 6.439665;
      geom.centralY = 105.782078;
      geom.centralZ = -134.504208;
      geom.pRotationAngleX = -2.482741;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }
    {
      // tower 548 based Row/Col = 54/8
      geom_tower geom;
      geom.id = 548;
      geom.pDz = 8.326697;
      geom.pTheta = 0.053578;
      geom.pPhi = -0.125683;
      geom.pAlp1 = -0.068143;
      geom.pAlp2 = -0.068127;
      geom.pDy1 = 1.116997;
      geom.pDx1 = 1.235953;
      geom.pDx2 = 1.214052;
      geom.pDy2 = 1.231781;
      geom.pDx3 = 1.363889;
      geom.pDx4 = 1.339743;
      geom.centralX = 9.019696;
      geom.centralY = 105.782078;
      geom.centralZ = -134.504208;
      geom.pRotationAngleX = -2.482741;
      geom.ModuleSkinThickness = 0.010000;
      geom.NFiberX = 30;
      geom.NFiberY = 48;
      sector_tower_map[geom.id] = geom;
    }

}
