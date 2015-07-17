// $$Id: PHG4CylinderGeom_Spacalv3.cc,v 1.3 2014/08/28 22:18:35 jinhuang Exp $$

/*!
 * \file ${file_name}
 * \brief
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $$Revision: 1.3 $$
 * \date $$Date: 2014/08/28 22:18:35 $$
 */

#include "PHG4CylinderGeom_Spacalv3.h"

#include <Geant4/globals.hh>
#include <Geant4/G4PhysicalConstants.hh>

#include <cmath>
#include <limits>       // std::numeric_limits
#include <iostream>

ClassImp(PHG4CylinderGeom_Spacalv3)
ClassImp(PHG4CylinderGeom_Spacalv3::geom_super_tower)

using namespace std;

PHG4CylinderGeom_Spacalv3::PHG4CylinderGeom_Spacalv3()
{
  SetDefault();
}

void
PHG4CylinderGeom_Spacalv3::identify(std::ostream& os) const
{
  os << "PHG4CylinderGeom_Spacalv3: layer: " << layer //
      << ", radius: " << radius //
      << ", thickness: " << thickness //
      << ", zmin: " << zmin //
      << ", zmax: " << zmax << //
      ", num scint: " << nscint

      << endl;
  return;
}

void
PHG4CylinderGeom_Spacalv3::Print(Option_t *opt) const
{
  PHG4CylinderGeom_Spacalv2::Print(opt);

}

void
PHG4CylinderGeom_Spacalv3::SetDefault()
{
  PHG4CylinderGeom_Spacalv2::SetDefault();

  radius = 90.000000;
  thickness = 26.130000;
  zmin = 149.470000;
  zmax = -zmin;
  azimuthal_n_sec = 32;
  polar_taper_ratio = 1.;
  assembly_spacing = 0.002500;
  sidewall_thickness = 0.075000;
  sidewall_outer_torr = 0.030000;

}

PHG4CylinderGeom_Spacalv3::geom_super_tower::geom_super_tower() :
    id(numeric_limits<int>::min()), pDz(
        numeric_limits<double>::signaling_NaN()), //
    pDy1(numeric_limits<double>::signaling_NaN()), //
    pDx1(numeric_limits<double>::signaling_NaN()), //
    pDx2(numeric_limits<double>::signaling_NaN()), //
    pDy2(numeric_limits<double>::signaling_NaN()), //
    pDx3(numeric_limits<double>::signaling_NaN()), //
    pDx4(numeric_limits<double>::signaling_NaN()), //
    pRotationAngleX(numeric_limits<double>::signaling_NaN()), //
    centralY(numeric_limits<double>::signaling_NaN()), //
    centralZ(numeric_limits<double>::signaling_NaN())
{
}

void
PHG4CylinderGeom_Spacalv3::load_demo_geom_super_tower_map()
{
  // Chris Cullen 2D spacal design July 2015
  radius = 90.000000;
  thickness = 26.130000;
  zmin = 149.470000;
  zmax = -zmin;
  azimuthal_n_sec = 32;
  azimuthal_tilt = 0;
  azimuthal_seg_visible = false;
  polar_taper_ratio = 1.;
  assembly_spacing = 0.002500;
  sidewall_thickness = 0.075000;
  sidewall_outer_torr = 0.030000;
  geom_super_tower_map.clear();
    {
      // Super tower 1 based on Chris Cullen tower ID 1 1
      geom_super_tower geom;
      geom.id = 1.000000;
      geom.pDz = 6.749451;
      geom.pDy1 = 2.244890;
      geom.pDx1 = 9.569762;
      geom.pDx2 = 9.580269;
      geom.pDy2 = 2.565402;
      geom.pDx3 = 10.898655;
      geom.pDx4 = 10.910663;
      geom.centralY = 105.088944;
      geom.centralZ = 2.483168;
      geom.pRotationAngleX = -1.547066;
      geom_super_tower_map[geom.id] = geom;
    }
    {
      // Super tower 3 based on Chris Cullen tower ID 3 1
      geom_super_tower geom;
      geom.id = 3.000000;
      geom.pDz = 6.749476;
      geom.pDy1 = 2.246020;
      geom.pDx1 = 9.535014;
      geom.pDx2 = 9.585823;
      geom.pDy2 = 2.564232;
      geom.pDx3 = 10.852649;
      geom.pDx4 = 10.910663;
      geom.centralY = 104.898063;
      geom.centralZ = 12.114086;
      geom.pRotationAngleX = -1.455797;
      geom_super_tower_map[geom.id] = geom;
    }
    {
      // Super tower 5 based on Chris Cullen tower ID 5 1
      geom_super_tower geom;
      geom.id = 5.000000;
      geom.pDz = 6.749471;
      geom.pDy1 = 2.247152;
      geom.pDx1 = 9.511848;
      geom.pDx2 = 9.602284;
      geom.pDy2 = 2.562931;
      geom.pDx3 = 10.807519;
      geom.pDx4 = 10.910663;
      geom.centralY = 104.766491;
      geom.centralZ = 21.838412;
      geom.pRotationAngleX = -1.365252;
      geom_super_tower_map[geom.id] = geom;
    }
    {
      // Super tower 7 based on Chris Cullen tower ID 7 1
      geom_super_tower geom;
      geom.id = 7.000000;
      geom.pDz = 6.749418;
      geom.pDy1 = 2.249126;
      geom.pDx1 = 9.500365;
      geom.pDx2 = 9.629027;
      geom.pDy2 = 2.557714;
      geom.pDx3 = 10.764340;
      geom.pDx4 = 10.910663;
      geom.centralY = 104.695625;
      geom.centralZ = 31.745365;
      geom.pRotationAngleX = -1.276437;
      geom_super_tower_map[geom.id] = geom;
    }
    {
      // Super tower 9 based on Chris Cullen tower ID 9 1
      geom_super_tower geom;
      geom.id = 9.000000;
      geom.pDz = 6.749378;
      geom.pDy1 = 2.250932;
      geom.pDx1 = 9.500015;
      geom.pDx2 = 9.665026;
      geom.pDy2 = 2.550134;
      geom.pDx3 = 10.723738;
      geom.pDx4 = 10.910663;
      geom.centralY = 104.683052;
      geom.centralZ = 41.925557;
      geom.pRotationAngleX = -1.189890;
      geom_super_tower_map[geom.id] = geom;
    }
    {
      // Super tower 11 based on Chris Cullen tower ID 11 1
      geom_super_tower geom;
      geom.id = 11.000000;
      geom.pDz = 6.749340;
      geom.pDy1 = 2.254011;
      geom.pDx1 = 9.509672;
      geom.pDx2 = 9.708705;
      geom.pDy2 = 2.543841;
      geom.pDx3 = 10.686087;
      geom.pDx4 = 10.910663;
      geom.centralY = 104.722867;
      geom.centralZ = 52.466304;
      geom.pRotationAngleX = -1.106480;
      geom_super_tower_map[geom.id] = geom;
    }
    {
      // Super tower 13 based on Chris Cullen tower ID 13 1
      geom_super_tower geom;
      geom.id = 13.000000;
      geom.pDz = 6.749279;
      geom.pDy1 = 2.256703;
      geom.pDx1 = 9.528309;
      geom.pDx2 = 9.758764;
      geom.pDy2 = 2.534883;
      geom.pDx3 = 10.651839;
      geom.pDx4 = 10.910663;
      geom.centralY = 104.810306;
      geom.centralZ = 63.468504;
      geom.pRotationAngleX = -1.026390;
      geom_super_tower_map[geom.id] = geom;
    }
    {
      // Super tower 15 based on Chris Cullen tower ID 15 1
      geom_super_tower geom;
      geom.id = 15.000000;
      geom.pDz = 6.825456;
      geom.pDy1 = 2.257814;
      geom.pDx1 = 9.542269;
      geom.pDx2 = 9.801118;
      geom.pDy2 = 2.524663;
      geom.pDx3 = 10.621269;
      geom.pDx4 = 10.910663;
      geom.centralY = 104.875648;
      geom.centralZ = 74.978133;
      geom.pRotationAngleX = -0.950304;
      geom_super_tower_map[geom.id] = geom;
    }
    {
      // Super tower 17 based on Chris Cullen tower ID 17 1
      geom_super_tower geom;
      geom.id = 17.000000;
      geom.pDz = 7.047629;
      geom.pDy1 = 2.253643;
      geom.pDx1 = 9.541893;
      geom.pDx2 = 9.825634;
      geom.pDy2 = 2.511961;
      geom.pDx3 = 10.594451;
      geom.pDx4 = 10.910663;
      geom.centralY = 104.868853;
      geom.centralZ = 87.029925;
      geom.pRotationAngleX = -0.878224;
      geom_super_tower_map[geom.id] = geom;
    }
    {
      // Super tower 19 based on Chris Cullen tower ID 19 1
      geom_super_tower geom;
      geom.id = 19.000000;
      geom.pDz = 7.326995;
      geom.pDy1 = 2.248571;
      geom.pDx1 = 9.542018;
      geom.pDx2 = 9.847599;
      geom.pDy2 = 2.499107;
      geom.pDx3 = 10.571060;
      geom.pDx4 = 10.910663;
      geom.centralY = 104.865551;
      geom.centralZ = 99.766946;
      geom.pRotationAngleX = -0.810340;
      geom_super_tower_map[geom.id] = geom;
    }
    {
      // Super tower 21 based on Chris Cullen tower ID 21 1
      geom_super_tower geom;
      geom.id = 21.000000;
      geom.pDz = 7.682634;
      geom.pDy1 = 2.245016;
      geom.pDx1 = 9.539692;
      geom.pDx2 = 9.864485;
      geom.pDy2 = 2.488919;
      geom.pDx3 = 10.550621;
      geom.pDx4 = 10.910663;
      geom.centralY = 104.850629;
      geom.centralZ = 113.268570;
      geom.pRotationAngleX = -0.746991;
      geom_super_tower_map[geom.id] = geom;
    }
    {
      // Super tower 23 based on Chris Cullen tower ID 23 1
      geom_super_tower geom;
      geom.id = 23.000000;
      geom.pDz = 8.101892;
      geom.pDy1 = 2.235546;
      geom.pDx1 = 9.538416;
      geom.pDx2 = 9.879070;
      geom.pDy2 = 2.478731;
      geom.pDx3 = 10.532984;
      geom.pDx4 = 10.910663;
      geom.centralY = 104.839643;
      geom.centralZ = 127.643953;
      geom.pRotationAngleX = -0.687482;
      geom_super_tower_map[geom.id] = geom;
    }
    {
      // Super tower 2 based on Chris Cullen tower ID 2 1
      geom_super_tower geom;
      geom.id = 2.000000;
      geom.pDz = 6.749477;
      geom.pDy1 = 2.251979;
      geom.pDx1 = 9.550999;
      geom.pDx2 = 9.581820;
      geom.pDy2 = 2.548983;
      geom.pDx3 = 10.875790;
      geom.pDx4 = 10.910663;
      geom.centralY = 104.987217;
      geom.centralZ = 7.292911;
      geom.pRotationAngleX = -1.501344;
      geom_super_tower_map[geom.id] = geom;
    }
    {
      // Super tower 4 based on Chris Cullen tower ID 4 1
      geom_super_tower geom;
      geom.id = 4.000000;
      geom.pDz = 6.749468;
      geom.pDy1 = 2.253132;
      geom.pDx1 = 9.522130;
      geom.pDx2 = 9.593078;
      geom.pDy2 = 2.547734;
      geom.pDx3 = 10.830459;
      geom.pDx4 = 10.910663;
      geom.centralY = 104.827451;
      geom.centralZ = 16.959834;
      geom.pRotationAngleX = -1.410438;
      geom_super_tower_map[geom.id] = geom;
    }
    {
      // Super tower 6 based on Chris Cullen tower ID 6 1
      geom_super_tower geom;
      geom.id = 6.000000;
      geom.pDz = 6.749473;
      geom.pDy1 = 2.255103;
      geom.pDx1 = 9.504943;
      geom.pDx2 = 9.615042;
      geom.pDy2 = 2.542583;
      geom.pDx3 = 10.786530;
      geom.pDx4 = 10.910663;
      geom.centralY = 104.728073;
      geom.centralZ = 26.765695;
      geom.pRotationAngleX = -1.320578;
      geom_super_tower_map[geom.id] = geom;
    }
    {
      // Super tower 8 based on Chris Cullen tower ID 8 1
      geom_super_tower geom;
      geom.id = 8.000000;
      geom.pDz = 6.749421;
      geom.pDy1 = 2.257010;
      geom.pDx1 = 9.499165;
      geom.pDx2 = 9.646739;
      geom.pDy2 = 2.537410;
      geom.pDx3 = 10.744777;
      geom.pDx4 = 10.910663;
      geom.centralY = 104.687878;
      geom.centralZ = 36.799456;
      geom.pRotationAngleX = -1.232816;
      geom_super_tower_map[geom.id] = geom;
    }
    {
      // Super tower 10 based on Chris Cullen tower ID 10 1
      geom_super_tower geom;
      geom.id = 10.000000;
      geom.pDz = 6.749389;
      geom.pDy1 = 2.258607;
      geom.pDx1 = 9.504193;
      geom.pDx2 = 9.687016;
      geom.pDy2 = 2.527313;
      geom.pDx3 = 10.706151;
      geom.pDx4 = 10.910663;
      geom.centralY = 104.704833;
      geom.centralZ = 47.150147;
      geom.pRotationAngleX = -1.147836;
      geom_super_tower_map[geom.id] = geom;
    }
    {
      // Super tower 12 based on Chris Cullen tower ID 12 1
      geom_super_tower geom;
      geom.id = 12.000000;
      geom.pDz = 6.749261;
      geom.pDy1 = 2.262680;
      geom.pDx1 = 9.518553;
      geom.pDx2 = 9.734323;
      geom.pDy2 = 2.519675;
      geom.pDx3 = 10.670402;
      geom.pDx4 = 10.910663;
      geom.centralY = 104.770619;
      geom.centralZ = 57.912699;
      geom.pRotationAngleX = -1.065998;
      geom_super_tower_map[geom.id] = geom;
    }
    {
      // Super tower 14 based on Chris Cullen tower ID 14 1
      geom_super_tower geom;
      geom.id = 14.000000;
      geom.pDz = 6.749218;
      geom.pDy1 = 2.266766;
      geom.pDx1 = 9.541168;
      geom.pDx2 = 9.787208;
      geom.pDy2 = 2.509530;
      geom.pDx3 = 10.638305;
      geom.pDx4 = 10.910663;
      geom.centralY = 104.880791;
      geom.centralZ = 69.182107;
      geom.pRotationAngleX = -0.987818;
      geom_super_tower_map[geom.id] = geom;
    }
    {
      // Super tower 16 based on Chris Cullen tower ID 16 1
      geom_super_tower geom;
      geom.id = 16.000000;
      geom.pDz = 6.920622;
      geom.pDy1 = 2.262448;
      geom.pDx1 = 9.544195;
      geom.pDx2 = 9.816703;
      geom.pDy2 = 2.499247;
      geom.pDx3 = 10.609686;
      geom.pDx4 = 10.910663;
      geom.centralY = 104.890697;
      geom.centralZ = 80.944022;
      geom.pRotationAngleX = -0.913654;
      geom_super_tower_map[geom.id] = geom;
    }
    {
      // Super tower 18 based on Chris Cullen tower ID 18 1
      geom_super_tower geom;
      geom.id = 18.000000;
      geom.pDz = 7.174582;
      geom.pDy1 = 2.258718;
      geom.pDx1 = 9.543169;
      geom.pDx2 = 9.839218;
      geom.pDy2 = 2.491528;
      geom.pDx3 = 10.584144;
      geom.pDx4 = 10.910663;
      geom.centralY = 104.880410;
      geom.centralZ = 93.324172;
      geom.pRotationAngleX = -0.843676;
      geom_super_tower_map[geom.id] = geom;
    }
    {
      // Super tower 20 based on Chris Cullen tower ID 20 1
      geom_super_tower geom;
      geom.id = 20.000000;
      geom.pDz = 7.504837;
      geom.pDy1 = 2.253489;
      geom.pDx1 = 9.540217;
      geom.pDx2 = 9.856705;
      geom.pDy2 = 2.481413;
      geom.pDx3 = 10.562204;
      geom.pDx4 = 10.910663;
      geom.centralY = 104.861614;
      geom.centralZ = 106.422317;
      geom.pRotationAngleX = -0.778053;
      geom_super_tower_map[geom.id] = geom;
    }
    {
      // Super tower 22 based on Chris Cullen tower ID 22 1
      geom_super_tower geom;
      geom.id = 22.000000;
      geom.pDz = 7.873254;
      geom.pDy1 = 2.245851;
      geom.pDx1 = 9.540868;
      geom.pDx2 = 9.874717;
      geom.pDy2 = 2.471151;
      geom.pDx3 = 10.543366;
      geom.pDx4 = 10.910663;
      geom.centralY = 104.861169;
      geom.centralZ = 120.364250;
      geom.pRotationAngleX = -0.716799;
      geom_super_tower_map[geom.id] = geom;
    }
    {
      // Super tower 24 based on Chris Cullen tower ID 24 1
      geom_super_tower geom;
      geom.id = 24.000000;
      geom.pDz = 8.324211;
      geom.pDy1 = 2.236494;
      geom.pDx1 = 9.540117;
      geom.pDx2 = 9.888851;
      geom.pDy2 = 2.466062;
      geom.pDx3 = 10.526180;
      geom.pDx4 = 10.910663;
      geom.centralY = 104.851518;
      geom.centralZ = 135.224711;
      geom.pRotationAngleX = -0.658682;
      geom_super_tower_map[geom.id] = geom;
    }

}
