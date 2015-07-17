// $$Id: PHG4CylinderGeom_Spacalv3.h,v 1.3 2014/08/28 22:18:35 jinhuang Exp $$

/*!
 * \file ${file_name}
 * \brief
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $$Revision: 1.3 $$
 * \date $$Date: 2014/08/28 22:18:35 $$
 */
#ifndef PHG4CylinderGeom_Spacalv3_H__
#define PHG4CylinderGeom_Spacalv3_H__

#include "PHG4CylinderGeom_Spacalv2.h"
#include <string>
#include <map>

class PHG4CylinderGeom_Spacalv3 : public PHG4CylinderGeom_Spacalv2
{
public:
  PHG4CylinderGeom_Spacalv3();

  virtual
  ~PHG4CylinderGeom_Spacalv3()
  {
  }

  virtual void
  identify(std::ostream& os = std::cout) const;
  virtual void
  Print(Option_t* option = "") const;
  virtual void
  SetDefault();


  double
  get_sidewall_outer_torr() const
  {
    return sidewall_outer_torr;
  }

  void
  set_sidewall_outer_torr(double sidewallOuterTorr)
  {
    sidewall_outer_torr = sidewallOuterTorr;
  }

  double
  get_sidewall_thickness() const
  {
    return sidewall_thickness;
  }

  void
  set_sidewall_thickness(double sidewallThickness)
  {
    sidewall_thickness = sidewallThickness;
  }

  class geom_super_tower
  {
  public:
    int id;
    double pDz;
    double pDy1;
    double pDx1;
    double pDx2;
    double pDy2;
    double pDx3;
    double pDx4;
    double pRotationAngleX;
    double centralY;
    double centralZ;
    geom_super_tower();

  ClassDef(PHG4CylinderGeom_Spacalv3::geom_super_tower,1)

  };
  typedef std::map<int, geom_super_tower> geom_super_tower_map_t;

  void
  load_demo_geom_super_tower_map();

//  geom_super_tower_map_t
//  get_geom_super_tower_map() const
//  {
//    return geom_super_tower_map;
//  }
//
//  void
//  set_geom_super_tower_map(geom_super_tower_map_t geomSuperTowerMap)
//  {
//    geom_super_tower_map = geomSuperTowerMap;
//  }
protected:
  double sidewall_thickness;
  double sidewall_outer_torr;
  geom_super_tower_map_t geom_super_tower_map;

ClassDef(PHG4CylinderGeom_Spacalv3,2)

};

#endif
