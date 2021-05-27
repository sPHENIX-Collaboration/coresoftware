#ifndef CALOBASE_RAWTOWERGEOMCONTAINERV1_H
#define CALOBASE_RAWTOWERGEOMCONTAINERV1_H

#include "RawTowerGeomContainer.h"

#include "RawTowerDefs.h"

#include <iostream>

class RawTowerGeom;

/*! \class RawTowerGeomContainerv1
    \brief Generic tower geometry class, store each tower's geometry individually
*/
class RawTowerGeomContainerv1 : public RawTowerGeomContainer
{
 public:
  RawTowerGeomContainerv1(RawTowerDefs::CalorimeterId caloid = RawTowerDefs::NONE);
  ~RawTowerGeomContainerv1() override;

  void Reset() override;
  int isValid() const override;
  void identify(std::ostream &os = std::cout) const override;

  void set_calorimeter_id(RawTowerDefs::CalorimeterId caloid) override { _caloid = caloid; }
  RawTowerDefs::CalorimeterId get_calorimeter_id() const override { return _caloid; }

  ConstIterator add_tower_geometry(RawTowerGeom *geo) override;
  RawTowerGeom *get_tower_geometry(RawTowerDefs::keytype key) override;

  //! return all tower geometries
  ConstRange get_tower_geometries(void) const override;
  Range get_tower_geometries(void) override;

  unsigned int size() const override { return _geoms.size(); }

 protected:
  RawTowerDefs::CalorimeterId _caloid;
  Map _geoms;

  ClassDefOverride(RawTowerGeomContainerv1, 1)
};

#endif
