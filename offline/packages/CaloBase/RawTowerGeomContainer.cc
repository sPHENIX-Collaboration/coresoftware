#include "RawTowerGeomContainer.h"

#include <phool/phool.h>  // for PHOOL_VIRTUAL_WARN

#include <iostream>

RawTowerGeomContainer::Map DummyMap;

void RawTowerGeomContainer::identify(std::ostream& os) const
{
  os << "Base class RawTowerGeomContainer." << std::endl;
}

RawTowerGeomContainer::ConstIterator RawTowerGeomContainer::add_tower_geometry(RawTowerGeom* /*geo*/)
{
  PHOOL_VIRTUAL_WARN("add_tower_geometry()");
  return DummyMap.end();
}

RawTowerGeomContainer::ConstRange RawTowerGeomContainer::get_tower_geometries() const
{
  PHOOL_VIRTUAL_WARN("get_tower_geometries()");
  return ConstRange(DummyMap.begin(), DummyMap.end());
}

RawTowerGeomContainer::Range RawTowerGeomContainer::get_tower_geometries()
{
  PHOOL_VIRTUAL_WARN("get_tower_geometries()");
  return Range(DummyMap.begin(), DummyMap.end());
}
