#include "RawTower.h"

#include <phool/phool.h>  // for PHOOL_VIRTUAL_WARN

#include <cstdlib>  // for exit

RawTower::CellMap DummyCellMap;
RawTower::ShowerMap DummyShowerMap;

RawTower::CellIterator RawTower::find_g4cell(CellKeyType /*id*/)
{
  return DummyCellMap.end();
}

RawTower::CellConstIterator RawTower::find_g4cell(CellKeyType /*id*/) const
{
  return DummyCellMap.end();
}

RawTower::CellConstRange RawTower::get_g4cells() const
{
  PHOOL_VIRTUAL_WARN("get_g4cells()");
  return CellConstRange(DummyCellMap.begin(), DummyCellMap.end());
}

RawTower::ShowerConstRange RawTower::get_g4showers() const
{
  PHOOL_VIRTUAL_WARN("get_g4showers()");
  return ShowerConstRange(DummyShowerMap.begin(), DummyShowerMap.end());
}

RawTower::ShowerIterator RawTower::find_g4shower(int /*id*/)
{
  return DummyShowerMap.end();
}

RawTower::ShowerConstIterator RawTower::find_g4shower(int /*id*/) const
{
  return DummyShowerMap.end();
}

const std::string RawTower::get_property_info(RawTower::PROPERTY prop_id)
{
  switch (prop_id)
  {
  case prop_scint_gammas:
    return "Scintillation photon count or energy";
  case prop_cerenkov_gammas:
    return "Cherenkov photon count or energy";

  default:
    std::cout << "RawTower::get_property_info - Fatal Error - unknown index " << prop_id << std::endl;
    exit(1);
  }
}
