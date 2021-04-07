#include "RawTower.h"

RawTower::CellMap DummyCellMap;
RawTower::ShowerMap DummyShowerMap;

RawTower::CellIterator RawTower::find_g4cell(CellKeyType id)
{
  return DummyCellMap.end();
}

RawTower::CellConstIterator RawTower::find_g4cell(CellKeyType id) const
{
  return DummyCellMap.end();
}

RawTower::CellConstRange RawTower::get_g4cells() const
{
    PHOOL_VIRTUAL_WARN("get_g4cells()");
    return CellConstRange(DummyCellMap.begin(),DummyCellMap.end());
}

RawTower::ShowerConstRange RawTower::get_g4showers() const
{
    PHOOL_VIRTUAL_WARN("get_g4showers()");
    return ShowerConstRange(DummyShowerMap.begin(), DummyShowerMap.end());
}

RawTower::ShowerIterator RawTower::find_g4shower(int id)
{
  return DummyShowerMap.end();
}

RawTower::ShowerConstIterator RawTower::find_g4shower(int id) const
{
  return DummyShowerMap.end();
}
