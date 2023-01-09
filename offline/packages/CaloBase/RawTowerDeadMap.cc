#include "RawTowerDeadMap.h"

#include <iostream>

const RawTowerDeadMap::Map&
RawTowerDeadMap::getDeadTowers() const
{
  static Map tmp_map;
  return tmp_map;
}

RawTowerDeadMap::Map&
RawTowerDeadMap::getDeadTowers()
{
  static Map tmp_map;
  return tmp_map;
}

void RawTowerDeadMap::addDeadTower(const unsigned int /*ieta*/, const int unsigned /*iphi*/)
{
}

void RawTowerDeadMap::addDeadTower(RawTowerDefs::keytype /*key*/)
{
}

bool RawTowerDeadMap::isDeadTower(RawTowerDefs::keytype /*key*/)
{
  return false;
}

bool RawTowerDeadMap::isDeadTower(const unsigned int /*ieta*/, const unsigned int /*iphi*/)
{
  return false;
}

int RawTowerDeadMap::isValid() const
{
  return size() > 0;
}

void RawTowerDeadMap::Reset()
{
}

void RawTowerDeadMap::identify(std::ostream& os) const
{
  os << "RawTowerDeadMap" << std::endl;
}
