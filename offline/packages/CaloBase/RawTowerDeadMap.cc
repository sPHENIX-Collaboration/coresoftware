#include "RawTowerDeadMap.h"

#include <iostream>

const RawTowerDeadMap::Map&
RawTowerDeadMap::getDeadTowers(void) const
{
  static Map tmp_map;
  return tmp_map;
}

RawTowerDeadMap::Map&
RawTowerDeadMap::getDeadTowers(void)
{
  static Map tmp_map;
  return tmp_map;
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
