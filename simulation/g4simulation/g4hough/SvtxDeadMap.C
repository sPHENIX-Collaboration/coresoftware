#include "SvtxDeadMap.h"
#include <iostream>

using namespace std;

const SvtxDeadMap::Map&
SvtxDeadMap::getDeadChannels(void) const
{
  static Map tmp_map;
  return tmp_map;
}

SvtxDeadMap::Map&
SvtxDeadMap::getDeadChannels(void)
{
  static Map tmp_map;
  return tmp_map;
}

void SvtxDeadMap::addDeadChannel(const unsigned int layer, const unsigned int ieta, const int unsigned iphi)
{
}

void SvtxDeadMap::addDeadChannel(PHG4CellDefs::keytype key)
{
}

bool
SvtxDeadMap::isDeadChannel(PHG4CellDefs::keytype key)
{
  return false;
}

bool
SvtxDeadMap::isDeadChannel(const unsigned int layer, const unsigned int ieta, const unsigned int iphi)
{
  return false;
}

int SvtxDeadMap::isValid() const
{
  return size()>0;
}

void SvtxDeadMap::Reset()
{
}

void SvtxDeadMap::identify(std::ostream& os) const
{
  os << "SvtxDeadMap" << std::endl;
}

