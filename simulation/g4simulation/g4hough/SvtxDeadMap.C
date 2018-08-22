#include "SvtxDeadMap.h"

#include <cassert>
#include <limits>
#include <iostream>

using namespace std;

PHG4CellDefs::keytype s_wildCardID = std::numeric_limits<PHG4CellDefs::keytype>::max();

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
  const PHG4CellDefs::keytype key = PHG4CellDefs::EtaPhiBinning::genkey(layer, ieta, iphi);
  addDeadChannel(key);
}

void SvtxDeadMap::addDeadChannel(PHG4CellDefs::keytype key)
{
}

void SvtxDeadMap::addDeadChannelINTT(const unsigned int layer,
                                     const unsigned int ladder_z, const unsigned int ladder_phi,
                                     const unsigned int strip_z, const unsigned int strip_phi)
{
  addDeadChannel(getINTTKey(layer,
                            ladder_z, ladder_phi,
                            strip_z, strip_phi));
}

bool SvtxDeadMap::isDeadChannel(PHG4CellDefs::keytype key)
{
  return false;
}

bool SvtxDeadMap::isDeadChannel(const unsigned int layer, const unsigned int ieta, const unsigned int iphi)
{
  const PHG4CellDefs::keytype key = PHG4CellDefs::EtaPhiBinning::genkey(layer, ieta, iphi);
  return isDeadChannel(key);
}
bool SvtxDeadMap::isDeadChannel(const unsigned int layer,
                                const unsigned int ladder_z, const unsigned int ladder_phi,
                                const unsigned int strip_z, const unsigned int strip_phi)
{
  return isDeadChannel(getINTTKey(layer,
                                  ladder_z, ladder_phi,
                                  strip_z, strip_phi));
}

int SvtxDeadMap::isValid() const
{
  return size() > 0;
}

void SvtxDeadMap::Reset()
{
}

void SvtxDeadMap::identify(std::ostream& os) const
{
  os << "SvtxDeadMap" << std::endl;
}

PHG4CellDefs::keytype SvtxDeadMap::getINTTKey(const unsigned int layer,
                                              const unsigned int ladder_z, const unsigned int ladder_phi,
                                              const unsigned int strip_z, const unsigned int strip_phi)
{
  static const unsigned int layer_bit = 8;
  static const unsigned int ladder_phi_bit = 16;
  static const unsigned int ladder_z_bit = 8;
  static const unsigned int strip_z_bit = 16;
  static const unsigned int strip_phi_bit = 16;

  assert(layer_bit + ladder_phi_bit + ladder_z_bit + strip_z_bit + strip_phi_bit == numeric_limits<PHG4CellDefs::keytype>::digits());
  assert(ladder_phi < (1 << ladder_phi_bit));
  assert(ladder_z < (1 << ladder_z_bit));
  assert(strip_z < (1 << strip_z_bit));
  assert(strip_phi < (1 << strip_phi_bit));

  PHG4CellDefs::keytype key = 0;

  key += layer;

  key <<= ladder_phi_bit;
  key += ladder_phi;

  key <<= ladder_z_bit;
  key += ladder_z;

  key <<= strip_z_bit;
  key += strip_z;

  key <<= strip_phi_bit;
  key += strip_phi;

  return key;
}
