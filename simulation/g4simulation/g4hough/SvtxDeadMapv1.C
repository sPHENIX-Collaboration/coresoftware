#include "SvtxDeadMapv1.h"

#include <iostream>

using namespace std;

const SvtxDeadMapv1::Map&
SvtxDeadMapv1::getDeadChannels(void) const
{
  return m_DeadChannels;
}

SvtxDeadMapv1::Map&
SvtxDeadMapv1::getDeadChannels(void)
{
  return m_DeadChannels;
}

void SvtxDeadMapv1::addDeadChannel(const unsigned int layer, const unsigned int ieta, const int unsigned iphi)
{
  PHG4CellDefs::keytype key = PHG4CellDefs::EtaPhiBinning::genkey(layer, ieta, iphi);
  m_DeadChannels.insert(key);
}

void SvtxDeadMapv1::addDeadChannel(PHG4CellDefs::keytype key)
{
  m_DeadChannels.insert(key);
}

bool
SvtxDeadMapv1::isDeadChannel(PHG4CellDefs::keytype key)
{
  auto it = m_DeadChannels.find(key);
  if (it != m_DeadChannels.end())
  {
    return true;
  }
  return false;
}

bool
SvtxDeadMapv1::isDeadChannel(const unsigned int layer, const unsigned int ieta, const unsigned int iphi)
{
  PHG4CellDefs::keytype key = PHG4CellDefs::EtaPhiBinning::genkey(layer, ieta, iphi);
  return isDeadChannel(key);
}

int SvtxDeadMapv1::isValid() const
{
  return size()>0;
}

void SvtxDeadMapv1::Reset()
{
  m_DeadChannels.clear();
}

void SvtxDeadMapv1::identify(std::ostream& os) const
{
  os << "SvtxDeadMapv1, number of towers: " << size() << std::endl;
}

