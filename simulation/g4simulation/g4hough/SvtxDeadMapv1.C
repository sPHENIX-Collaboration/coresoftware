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

void SvtxDeadMapv1::addDeadChannel(PHG4CellDefs::keytype key)
{
  m_DeadChannels.insert(key);
}

bool SvtxDeadMapv1::isDeadChannel(PHG4CellDefs::keytype key) const
{
  auto it = m_DeadChannels.find(key);
  if (it != m_DeadChannels.end())
  {
    return true;
  }
  return false;
}

int SvtxDeadMapv1::isValid() const
{
  return size() > 0;
}

void SvtxDeadMapv1::Reset()
{
  m_DeadChannels.clear();
}

void SvtxDeadMapv1::identify(std::ostream& os) const
{
  os << "SvtxDeadMapv1, number of towers: " << size() << std::endl;
}
