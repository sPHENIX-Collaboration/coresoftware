#include "INTTDeadMapv1.h"

#include <iostream>

using namespace std;

const INTTDeadMapv1::Map&
INTTDeadMapv1::getDeadChannels(void) const
{
  return m_DeadChannels;
}

INTTDeadMapv1::Map&
INTTDeadMapv1::getDeadChannels(void)
{
  return m_DeadChannels;
}

void INTTDeadMapv1::addDeadChannel(PHG4CellDefs::keytype key)
{
  m_DeadChannels.insert(key);
}

bool INTTDeadMapv1::isDeadChannel(PHG4CellDefs::keytype key) const
{
  auto it = m_DeadChannels.find(key);
  if (it != m_DeadChannels.end())
  {
    return true;
  }
  return false;
}

int INTTDeadMapv1::isValid() const
{
  return size() > 0;
}

void INTTDeadMapv1::Reset()
{
  m_DeadChannels.clear();
}

void INTTDeadMapv1::identify(std::ostream& os) const
{
  os << "INTTDeadMapv1, number of dead channel & sensors: " << size() << std::endl;
}
