#include "InttDeadMapv1.h"

#include <g4detectors/PHG4CellDefs.h>  // for keytype

#include <iostream>
#include <map>  // for _Rb_tree_const_iterator

const InttDeadMapv1::Map&
InttDeadMapv1::getDeadChannels(void) const
{
  return m_DeadChannels;
}

InttDeadMapv1::Map&
InttDeadMapv1::getDeadChannels(void)
{
  return m_DeadChannels;
}

void InttDeadMapv1::addDeadChannel(PHG4CellDefs::keytype key)
{
  m_DeadChannels.insert(key);
}

bool InttDeadMapv1::isDeadChannel(PHG4CellDefs::keytype key) const
{
  auto it = m_DeadChannels.find(key);
  if (it != m_DeadChannels.end())
  {
    return true;
  }
  return false;
}

int InttDeadMapv1::isValid() const
{
  return size() > 0;
}

void InttDeadMapv1::Reset()
{
  m_DeadChannels.clear();
}

void InttDeadMapv1::identify(std::ostream& os) const
{
  os << "InttDeadMapv1, number of dead channel & sensors: " << size() << std::endl;
}
