
#include "TriggerPrimitiveContainer.h"

#include <cmath>
#include <iostream>

namespace
{
  TriggerPrimitiveContainer::Map dummy_map;
}

//______________________________________
void TriggerPrimitiveContainer::identify(std::ostream& out) const
{
  out << "identify yourself: I am a TriggerPrimitiveContainer object" << std::endl;
}

TriggerPrimitiveContainer::ConstRange TriggerPrimitiveContainer::getTriggerPrimitives() const
{
  return std::make_pair(dummy_map.begin(), dummy_map.end());
}

TriggerPrimitiveContainer::Range TriggerPrimitiveContainer::getTriggerPrimitives()
{
  return std::make_pair(dummy_map.begin(), dummy_map.end());
}
