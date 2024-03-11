
#include "TriggerPrimitiveContainer.h"

#include <cmath>
#include <iostream>

ClassImp(TriggerPrimitiveContainer)

namespace
{
  TriggerPrimitiveContainer::Map dummy_map;
}

TriggerPrimitiveContainer::TriggerPrimitiveContainer()
{
}
TriggerPrimitiveContainer::~TriggerPrimitiveContainer()
{

}

//______________________________________
void TriggerPrimitiveContainer::Reset()
{

}


//______________________________________
void TriggerPrimitiveContainer::identify(std::ostream& out) const
{
  out << "identify yourself: I am a TriggerPrimitiveContainer object" << std::endl;

}

int TriggerPrimitiveContainer::isValid() const
{
  return 1;
}
  

TriggerPrimitiveContainer::ConstRange TriggerPrimitiveContainer::getTriggerPrimitives() const
{
  return std::make_pair(dummy_map.begin(), dummy_map.end());
}

TriggerPrimitiveContainer::Range TriggerPrimitiveContainer::getTriggerPrimitives()
{
  return std::make_pair(dummy_map.begin(), dummy_map.end());
}
