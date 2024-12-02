#include "TriggerPrimitive.h"
#include "TriggerDefs.h"

#include <cmath>
#include <iostream>

namespace
{
  TriggerPrimitive::Map dummy_map;
}

//______________________________________
void TriggerPrimitive::identify(std::ostream& out) const
{
  out << "identify yourself: I am a TriggerPrimitive object" << std::endl;
}

//__________________________________________________________
TriggerPrimitive::ConstRange TriggerPrimitive::getSums() const
{
  return std::make_pair(dummy_map.begin(), dummy_map.end());
}
TriggerPrimitive::Range TriggerPrimitive::getSums()
{
  return std::make_pair(dummy_map.begin(), dummy_map.end());
}
