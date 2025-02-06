#include "LL1Out.h"

#include <iostream>

namespace
{
  LL1Out::Map dummy_map;
}

void LL1Out::identify(std::ostream& os) const
{
  os << "virtual LL1Out object" << std::endl;
}

LL1Out::ConstRange LL1Out::getTriggerWords() const
{
  return std::make_pair(dummy_map.begin(), dummy_map.end());
}
LL1Out::Range LL1Out::getTriggerWords()
{
  return std::make_pair(dummy_map.begin(), dummy_map.end());
}

std::vector<std::pair<TriggerDefs::TriggerSumKey, unsigned short>> LL1Out::getTriggeredSums()
{
  return std::vector<std::pair<TriggerDefs::TriggerSumKey, unsigned short>>();
}

std::vector<TriggerDefs::TriggerSumKey> LL1Out::getTriggeredSumKeys(int /*ith*/)
{
  return std::vector<TriggerDefs::TriggerSumKey>();
}

std::vector<TriggerDefs::TriggerPrimKey> LL1Out::getTriggeredPrimitives()
{
  return std::vector<TriggerDefs::TriggerPrimKey>();
}
