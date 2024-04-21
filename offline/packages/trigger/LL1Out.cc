#include "LL1Out.h"

#include "LL1ReturnCodes.h"

#include <cmath>
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
