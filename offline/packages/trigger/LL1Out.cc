
#include "LL1ReturnCodes.h"
#include "LL1Out.h"

#include <cmath>
#include <iostream>

ClassImp(LL1Out)
namespace
{
  LL1Out::Map dummy_map;
}

LL1Out::LL1Out()
= default;

LL1Out::~LL1Out()
= default;

//______________________________________
void LL1Out::Reset()
{

}

void LL1Out::identify(std::ostream& os) const
{
  os << "virtual LL1Out object" << std::endl;
}

int LL1Out::isValid() const
{
  return 1;
}

LL1Out::ConstRange LL1Out::getTriggerWords() const
{
  return std::make_pair(dummy_map.begin(), dummy_map.end());
}
LL1Out::Range LL1Out::getTriggerWords()
{
  return std::make_pair(dummy_map.begin(), dummy_map.end());
}



