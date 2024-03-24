
#include "TriggerPrimitive.h"
#include "TriggerDefs.h"
#include <cmath>
#include <iostream>

ClassImp(TriggerPrimitive)

namespace
{
  TriggerPrimitive::Map dummy_map;
}

TriggerPrimitive::TriggerPrimitive()
= default;

TriggerPrimitive::TriggerPrimitive(TriggerDefs::TriggerPrimKey key)
{
  m_triggerprimkey = key;
}

TriggerPrimitive::~TriggerPrimitive()
= default;

//______________________________________
void TriggerPrimitive::Reset()
{

}

void TriggerPrimitive::add_sum(TriggerDefs::TriggerSumKey /*key*/, std::vector<unsigned int> */*sum*/)
{

}

//______________________________________
void TriggerPrimitive::identify(std::ostream& out) const
{
  out << "identify yourself: I am a TriggerPrimitive object" << std::endl;
}

int TriggerPrimitive::isValid() const
{
  return (!_sums.empty());
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
