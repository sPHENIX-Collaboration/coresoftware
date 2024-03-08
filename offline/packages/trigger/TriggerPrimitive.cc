
#include "TriggerPrimitive.h"
#include "TriggerDefs.h"
#include <cmath>
#include <iostream>

ClassImp(TriggerPrimitive)

TriggerPrimitive::TriggerPrimitive()
{

}

TriggerPrimitive::TriggerPrimitive(TriggerDefs::TriggerPrimKey key)
{
  m_triggerprimkey = key;
}

TriggerPrimitive::~TriggerPrimitive()
{

}

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
