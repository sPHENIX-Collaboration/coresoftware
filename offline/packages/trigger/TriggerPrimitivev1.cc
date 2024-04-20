#include "TriggerPrimitivev1.h"
#include "TriggerDefs.h"

#include <cmath>
#include <iostream>

TriggerPrimitivev1::TriggerPrimitivev1(TriggerDefs::TriggerPrimKey key)
{
  m_triggerprimkey = key;
}

//______________________________________
void TriggerPrimitivev1::Reset()
{
  _sums.clear();
}

void TriggerPrimitivev1::add_sum(TriggerDefs::TriggerSumKey key, std::vector<unsigned int>* sum)
{
  _sums[key] = sum;
}
std::vector<unsigned int>* TriggerPrimitivev1::get_sum_at_key(TriggerDefs::TriggerSumKey key)
{
  if (!_sums[key])
  {
    return nullptr;
  }

  return _sums[key];
}

TriggerPrimitivev1::ConstRange TriggerPrimitivev1::getSums() const
{
  return make_pair(_sums.begin(), _sums.end());
}

TriggerPrimitivev1::Range TriggerPrimitivev1::getSums()
{
  return make_pair(_sums.begin(), _sums.end());
}

//______________________________________
void TriggerPrimitivev1::identify(std::ostream& out) const
{
  out << __FILE__ << __FUNCTION__ << ":: primitive key : " << std::hex << m_triggerprimkey << std::endl;
  out << " TriggerId: " << TriggerDefs::getTriggerId_from_TriggerPrimKey(m_triggerprimkey) << std::endl;
  out << " DetectorId: " << TriggerDefs::getDetectorId_from_TriggerPrimKey(m_triggerprimkey) << std::endl;
  out << " PrimitiveId: " << TriggerDefs::getPrimitiveId_from_TriggerPrimKey(m_triggerprimkey) << std::endl;

  ConstRange range = getSums();
  for (auto iter = range.first; iter != range.second; ++iter)
  {
    out << "Sumkey " << (*iter).first << " :";
    for (unsigned int& i : *(*iter).second)
    {
      out << i << " ";
    }
    out << " " << std::endl;
  }
}

int TriggerPrimitivev1::isValid() const
{
  return (!_sums.empty());
}
