
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
  _sums.clear();
}

void TriggerPrimitive::add_sum(TriggerDefs::TriggerSumKey key, std::vector<unsigned int> *sum)
{
  _sums[key] = sum;
}
std::vector<unsigned int>*  TriggerPrimitive::get_sum_at_key(TriggerDefs::TriggerSumKey key)
{
  if (!_sums[key]) return nullptr;

  return _sums[key];
}

TriggerPrimitive::ConstRange TriggerPrimitive::getSums() const
{
  return make_pair(_sums.begin(), _sums.end());
}

TriggerPrimitive::Range TriggerPrimitive::getSums()
{
  return make_pair(_sums.begin(), _sums.end());
}


//______________________________________
void TriggerPrimitive::identify(std::ostream& out)
{
  out << __FILE__<<__FUNCTION__<<":: primitive key : "<< m_triggerprimkey << std::endl;
  Range range = getSums();
  for (auto iter = range.first; iter != range.second ; ++iter )
    {
      out << "Sumkey "<<(*iter).first <<" :";
      for (auto i = (*iter).second->begin(); i != (*iter).second->end(); i++) out << (*i) << " "; 
      out <<" "<<std::endl;
    }
}

int TriggerPrimitive::isValid()
{
  return (!_sums.empty());
}
