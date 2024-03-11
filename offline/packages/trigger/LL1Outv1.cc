
#include "LL1ReturnCodes.h"
#include "LL1Outv1.h"
#include "TriggerDefs.h"
#include "TriggerPrimitivev1.h"
#include "TriggerPrimitiveContainerv1.h"

#include <cmath>
#include <iostream>

ClassImp(LL1Outv1)

LL1Outv1::LL1Outv1()
: _ll1_type("NONE")
  ,_trigger_type("NONE")
{
  _trigger_key = TriggerDefs::getTriggerKey(TriggerDefs::GetTriggerId(_trigger_type));

  _trigger_bits = new std::vector<unsigned int>();

}


LL1Outv1::LL1Outv1(const std::string triggertype, const std::string ll1type)
{
  _trigger_type = triggertype;

  _trigger_key = TriggerDefs::getTriggerKey(TriggerDefs::GetTriggerId(triggertype));

  _ll1_type = ll1type;

  _trigger_bits = new std::vector<unsigned int>();
}

LL1Outv1::~LL1Outv1()
{
  Reset();
}

//______________________________________
void LL1Outv1::Reset()
{
  _trigger_bits->clear();

  while (_trigger_words.begin() != _trigger_words.end())
    {
      delete _trigger_words.begin()->second;
      _trigger_words.erase(_trigger_words.begin());
    } 
  

}

LL1Outv1::ConstRange LL1Outv1::getTriggerWords() const
{
  return make_pair(_trigger_words.begin(), _trigger_words.end());
}  

LL1Outv1::Range LL1Outv1::getTriggerWords()
{
  return make_pair(_trigger_words.begin(), _trigger_words.end());
}

//______________________________________
void LL1Outv1::identify(std::ostream& out) const
{
  out << __FILE__ << __FUNCTION__ <<" LL1Out: Triggertype = " << _trigger_type << " LL1 type = " << _ll1_type << std::endl;
  out << __FILE__ << __FUNCTION__ <<" Event number: "<< _event_number <<"    Clock: "<< _clock_number << std::endl;
  out << __FILE__ << __FUNCTION__ <<" Trigger bits    "<<_trigger_bits->size()<<"         : ";
  for (int i = 0; i < static_cast<int>(_trigger_bits->size()) ; i++) cout <<" "<< _trigger_bits->at(i);
  out << " " <<std::endl;
  out << __FILE__ << __FUNCTION__ <<" Trigger words    "<<std::endl;
  
  for (auto i = _trigger_words.begin(); i != _trigger_words.end() ; ++i) 
    {
      for (auto j = (*i).second->begin(); j != (*i).second->end(); ++j)cout <<" "<< (*j);
      
      out << " " <<std::endl;
    }
}

int LL1Outv1::isValid() const
{
  return 0;
}

