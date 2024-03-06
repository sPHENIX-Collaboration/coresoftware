
#include "LL1ReturnCodes.h"
#include "LL1Outv2.h"
#include "TriggerDefs.h"
#include "TriggerPrimitive.h"
#include "TriggerPrimitiveContainerv1.h"

#include <cmath>
#include <iostream>

ClassImp(LL1Outv2)

LL1Outv2::LL1Outv2()
{
  _ll1_type = "NONE";
  _trigger_type = "NONE";

  setTriggerKey(TriggerDefs::getTriggerKey(TriggerDefs::GetTriggerId(_trigger_type)));

  Init();
}


LL1Outv2::LL1Outv2(std::string triggertype, std::string ll1type)
{

  setTriggerType(triggertype);
  setTriggerKey(TriggerDefs::getTriggerKey(TriggerDefs::GetTriggerId(triggertype)));

  setLL1Type(ll1type);

  Init();
}

LL1Outv2::~LL1Outv2()
{
  Reset();
  delete _trigger_primitives;
  delete _trigger_bits;
}

//______________________________________
void LL1Outv2::Init()
{
  _trigger_primitives = new TriggerPrimitiveContainerv1();
  _trigger_primitives->setTriggerType(TriggerDefs::GetTriggerId(_trigger_type));
  _trigger_bits = new std::vector<unsigned int>();
}


//______________________________________
void LL1Outv2::Reset()
{
  _trigger_bits->clear();
  _trigger_primitives->Reset();

  while (_trigger_words.begin() != _trigger_words.end())
    {
      delete _trigger_words.begin()->second;
      _trigger_words.erase(_trigger_words.begin());
    } 
  

}

LL1Outv2::ConstRange LL1Outv2::getTriggerWords() const
{
  return make_pair(_trigger_words.begin(), _trigger_words.end());
}  

LL1Outv2::Range LL1Outv2::getTriggerWords()
{
  return make_pair(_trigger_words.begin(), _trigger_words.end());
}

//______________________________________
void LL1Outv2::identify(std::ostream& out) const
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
  out << __FILE__ << __FUNCTION__ <<" Trigger Primitive Container: "<< _trigger_primitives->size() << std::endl;
  if (_trigger_primitives->size()) _trigger_primitives->identify();

}

int LL1Outv2::isValid() const
{
  return _trigger_primitives->isValid();
}

