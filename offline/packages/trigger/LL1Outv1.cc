#include "LL1Outv1.h"
#include "LL1ReturnCodes.h"
#include "TriggerDefs.h"
#include "TriggerPrimitiveContainerv1.h"
#include "TriggerPrimitivev1.h"

#include <cmath>
#include <iostream>

LL1Outv1::LL1Outv1()
  : _trigger_key(TriggerDefs::getTriggerKey(TriggerDefs::GetTriggerId(_trigger_type)))
{
  _trigger_bits = new std::vector<unsigned int>();
}

LL1Outv1::LL1Outv1(const std::string& triggertype, const std::string& ll1type)
  : _trigger_key(TriggerDefs::getTriggerKey(TriggerDefs::GetTriggerId(triggertype)))
  , _ll1_type(ll1type)
  , _trigger_type(triggertype)
{
  _trigger_bits = new std::vector<unsigned int>();
}

LL1Outv1::~LL1Outv1()
{
  // you cannot call a virtual method in the ctor, you need to be specific
  // which one you want to call
  LL1Outv1::Reset();
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
  out << __FILE__ << __FUNCTION__ << " LL1Out: Triggertype = " << _trigger_type << " LL1 type = " << _ll1_type << std::endl;
  out << __FILE__ << __FUNCTION__ << " Event number: " << _event_number << "    Clock: " << _clock_number << std::endl;
  out << __FILE__ << __FUNCTION__ << " Trigger bits    " << _trigger_bits->size() << " : ";
  for (unsigned int _trigger_bit : *_trigger_bits)
  {
    std::cout << " " << _trigger_bit;
  }
  out << " " << std::endl;
  out << __FILE__ << __FUNCTION__ << " Trigger words    " << std::endl;

  for (auto _trigger_word : _trigger_words)
  {
    for (unsigned int& j : *_trigger_word.second)
    {
      std::cout << " " << j;
    }

    std::cout << " " << std::endl;
  }
}

int LL1Outv1::isValid() const
{
  return 0;
}

bool LL1Outv1::passesTrigger()
{
  for (unsigned int& _trigger_bit : *_trigger_bits)
  {
    if (_trigger_bit)
    {
      return true;
    }
  }
  return false;
}
