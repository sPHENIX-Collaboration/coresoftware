#include "LL1Outv1.h"
#include "LL1ReturnCodes.h"
#include "TriggerDefs.h"
#include "TriggerPrimitiveContainerv1.h"
#include "TriggerPrimitivev1.h"

#include <cmath>
#include <iostream>
#include <algorithm>

LL1Outv1::LL1Outv1()
  : m_trigger_key(TriggerDefs::getTriggerKey(TriggerDefs::GetTriggerId(m_trigger_type)))
{
  m_trigger_bits = new std::vector<unsigned int>();
}

LL1Outv1::LL1Outv1(const std::string& triggertype, const std::string& ll1type)
  : m_trigger_key(TriggerDefs::getTriggerKey(TriggerDefs::GetTriggerId(triggertype)))
  , m_triggerid(TriggerDefs::GetTriggerId(triggertype))
  , m_ll1_type(ll1type)
  , m_trigger_type(triggertype)
{
  m_trigger_bits = new std::vector<unsigned int>();

  int ntriggerwords = 0;
  if (m_triggerid == TriggerDefs::TriggerId::jetTId || m_triggerid == TriggerDefs::TriggerId::photonTId)
    {
      ntriggerwords = 384;
    }
  else if (m_triggerid == TriggerDefs::TriggerId::pairTId )
    {
      ntriggerwords = 0;
    }
  else if (m_triggerid == TriggerDefs::TriggerId::mbdTId )
    {
      ntriggerwords = 8;
    }

  for (int channel = 0; channel < ntriggerwords; channel++)
    {
      std::vector<unsigned int>* sum = new std::vector<unsigned int>();
      if (m_triggerid == TriggerDefs::TriggerId::jetTId || m_triggerid == TriggerDefs::TriggerId::photonTId)
	{
	  add_word(((unsigned int) (channel % 32) & 0xffffU) + (((unsigned int) (channel / 32) & 0xffffU) << 16U), sum);
	}
      if (m_triggerid == TriggerDefs::TriggerId::mbdTId)
	{
	  add_word(channel, sum);
	}

    }
}

LL1Outv1::~LL1Outv1()
{
  // you cannot call a virtual method in the ctor, you need to be specific
  // which one you want to call

  m_trigger_bits->clear();
  m_triggered_sums.clear();
  m_triggered_primitives.clear();
  for (auto &word : m_trigger_words)
  {
    word.second->clear();
  }
  
}

//______________________________________
void LL1Outv1::Reset()
{
  m_trigger_bits->clear();
  m_triggered_sums.clear();
  m_triggered_primitives.clear();
  for (auto &word : m_trigger_words)
  {
    word.second->clear();
  }
}

LL1Outv1::ConstRange LL1Outv1::getTriggerWords() const
{
  return make_pair(m_trigger_words.begin(), m_trigger_words.end());
}

LL1Outv1::Range LL1Outv1::getTriggerWords()
{
  return make_pair(m_trigger_words.begin(), m_trigger_words.end());
}

//______________________________________
void LL1Outv1::identify(std::ostream& out) const
{
  out << __FILE__ << __FUNCTION__ << " LL1Out: Triggertype = " << m_trigger_type << " LL1 type = " << m_ll1_type << std::endl;
  out << __FILE__ << __FUNCTION__ << " Event number: " << m_event_number << "    Clock: " << m_clock_number << std::endl;
  out << __FILE__ << __FUNCTION__ << " Trigger bits    " << m_trigger_bits->size() << " : ";
  for (unsigned int trigger_bit : *m_trigger_bits)
  {
    std::cout << " " << trigger_bit;
  }
  out << " " << std::endl;
  out << __FILE__ << __FUNCTION__ << " Trigger words    " << std::endl;

  for (auto trigger_word : m_trigger_words)
  {
    for (unsigned int& j : *trigger_word.second)
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
  for (unsigned int& trigger_bit : *m_trigger_bits)
  {
    if (trigger_bit)
    {
      return true;
    }
  }
  return false;
}

bool LL1Outv1::passesThreshold(int ith)
{
  for (unsigned int& trigger_bit : *m_trigger_bits)
  {
    if (((trigger_bit >> (uint16_t) (ith - 1)) & 0x1U) == 0x1U)
    {
      return true;
    }
  }
  return false;
}

void LL1Outv1::addTriggeredSum(TriggerDefs::TriggerSumKey sk) 
{
  unsigned int sumk = sk;
  if (!m_triggered_sums.size())
    {
      m_triggered_sums.push_back(sumk);
      return;
    }
  if (std::find(m_triggered_sums.begin(), m_triggered_sums.end(), sumk) == std::end(m_triggered_sums))
    {
      m_triggered_sums.push_back(sumk);
    }
}
void LL1Outv1::addTriggeredPrimitive(TriggerDefs::TriggerPrimKey pk)
{
  unsigned int primk = pk;
  if (!m_triggered_primitives.size())
    {
      m_triggered_primitives.push_back(primk);
      return;
    }
  if (std::find(m_triggered_primitives.begin(), m_triggered_primitives.end(), primk) == std::end(m_triggered_primitives))
    {
      m_triggered_primitives.push_back(primk);
    }
  return;
}
