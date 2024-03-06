
#include "LL1ReturnCodes.h"
#include "LL1Outv1.h"

#include <cmath>
#include <iostream>

ClassImp(LL1Outv1)

LL1Outv1::LL1Outv1()
{
  Init();
}

LL1Outv1::~LL1Outv1()
{

}

//______________________________________
void LL1Outv1::Init()
{
  _trigger_type = "NONE";

}


//______________________________________
void LL1Outv1::Reset()
{
  Init();
}

//______________________________________
void LL1Outv1::identify(std::ostream& out) const
{
  out << "identify yourself: I am a LL1Out object" << std::endl;
  out << "triggertype: " << _trigger_type << std::endl;
}

int LL1Outv1::isValid() const 
{

  return 1;
}


void LL1Outv1::AddTriggerBits(unsigned int t_bits)
{
  _trigger_bits.push_back(t_bits);
}
unsigned int LL1Outv1::GetTriggerBits(unsigned int i_sample)
{
  if (_trigger_bits.size() <= i_sample) return 0;
  return _trigger_bits[i_sample];
}

void LL1Outv1::AddTriggerWords(vector<unsigned int> t_words)
{
  _trigger_words.push_back(t_words);
}

unsigned int LL1Outv1::GetTriggerWord(unsigned int i_word, unsigned int i_sample)
{
  if (_trigger_words.size() <= i_sample) return 0;
  if (_trigger_words[i_sample].size() <= i_word) return 0;
  return _trigger_words[i_sample][i_word];
}

void LL1Outv1::AddTriggerPrimitives(vector<unsigned int> t_prims)
{
  _trigger_prims.push_back(t_prims);
}

unsigned int LL1Outv1::GetTriggerPrimitive( unsigned int i_adc, unsigned int i_sample)
{
  if (_trigger_prims.size() <= i_adc) return 0;
  if (_trigger_prims[i_sample].size() <= i_sample) return 0;
  return _trigger_prims[i_sample][i_adc];
}
