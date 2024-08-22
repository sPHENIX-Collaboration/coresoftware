#include "TriggerRunInfov1.h"

TriggerRunInfov1::TriggerRunInfov1()
  : trigger_names{}
  , trigger_bits{}
  , trigger_prescales{}
{
}

void TriggerRunInfov1::setTrigger(int index, const std::string& name, int bit, int prescale)
{
  if (index >= 0 && index < 64)
  {
    trigger_names[index] = name;
    trigger_bits[index] = bit;
    trigger_prescales[index] = prescale;
  }
  else
  {
    std::cerr << "Index out of bounds: " << index << std::endl;
  }
}

int TriggerRunInfov1::getPrescaleByName(const std::string& name) const
{
  for (int i = 0; i < 64; ++i)
  {
    if (trigger_names[i] == name)
    {
      return trigger_prescales[i];
    }
  }
  std::cerr << "Trigger name not found: " << name << std::endl;
  return 0;
}

int TriggerRunInfov1::getPrescaleByBit(int triggerbit) const {
  if (triggerbit >= 0 && triggerbit < 64) {
    return trigger_prescales[triggerbit];
  }
  return -1;
}
uint32_t TriggerRunInfov1::getTriggerBitByName(const std::string& name) const {
  for (uint32_t i = 0; i < 64; ++i) {
    if (trigger_names[i] == name) {
      return i;
    }
  }
  std::cerr << "Trigger name not found: " << name << std::endl;
  return 0;
}

void TriggerRunInfov1::identify(std::ostream& os) const 
{
  for (int i = 0; i < 64; ++i)
  {
    os << "Trigger " << i << ": Name = " << trigger_names[i]
       << ", Bit = " << trigger_bits[i]
       << ", Prescale = " << trigger_prescales[i] << std::endl;
  }
}
