#include "TriggerRunInfov1.h"

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

void TriggerRunInfov1::setTriggerScalers(int index, int scalertype, uint64_t scalers)
{
  if (index >= 0 && index < 64 && scalertype >= 0 && scalertype < 3)
  {
    trigger_scalers[index][scalertype] = scalers;
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

std::string TriggerRunInfov1::getTriggerName(int triggerbit) const {
  if (triggerbit >= 0 && triggerbit < 64) {
    return trigger_names[triggerbit];
  }
  return "unknown";
}

uint64_t TriggerRunInfov1::getScalersByName(const std::string& name) const
{
  for (int i = 0; i < 64; ++i)
  {
    if (trigger_names[i] == name)
    {
      return trigger_scalers[i][0];
    }
  }
  std::cerr << "Trigger name not found: " << name << std::endl;
  return 0;
}

uint64_t TriggerRunInfov1::getScalersByBit(int triggerbit) const {
  if (triggerbit >= 0 && triggerbit < 64) {
    return trigger_scalers[triggerbit][0];
  }
  return -1;
}
uint64_t TriggerRunInfov1::getLiveScalersByName(const std::string& name) const
{
  for (int i = 0; i < 64; ++i)
  {
    if (trigger_names[i] == name)
    {
      return trigger_scalers[i][1];
    }
  }
  std::cerr << "Trigger name not found: " << name << std::endl;
  return 0;
}

uint64_t TriggerRunInfov1::getLiveScalersByBit(int triggerbit) const {
  if (triggerbit >= 0 && triggerbit < 64) {
    return trigger_scalers[triggerbit][1];
  }
  return -1;
}


uint64_t TriggerRunInfov1::getRawScalersByName(const std::string& name) const
{
  for (int i = 0; i < 64; ++i)
  {
    if (trigger_names[i] == name)
    {
      return trigger_scalers[i][2];
    }
  }
  std::cerr << "Trigger name not found: " << name << std::endl;
  return 0;
}

uint64_t TriggerRunInfov1::getRawScalersByBit(int triggerbit) const {
  if (triggerbit >= 0 && triggerbit < 64) {
    return trigger_scalers[triggerbit][2];
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
       << ", Prescale = " << trigger_prescales[i]
       << ", Scalers = " << trigger_scalers[i][0] << "/" << trigger_scalers[i][1] << "/" << trigger_scalers[i][2] << std::endl;
  }
}
