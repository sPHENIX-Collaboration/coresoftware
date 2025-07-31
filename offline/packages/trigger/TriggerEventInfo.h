#include <phool/PHObject.h>

#include <array>
#include <cstdint>
#include <iostream>
#include <string>

class TriggerRunInfo : public PHObject
{
 public:
  TriggerRunInfo()
    : trigger_names{}
    , trigger_bits{}
    , trigger_prescales{}
  {
  }

  void set(int index, const std::string& name, int bit, int prescale)
  {
    if (index >= 0 && index < 64)
    {
      trigger_names[index] = name;
      trigger_bits[index] = bit;
      trigger_prescales[index] = prescale;
    }
    else
    {
      std::cout << "Index out of bounds: " << index << std::endl;
    }
  }

  int getPrescaleByName(const std::string& name) const
  {
    for (int i = 0; i < 64; ++i)
    {
      if (trigger_names[i] == name)
      {
        return trigger_prescales[i];
      }
    }
    std::cout << "Trigger name not found: " << name << std::endl;
    return 0;
  }

  void printTriggers() const
  {
    for (int i = 0; i < 64; ++i)
    {
      std::cout << "Trigger " << i << ": Name = " << trigger_names[i]
                << ", Bit = " << trigger_bits[i]
                << ", Prescale = " << trigger_prescales[i] << std::endl;
    }
  }

 private:
  std::array<int, 64> trigger_live_bits;
  std::array<int, 64> trigger_live_bits;
};
