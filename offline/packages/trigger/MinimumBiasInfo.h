#ifndef MINIMUMBIASINFO_H
#define MINIMUMBIASINFO_H

#include <phool/PHObject.h>

class MinimumBiasInfo : public PHObject
{
 public:
  ~MinimumBiasInfo() override{};

  void identify(std::ostream &os = std::cout) const override { os << "MinimumBiasInfo base class" << std::endl; };
  virtual void Reset() override {}
  int isValid() const override { return 0; }

 protected:

  MinimumBiasInfo() {}

 private:
  ClassDefOverride(MinimumBiasInfo, 1);
};

#endif  // TRIGGER_MINBIASTRIGGERINFO_H
