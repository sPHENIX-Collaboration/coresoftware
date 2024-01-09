#ifndef TRIGGER_MINIMUMBIASINFO_H
#define TRIGGER_MINIMUMBIASINFO_H

#include <phool/PHObject.h>

class MinimumBiasInfo : public PHObject
{
 public:
  ~MinimumBiasInfo() override{};

  void identify(std::ostream &os = std::cout) const override { os << "MinimumBiasInfo base class" << std::endl; };
  void Reset() override {}
  int isValid() const override { return 0; }
  virtual void setIsAuAuMinimumBias(bool) { return; }
  virtual bool isAuAuMinimumBias() const { return false; }

 protected:
  MinimumBiasInfo() {}

 private:
  ClassDefOverride(MinimumBiasInfo, 1);
};

#endif  // TRIGGER_MINIMUMBIASINFO_H
