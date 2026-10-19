#ifndef TRIGGER_MINIMUMBIASINFO_H
#define TRIGGER_MINIMUMBIASINFO_H

#include <phool/PHObject.h>

#include <limits>

class MinimumBiasInfo : public PHObject
{
 public:

  enum SPECIES
    {
      AUAU = 0,
      OO = 1,
      PP = 2
    };

  //! Au+Au bunch class from the MBD t0 relative to the run's main-bunch offset (centrality is calibrated per
  //! class); NONE = outside the three bunch windows (no centrality)
  enum BUNCH
  {
    NONE = -1,
    MAIN = 0,
    M5 = 1,
    P5 = 2
  };

  ~MinimumBiasInfo() override {};

  void identify(std::ostream &os = std::cout) const override { os << "MinimumBiasInfo base class" << std::endl; };
  void Reset() override {}
  int isValid() const override { return 0; }
  virtual void CopyTo(MinimumBiasInfo *) { return; }
  virtual void setIsAuAuMinimumBias(bool) { return; }
  virtual bool isAuAuMinimumBias() const { return false; }
  //! v2 only; MinimumBiasInfov1 has no bunch class, so its events read as the main bunch
  virtual void setBunchClass(int) { return; }
  virtual int getBunchClass() const { return MAIN; }
  virtual void setT0Corrected(float) { return; }
  virtual float getT0Corrected() const { return std::numeric_limits<float>::quiet_NaN(); }

 protected:
  MinimumBiasInfo() {}

 private:
  ClassDefOverride(MinimumBiasInfo, 1);
};

#endif  // TRIGGER_MINIMUMBIASINFO_H
