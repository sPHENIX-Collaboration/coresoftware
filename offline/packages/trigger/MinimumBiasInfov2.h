#ifndef TRIGGER_MINIMUMBIASINFOV2_H
#define TRIGGER_MINIMUMBIASINFOV2_H

#include "MinimumBiasInfo.h"

#include <iostream>
#include <limits>

class MinimumBiasInfov2 : public MinimumBiasInfo
{
 public:
  MinimumBiasInfov2() = default;
  ~MinimumBiasInfov2() override = default;

  void identify(std::ostream &os = std::cout) const override;

  void Reset() override;

  int isValid() const override { return 1; }

  PHObject *CloneMe() const override { return new MinimumBiasInfov2(*this); }

  void CopyTo(MinimumBiasInfo *mbinfo) override;

  void setIsAuAuMinimumBias(bool is_min_bias) override { _isMinBias = is_min_bias; }
  bool isAuAuMinimumBias() const override { return _isMinBias; }

  void setBunchClass(int bunch) override { _bunchClass = bunch; }
  int getBunchClass() const override { return _bunchClass; }

  //! MBD t0 minus the run's main-bunch offset (NaN without an MBD vertex)
  void setT0Corrected(float t0) override { _t0Corrected = t0; }
  float getT0Corrected() const override { return _t0Corrected; }

 private:
  bool _isMinBias{false};
  int _bunchClass{MAIN};
  float _t0Corrected{std::numeric_limits<float>::quiet_NaN()};

  ClassDefOverride(MinimumBiasInfov2, 1)
};

#endif
