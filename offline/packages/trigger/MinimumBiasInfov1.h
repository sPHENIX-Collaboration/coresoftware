#ifndef MINIMUMBIASINFOV1_H
#define MINIMUMBIASINFOV1_H

#include "MinimumBiasInfo.h"

#include <iostream>

class MinimumBiasInfov1 : public MinimumBiasInfo
{
 public:
  MinimumBiasInfov1() {}
  virtual ~MinimumBiasInfov1() override {}

  void identify(std::ostream &os = std::cout) const override;

  virtual void Reset() override {};

  int isValid() const override { return 1; }


  void setIsAuAuMinimumBias(bool is_min_bias) {_isMinBias = is_min_bias;}
  bool isAuAuMinimumBias() {return _isMinBias;}

 private:

  bool _isMinBias = false;

  ClassDefOverride(MinimumBiasInfov1, 3)
};

#endif
