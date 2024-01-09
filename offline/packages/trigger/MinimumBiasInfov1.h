<<<<<<< HEAD
#ifndef MINIMUMBIASINFOV1_H
#define MINIMUMBIASINFOV1_H
=======
#ifndef TRIGGER_MINIMUMBIASINFOV1_H
#define TRIGGER_MINIMUMBIASINFOV1_H
>>>>>>> b183955abcd8650f7a6403b3814a7542343be7ac

#include "MinimumBiasInfo.h"

#include <iostream>

class MinimumBiasInfov1 : public MinimumBiasInfo
{
 public:
<<<<<<< HEAD
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
=======
  MinimumBiasInfov1() = default;
  virtual ~MinimumBiasInfov1() override = default;

  void identify(std::ostream &os = std::cout) const override;

  virtual void Reset() override{};

  int isValid() const override { return 1; }

  void setIsAuAuMinimumBias(bool is_min_bias) override { _isMinBias = is_min_bias; }
  bool isAuAuMinimumBias() const override { return _isMinBias; }

 private:
  bool _isMinBias{false};

  ClassDefOverride(MinimumBiasInfov1, 1)
>>>>>>> b183955abcd8650f7a6403b3814a7542343be7ac
};

#endif
