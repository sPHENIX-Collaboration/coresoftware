#ifndef CENTRALITY_CENTRALITYINFOV2_H
#define CENTRALITY_CENTRALITYINFOV2_H

#include "CentralityInfov1.h"

#include <iostream>
#include <map>

class CentralityInfov2 : public CentralityInfov1
{
 public:
  CentralityInfov2() = default;
  ~CentralityInfov2() override = default;

  void identify(std::ostream &os = std::cout) const override;
  void Reset() override {}

  bool isMinBias() const override { return _isMinBias; }
  void setMinBias(bool isminbias) override { _isMinBias = isminbias; }

 private:
  bool _isMinBias = false;

  ClassDefOverride(CentralityInfov2, 1);
};

#endif
