#ifndef CENTRALITY_IO_CENTRALITYINFOV2_H
#define CENTRALITY_IO_CENTRALITYINFOV2_H

#include "CentralityInfov1.h"

#include <iostream>
#include <map>

class CentralityInfov2 : public CentralityInfov1
{
 public:
  CentralityInfov2();
  ~CentralityInfov2() override {}

  void identify(std::ostream &os = std::cout) const override;
  void Reset() override {}

  bool isMinBias() { return _isMinBias; }
  void setMinBias(bool isminbias) { _isMinBias = isminbias;}
  
 private:
  bool _isMinBias();
  ClassDefOverride(CentralityInfov2, 1);
};

#endif
