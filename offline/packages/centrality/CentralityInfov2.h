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
  void setMinBias(bool isminbias) override { _isMinBias = isminbias;}

  float getVertex() const override { return _vertex;}
  void setVertex(float vertex) override {_vertex = vertex;}
  
 private:
  bool _isMinBias = false;
  float _vertex = std::numeric_limits<float>::signaling_NaN();

  ClassDefOverride(CentralityInfov2, 1);
};

#endif
