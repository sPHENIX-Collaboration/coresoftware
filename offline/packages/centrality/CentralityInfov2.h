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

  float getVertex() { return _vertex;}
  void setVertex(float vertex) {_vertex = vertex;}
  
 private:
  bool _isMinBias;
  float _vertex;

  ClassDefOverride(CentralityInfov2, 1);
};

#endif
