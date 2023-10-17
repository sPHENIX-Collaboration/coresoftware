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

 private:

  ClassDefOverride(CentralityInfov2, 1);
};

#endif
