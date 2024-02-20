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

  bool has_centrality_bin(const PROP prop_id) const override;
  int get_centrality_bin(const PROP prop_id) const override;
  void set_centrality_bin(const PROP prop_id, const int value) override;
 private:

  std::map<int, int> _centrality_bin_map;

  ClassDefOverride(CentralityInfov2, 1);
};

#endif
