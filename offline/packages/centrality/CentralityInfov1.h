#ifndef CENTRALITY_IO_CENTRALITYINFOV1_H
#define CENTRALITY_IO_CENTRALITYINFOV1_H

#include "CentralityInfo.h"

#include <iostream>
#include <map>

class CentralityInfov1 : public CentralityInfo
{
 public:
  CentralityInfov1();
  ~CentralityInfov1() override {}

  void identify(std::ostream &os = std::cout) const override;
  void Reset() override {}
  int isValid() const override { return 1; }

  bool has_quantity(const PROP prop_id) const override;
  float get_quantity(const PROP prop_id) const override;
  void set_quantity(const PROP prop_id, const float value) override;

  bool has_centile(const PROP prop_id) const override;
  float get_centile(const PROP prop_id) const override;
  void set_centile(const PROP prop_id, const float value) override;

 private:
  std::map<int, float> _quantity_map;
  std::map<int, float> _centile_map;

  ClassDefOverride(CentralityInfov1, 1);
};

#endif
