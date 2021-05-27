#ifndef CALOBASE_RawTowerv2_H
#define CALOBASE_RawTowerv2_H

#include "RawTowerv1.h"

#include "RawTowerDefs.h"

#include <cstdint>  // for uint8_t
#include <iostream>
#include <map>

//! RawTowerv1 but allow flexible tags
class RawTowerv2 : public RawTowerv1
{
 public:
  RawTowerv2();
  RawTowerv2(const RawTower& tower);
  RawTowerv2(RawTowerDefs::keytype id);
  RawTowerv2(const unsigned int ieta, const unsigned int iphi);
  RawTowerv2(const RawTowerDefs::CalorimeterId caloid, const unsigned int ieta,
             const unsigned int iphi);
  ~RawTowerv2() override {}

  void Reset() override;
  int isValid() const override;
  void identify(std::ostream& os = std::cout) const override;

  double get_scint_gammas() const override { return get_property(prop_scint_gammas); }
  void set_scint_gammas(const double e) override { set_property(prop_scint_gammas, e); }
  double get_cerenkov_gammas() const override { return get_property(prop_cerenkov_gammas); }
  void set_cerenkov_gammas(const double e) override { set_property(prop_cerenkov_gammas, e); }

  bool has_property(const PROPERTY prop_id) const override;
  double get_property(const PROPERTY prop_id) const override;
  void set_property(const PROPERTY prop_id, const double value) override;

 protected:
  typedef uint8_t prop_id_t;
  typedef std::map<prop_id_t, double> prop_map_t;
  //! container for additional property
  prop_map_t prop_map;

  ClassDefOverride(RawTowerv2, 1)
};

#endif
