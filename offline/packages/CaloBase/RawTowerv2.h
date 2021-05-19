#ifndef CALOBASE_RawTowerv2_H
#define CALOBASE_RawTowerv2_H

#include "RawTowerv1.h"

#include "RawTowerDefs.h"

#include <cstddef>
#include <iostream>
#include <map>
#include <utility>

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
  virtual ~RawTowerv2() {}

  void Reset();
  int isValid() const;
  void identify(std::ostream& os = std::cout) const;

  double get_scint_gammas() const { return get_property(prop_scint_gammas); }
  void set_scint_gammas(const double e) { set_property(prop_scint_gammas, e); }
  double get_cerenkov_gammas() const { return get_property(prop_cerenkov_gammas); }
  void set_cerenkov_gammas(const double e) { set_property(prop_cerenkov_gammas, e); }

  bool has_property(const PROPERTY prop_id) const;
  double get_property(const PROPERTY prop_id) const;
  void set_property(const PROPERTY prop_id, const double value);

 protected:
  typedef uint8_t prop_id_t;
  typedef std::map<prop_id_t, double> prop_map_t;
  //! container for additional property
  prop_map_t prop_map;

  ClassDef(RawTowerv2, 1)
};

#endif
