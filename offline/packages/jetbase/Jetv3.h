/*!
 * \file    Jetv3.h
 * \brief   Same as Jetv2, plus the corrected jet pT (MC JES calibration, in-situ JES on data, in-situ
 *          JER smearing on simulation), so that the correction modules write into the jet itself
 *          instead of creating a new jet node. Stored in JetContainerv2.
 */

#ifndef JETBASE_JETV3_H
#define JETBASE_JETV3_H

#include "Jet.h"

#include <array>
#include <cmath>
#include <cstddef>  // for size_t
#include <iostream>
#include <limits>
#include <map>
#include <utility>  // for pair, make_pair

class PHObject;

/*!
 * \brief Jetv3
 *
 * Derives from the Jet interface only (no inheritance from Jetv2: class inheritance has caused
 * readback problems before). get_pt() stays the reconstructed (raw) pT; the corrected values are
 * separate members, NaN until a module fills them:
 *   pt_calib               after the MC JES calibration (JetCalib)
 *   pt_insitu              after the in-situ JES calibration (data)
 *   pt_insitu_up / _down   its systematic variations
 *   pt_smear               after the in-situ JER smearing (simulation)
 *   pt_smear_up / _down    its systematic variations
 */
class Jetv3 : public Jet
{
 public:
  Jetv3() = default;
  explicit Jetv3(unsigned int);

  // method used to sort the consistuents of the jet -- to be used when generating the jets only
  struct CompareSRC
  {
    bool operator()(const std::pair<Jet::SRC, int>& lhs, const unsigned int rhs)
    {
      return static_cast<unsigned int>(lhs.first) < rhs;
    }
    bool operator()(const unsigned int lhs, const std::pair<Jet::SRC, int>& rhs)
    {
      return lhs < static_cast<unsigned int>(rhs.first);
    }
    bool operator()(const std::pair<Jet::SRC, int>& lhs, const std::pair<Jet::SRC, int>& rhs)
    {
      return static_cast<unsigned int>(lhs.first) < static_cast<unsigned int>(rhs.first);
    }
  };

  // PHObject virtual overloads

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;
  int isValid() const override;
  PHObject* CloneMe() const override;

  // jet info

  unsigned int get_id() const override { return _id; }
  void set_id(unsigned int id) override { _id = id; }

  float get_px() const override { return _mom[0]; }
  void set_px(float px) override { _mom[0] = px; }

  float get_py() const override { return _mom[1]; }
  void set_py(float py) override { _mom[1] = py; }

  float get_pz() const override { return _mom[2]; }
  void set_pz(float pz) override { _mom[2] = pz; }

  float get_e() const override { return _e; }
  void set_e(float e) override { _e = e; }

  // calorimeter energy fractions (layer energy / jet total energy)
  float get_emcal_frac() const override { return _emcal_frac; }
  void set_emcal_frac(float frac) override { _emcal_frac = frac; }

  float get_ihcal_frac() const override { return _ihcal_frac; }
  void set_ihcal_frac(float frac) override { _ihcal_frac = frac; }

  float get_ohcal_frac() const override { return _ohcal_frac; }
  void set_ohcal_frac(float frac) override { _ohcal_frac = frac; }

  // corrected jet pT (new in v3)
  float get_pt_calib() const override { return _pt_calib; }
  void set_pt_calib(float pt) override { _pt_calib = pt; }

  float get_pt_insitu() const override { return _pt_insitu; }
  void set_pt_insitu(float pt) override { _pt_insitu = pt; }

  float get_pt_insitu_up() const override { return _pt_insitu_up; }
  void set_pt_insitu_up(float pt) override { _pt_insitu_up = pt; }

  float get_pt_insitu_down() const override { return _pt_insitu_down; }
  void set_pt_insitu_down(float pt) override { _pt_insitu_down = pt; }

  float get_pt_smear() const override { return _pt_smear; }
  void set_pt_smear(float pt) override { _pt_smear = pt; }

  float get_pt_smear_up() const override { return _pt_smear_up; }
  void set_pt_smear_up(float pt) override { _pt_smear_up = pt; }

  float get_pt_smear_down() const override { return _pt_smear_down; }
  void set_pt_smear_down(float pt) override { _pt_smear_down = pt; }

  float get_eta() const override;

  int get_isCalib() override { return _isCalib; }
  void set_isCalib(int calib) override { _isCalib = calib; }

  float get_p() const override;
  float get_pt() const override;
  float get_et() const override;
  float get_phi() const override;
  float get_mass() const override;
  float get_mass2() const override;

  // Jet properties
  void resize_properties(size_t size) override { _properties.resize(size, std::numeric_limits<float>::quiet_NaN()); };
  std::vector<float>& get_property_vec() override { return _properties; }
  size_t size_properties() const override { return _properties.size(); }

  float get_property(Jet::PROPERTY index) const override { return _properties[static_cast<int>(index)]; };
  inline void set_property(Jet::PROPERTY index, float value) override
  {
    _properties[static_cast<int>(index)] = value;
  };

  // Jet components
  size_t size_comp() const override { return _comp_ids.size(); }
  void clear_comp() override { _comp_ids.clear(); }
  void insert_comp(SRC iSRC, unsigned int compid) override;
  void insert_comp(Jet::SRC, unsigned int compid, bool) override;  // skips setting _is_sorted flag
  void insert_comp(TYPE_comp_vec&) override;
  void insert_comp(TYPE_comp_vec&, bool) override;
  void set_comp_sort_flag(bool f = false) override { _is_sorted = f; };

  void print_comp(std::ostream& os = std::cout, bool single_line = false) override;
  size_t num_comp(SRC iSRC = Jet::SRC::VOID) override;
  std::vector<Jet::SRC> comp_src_vec() override;
  std::map<Jet::SRC, size_t> comp_src_sizemap() override;  // map of Jet::SRC to number of entries

  // FYI: ITER_comp_vec = vector<pair<Jet::SRC, unsigned int>>::iterator
  ITER_comp_vec comp_begin() override { return _comp_ids.begin(); }
  ITER_comp_vec comp_begin(Jet::SRC) override;
  ITER_comp_vec comp_end() override { return _comp_ids.end(); }
  ITER_comp_vec comp_end(Jet::SRC) override;
  TYPE_comp_vec& get_comp_vec() override { return _comp_ids; };

  inline void Clear(Option_t* = nullptr) override { Reset(); }

 private:
  /// unique identifier within container
  unsigned int _id = ~0x0;
  void ensure_sorted();
  bool _is_sorted{false};

  /// jet momentum vector (px,py,pz)
  std::array<float, 3> _mom{{std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN()}};

  /// jet energy
  float _e{std::numeric_limits<float>::quiet_NaN()};

  /// calorimeter energy fractions (layer energy / jet total energy)
  float _emcal_frac{std::numeric_limits<float>::quiet_NaN()};
  float _ihcal_frac{std::numeric_limits<float>::quiet_NaN()};
  float _ohcal_frac{std::numeric_limits<float>::quiet_NaN()};

  /// corrected jet pT
  float _pt_calib{std::numeric_limits<float>::quiet_NaN()};
  float _pt_insitu{std::numeric_limits<float>::quiet_NaN()};
  float _pt_insitu_up{std::numeric_limits<float>::quiet_NaN()};
  float _pt_insitu_down{std::numeric_limits<float>::quiet_NaN()};
  float _pt_smear{std::numeric_limits<float>::quiet_NaN()};
  float _pt_smear_up{std::numeric_limits<float>::quiet_NaN()};
  float _pt_smear_down{std::numeric_limits<float>::quiet_NaN()};

  /// calibration flag
  int _isCalib{0};
  /// source id -> component id
  std::vector<std::pair<Jet::SRC, unsigned int>> _comp_ids;

  std::vector<float> _properties{};

  bool empty_comp() const override;
  size_t count_comp(SRC source /**/) const override;

  // only in v1 msg
  static void msg_dep_fn(const std::string& fn_name);

  // functions deprecated in this Jet Version
  bool has_property(Jet::PROPERTY /*prop_id*/) const override;
  void print_property(std::ostream& /**/) const override;
  ConstIter begin_comp() const override;
  ConstIter lower_bound_comp(SRC source) const override;
  ConstIter upper_bound_comp(SRC source) const override;
  ConstIter find(Jet::SRC source) const override;
  ConstIter end_comp() const override;

  Iter begin_comp() override;
  Iter lower_bound_comp(SRC source) override;
  Iter upper_bound_comp(SRC source) override;
  Iter find(SRC source) override;
  Iter end_comp() override;

  size_t erase_comp(Jet::SRC /**/) override;
  void erase_comp(Iter /*iter*/) override;
  void erase_comp(Iter /*first*/, Iter /*last*/) override;

  ClassDefOverride(Jetv3, 1);
};

#endif  // JETBASE_JETV3_H
