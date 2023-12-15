/*!
 * \file    Jetv2.h
 * \brief   Version of Jet.h updated to include data in vector<float> instead of map<int,float>
 * \author  David Stewart <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date    $Date: $
 */

#ifndef G4JET_JETV2_H
#define G4JET_JETV2_H

#include "Jet.h"

#include <cmath>
#include <cstddef>  // for size_t
#include <iostream>
#include <map>
#include <array>
#include <utility>  // for pair, make_pair

class PHObject;

/*!
 * \brief Jetv2
 */
class Jetv2 : public Jet
{
 public:
  Jetv2() = default;
  Jetv2(unsigned int);
  /* ~Jetv2() override = default; //{} */
  /* Jetv2(const Jetv2&) = default; */
  /* Jetv2(Jetv2&) = default; */
  /* Jetv2& operator=(const Jetv2&) = default; */

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
    /* static void sort_comp_ids(Jetv2* jet); */
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

  float get_p() const override;
  float get_pt() const override;
  float get_et() const override;
  float get_eta() const override;
  float get_phi() const override;
  float get_mass() const override;
  float get_mass2() const override;

  // Jet properties
  void resize_properties(size_t size) override { _properties.resize(size, NAN); };
  std::vector<float>& get_property_vec() override { return _properties; } // new in v2
  size_t size_properties() const override { return _properties.size(); }  // implemented in v1 and v2
                                                                                                       
  float get_property(Jet::PROPERTY index) const override { return _properties[static_cast<int>(index)]; };
  inline void set_property(Jet::PROPERTY index, float value) override
  { _properties[static_cast<int>(index)] = value; };

  // Jet components
  size_t size_comp() const override { return _comp_ids.size(); }
  void clear_comp() override { _comp_ids.clear(); }
  void insert_comp(SRC iSRC, unsigned int compid) override;
  void insert_comp(Jet::SRC, unsigned int compid, bool) override; //skips setting _is_sorted flag
  void insert_comp(TYPE_comp_vec&) override;
  void insert_comp(TYPE_comp_vec&, bool) override;
  void set_comp_sort_flag(bool f=false) override {_is_sorted=f;}; 

  void print_comp(std::ostream& os = std::cout, bool single_line = false) override;
  size_t num_comp(SRC iSRC = Jet::SRC::VOID) override;
  std::vector<Jet::SRC> comp_src_vec() override;
  std::map<Jet::SRC, size_t> comp_src_sizemap() override;  // map of Jet::SRC to number of entries

  size_t n_clustered() const override { return _n_clustered; }
  virtual void  set_n_clustered(unsigned int n) override { _n_clustered =n;} ;

  // FYI: ITER_comp_vec = vector<pair<Jet::SRC, unsigned int>>::iterator
  ITER_comp_vec comp_begin() override { return _comp_ids.begin(); }  // new in v2
  ITER_comp_vec comp_begin(Jet::SRC) override;                       // new in v2
  ITER_comp_vec comp_end() override { return _comp_ids.end(); }      // new in v2
  ITER_comp_vec comp_end(Jet::SRC) override;                         // new in v2
  TYPE_comp_vec& get_comp_vec() override { return _comp_ids; };      // new in v2

  inline void Clear(Option_t* = nullptr) override { Reset(); }

 private:
  /// unique identifier within container
  unsigned int _id = ~0x0;
  size_t _n_clustered {0};
  void ensure_sorted();
  bool _is_sorted { false };

  /// jet momentum vector (px,py,pz)
  std::array<float,3> _mom {{ NAN, NAN, NAN}} ;

  /// jet energy
  float _e = NAN;

  /// source id -> component id
  /* typ_comp_ids _comp_ids; */
  std::vector<std::pair<Jet::SRC, unsigned int>> _comp_ids;

  std::vector<float> _properties{};

  bool empty_comp() const override;
  size_t count_comp(SRC source /**/) const override;

  // only in v1 msg
  void msg_dep_fn(const std::string& method_name) const;

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



  ClassDefOverride(Jetv2, 1);
};

#endif  // G4JET_JETV2_H
