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

#include <cstddef>  // for size_t
#include <cmath>
#include <iostream>
#include <map>
#include <utility>  // for pair, make_pair

class PHObject;

/*!
 * \brief Jetv2
 */
class Jetv2 : public Jet
{
 public:
  Jetv2();
  Jetv2(unsigned int);
  ~Jetv2() override {}
  Jetv2(const Jetv2&);

  // PHObject virtual overloads

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;
  int isValid() const override;
  PHObject* CloneMe() const override;

  // jet info

  unsigned int get_id() const override { return _id; }
  void  set_id(unsigned int id) override { _id = id; }

  float get_px() const override { return _mom[0]; }
  void  set_px(float px) override { _mom[0] = px; }

  float get_py() const override { return _mom[1]; }
  void  set_py(float py) override { _mom[1] = py; }

  float get_pz() const override { return _mom[2]; }
  void  set_pz(float pz) override { _mom[2] = pz; }

  float get_e() const override { return _e; }
  void  set_e(float e) override { _e = e; }

  float get_p() const override;
  float get_pt() const override;
  float get_et() const override;
  float get_eta() const override;
  float get_phi() const override;
  float get_mass() const override;
  float get_mass2() const override;

  // extended jet info

  bool  has_property(Jet::PROPERTY prop_id) const override;      // deprecated -- vector needs to sort
  float get_property(Jet::PROPERTY prop_id) const override;      // deprecated -- deprecated vec to sort
  void  set_property(Jet::PROPERTY index, float value) override; // deprecated -- index provided by JetContainer
  void  print_property(std::ostream& os=std::cout) const override;

  //
  // clustered component methods (multimap interface based)
  // source type id --> unique id within that storage
  //
  bool   empty_comp() const override { return _comp_ids.empty(); }
  size_t size_comp() const override { return _comp_ids.size(); }
  size_t count_comp(SRC source) const override;
  size_t cnt_comp(SRC iSRC); // new in v2

  void clear_comp() override { _comp_ids.clear(); }
  void insert_comp (SRC iSRC, unsigned int compid) override 
  { _comp_ids.push_back(std::make_pair(iSRC, compid)); } // 
  size_t erase_comp  (SRC iSRC) override; // THIS IS VERY EXPENSIVE ON VECTORS--shouldn't be used {return _comp_ids.erase(iSRC); }
  void erase_comp(Iter /*iter*/) override;                 // deprecated
  void erase_comp(Iter /*first*/, Iter /*last*/) override; // deprecated

  //----------------------------------------------------------------------------------------------
  // Begin methods deprecated from Jet.h (They are iterators for Maps instead of vectors)
  ConstIter begin_comp() const override;               // deprecated!
  ConstIter lower_bound_comp(Jet::SRC iSRC) const override; // deprecated!
  ConstIter upper_bound_comp(Jet::SRC iSRC) const override; // deprecated!
  ConstIter find(SRC iSRC) const override;             // deprecated!
  ConstIter end_comp() const override;                 // deprecated!

  Iter begin_comp() override;                     // deprecated!
  Iter lower_bound_comp(SRC iSRC) override;       // deprecated!
  Iter upper_bound_comp(SRC iSRC) override;       // deprecated!
  Iter find(SRC iSRC) override;                   // deprecated!
  Iter end_comp() override;                       // deprecated!
  //----------------------------------------------------------------------------------------------

  //
  void disable_v2_warning () { _print_v2_warning = false; } // new in v2

  // new data
  typedef std::vector<std::pair<Jet::SRC, unsigned int>> TYP_comp_vec;
  typedef TYP_comp_vec::iterator TYP_Iter_comp_vec;

  TYP_Iter_comp_vec  comp_begin() { return _comp_ids.begin(); } // new in v2
  TYP_Iter_comp_vec  comp_begin(Jet::SRC); // new in v2
  TYP_Iter_comp_vec  comp_end()   { return _comp_ids.end();   } // new in v2
  TYP_Iter_comp_vec  comp_end(Jet::SRC); // new in v2

  void resize_properties(size_t size) { _properties.resize(size, NAN); };

  std::vector<float>& get_vec_properties() { return _properties; } // new in v2
  unsigned int n_properties() { return _properties.size(); } // new in v2
  inline float get_prop_by_index(unsigned int index) const        { return _properties[index]; }; // new in v2
  inline void  set_prop_by_index(unsigned int index, float value) { _properties[index]=value; }; // new in v2

  Bool_t IsSortable() const override { return kTRUE; }

  TYP_comp_vec& get_comp_vec() { return _comp_ids; }; // new in v2

  void set_sort_criteria (Jet::SORT which_sort, bool large_to_small, unsigned int index=0);

 private:
  /// unique identifier within container
  unsigned int _id = ~0x0;

  /// jet momentum vector (px,py,pz)
  float _mom[3];

  /// jet energy
  float _e = NAN;

  /// keep track if sorted
  mutable bool _is_sorted { false }; // if the _comp_ids (component ids) are sorted

  /// source id -> component id
  /* typ_comp_ids _comp_ids; */
  std::vector<std::pair<Jet::SRC, unsigned int>> _comp_ids;
  void sort_comp_ids();

  //
  bool _print_v2_warning { true };

  struct CompareSRC {
      bool operator()(const std::pair<Jet::SRC,int> &lhs, const unsigned int rhs) {
          return lhs.first < rhs;
      }
      bool operator()(const unsigned int lhs, const std::pair<Jet::SRC,int> &rhs) {
          return lhs < rhs.first;
      }
     bool operator()(const std::pair<Jet::SRC,int> &lhs, const std::pair<Jet::SRC,int> &rhs) {
          return lhs.first < rhs.first;
      }
  };

  std::vector<float> _properties {};

  // member data to allow sorting and comparison
  Jet::SORT _which_sort { Jet::SORT::PT };
  unsigned int _isort_prop_index { 0 };
  int _sort_sign { -1 }; // sort sign -- for large to small or small to large

  bool IsEqual(const TObject* obj) const override;
  Int_t  Compare(const TObject* obj) const override;

  inline int intCompare (float a, float b) const {
    float val = _sort_sign*(a-b);
    if (val < 0) return -1;
    else if (val > 0) return 1;
    else return 0;
  }

  ClassDefOverride(Jetv2, 1);
};

#endif  // G4JET_JETV2_H
