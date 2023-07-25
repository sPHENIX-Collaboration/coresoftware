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

  float get_p()     const override;
  float get_pt()    const override;
  float get_et()    const override;
  float get_eta()   const override;
  float get_phi()   const override;
  float get_mass()  const override;
  float get_mass2() const override;

  // Jet properties
  void resize_properties(size_t size) override { _properties.resize(size, NAN); };
  std::vector<float>& get_vec_properties() override { return _properties; } // new in v2
  size_t n_properties() override { return _properties.size(); } // new in v2
  inline float get_prop_by_index(unsigned int index) const override { return _properties[index]; }; // new in v2
  inline void  set_prop_by_index(unsigned int index, float value) override { 
    _properties[index]=value; 
  }; 

  // Jet components
  void   clear_comp() override { _comp_ids.clear(); }
  void   insert_comp(SRC iSRC, unsigned int compid) override;
  void   print_comp(std::ostream& os = std::cout, bool single_line=false) override;
  size_t num_comp(SRC iSRC=Jet::SRC::VOID) override; 
  std::vector<Jet::SRC> comp_src_vec() override;
  std::map<Jet::SRC,size_t> comp_src_sizemap() override; // map of Jet::SRC to number of entries

  // FYI: ITER_comp_vec = vector<pair<Jet::SRC, unsigned int>>::iterator
  ITER_comp_vec  comp_begin() override { return _comp_ids.begin(); } // new in v2
  ITER_comp_vec  comp_begin(Jet::SRC) override;                      // new in v2
  ITER_comp_vec  comp_end()  override  { return _comp_ids.end();   } // new in v2
  ITER_comp_vec  comp_end(Jet::SRC) override;                        // new in v2
  TYPE_comp_vec& get_comp_vec()  override{ return _comp_ids; }; // new in v2
  void sort_comp_ids() override;

  Bool_t IsSortable() const override { return kTRUE; }
  void set_sort_selptr (Jet::SortSelections* _) override { _sortopt = _; }
  Jet::SortSelections* get_sort_selptr() override { return _sortopt; }

  inline void Clear(Option_t* =nullptr) override { Reset(); }

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

  //
  /* bool _print_v2_warning { true }; */

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
  Jet::SortSelections* _sortopt { nullptr }; // points to sorting options set in JetContainer

  bool  IsEqual(const TObject* obj) const override;
  Int_t Compare(const TObject* obj) const override;


  inline int intCompare (float a, float b) const {
    if (_sortopt->order == Jet::SORT_ORDER::ASCENDING) {
      if      (a<b) return -1; // will fail for either number NAN
      else if (b<a) return  1; // will fail for either number NAN
      else if (b==a) return 0;
      else if (isnan(a) && isnan(b)) return 0;
      else if (isnan(a)) return  1;
      else if (isnan(b)) return -1;
      else               return  0;
    } else {
      if      (a<b) return  1; // will fail for either number NAN
      else if (b<a) return -1; // will fail for either number NAN
      else if (b==a) return 0;
      else if (isnan(a) && isnan(b)) return 0;
      else if (isnan(a)) return  1;
      else if (isnan(b)) return -1;
      else               return  0;
    }
  }

  ClassDefOverride(Jetv2, 1);
};

#endif  // G4JET_JETV2_H
