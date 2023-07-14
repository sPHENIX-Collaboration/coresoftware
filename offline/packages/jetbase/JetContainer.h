#ifndef JETBASE_JETCONTAINER_H
#define JETBASE_JETCONTAINER_H

#include "Jet.h"
/* #include "Jetv2.h" */

#include <phool/PHObject.h>

#include <cmath>
#include <cfloat>
#include <cstddef>  // for size_t
#include <iostream>
#include <map>
#include <vector>
#include <functional>
#include <TClonesArray.h>
#include <set>
#include <limits.h>

class Jet;
class Jetv2;

// ---------------------------------------------------------------------------------------
// Convenience class for iterating over jets in TClonesArray in JetContainer
// ---------------------------------------------------------------------------------------
struct IterJetTCA {
  TClonesArray* tca   { nullptr };
  Jetv2*& current_jet ;
  int index { 0 };
  int size;

  // build Iterator -- capture reference to current_jet pointer from JetContainer
  IterJetTCA  (TClonesArray* _tca, Jetv2*& _in_jet) 
      : tca{_tca}, current_jet{_in_jet}, size { tca->GetEntriesFast() }
  {}

  IterJetTCA begin() { return *this; };

  IterJetTCA end() { return *this; };

  void operator++() {
      current_jet = (Jetv2*) tca->UncheckedAt(++index);
  };

  Jetv2* operator*() { return current_jet; };

  bool operator!=(const IterJetTCA& rhs) { return index!=rhs.size;};
};


// ---------------------------------------------------------------------------------------
// JetContainer class -- used to fill, update, and sort TClonesArray of jets
// ---------------------------------------------------------------------------------------
class JetContainer : public PHObject
{
public:
    ~JetContainer() = default;
    virtual void identify(std::ostream&/*-*/) const override;
    int isValid() const override { return 0; }
    PHObject* CloneMe() const override { return nullptr; }
    virtual TClonesArray* clone_data() const;

    // status of jet contents
    virtual bool   empty() const { return true; }
    virtual size_t size()  const { return 0;    }

    // adding/access jets
    virtual Jetv2* current_jet() ; // points to most recently accessed jet
    virtual Jetv2* add_jet()                  { return nullptr; }; // Add a new jet to the TClonesArray and return the pointer
    virtual Jetv2* get_jet(unsigned int /*index*/)         { return nullptr; }; // Get jet at loc. return nullptr is out of range 
    virtual Jetv2* get_UncheckedAt(unsigned int /*index*/) { return nullptr; }; // Get get at location; no range checking
                                                
    // convenience shortcuts of get_{jet,UncheckedAt}
    virtual Jetv2* operator()(int /*index*/) { return nullptr; }; // synonym for get_get()
    virtual Jetv2* operator[](int /*index*/) { return nullptr; }; // get jet, don't check for length

    // add compondent to the jet source components
    virtual void add_component(Jet::SRC /*-*/, unsigned int /*component id*/) {};


    // ---------------------------------------------------------------------------
    // Legacy functions (copied from JetMap) used to record which functions are used.
    // These options are not actually consulted when the clustering is done in the 
    // FastJetAlgo class; these are just informational.
    // ---------------------------------------------------------------------------
    virtual void set_algo(Jet::ALGO /*algo*/) { return; }
    virtual Jet::ALGO get_algo() const { return Jet::ALGO::NONE; }

    virtual void set_jetpar_R(float) { return; }
    virtual float get_jetpar_R() const { return NAN; }

    // ---------------------------------------------------------------------------
    //  Sources "src":
    // Keep a std::set of data listing the input sources
    // present in the jest
    // ---------------------------------------------------------------------------
    typedef std::set<Jet::SRC>::const_iterator ConstSrcIter;
    typedef std::set<Jet::SRC>::iterator SrcIter;

    virtual bool empty_src() const { return true; }
    virtual void insert_src(Jet::SRC /*src*/) { return; }

    virtual ConstSrcIter begin_src() const;
    virtual ConstSrcIter find_src(Jet::SRC src) const;
    virtual ConstSrcIter end_src() const;

    virtual SrcIter begin_src();
    virtual SrcIter find_src(Jet::SRC src);
    virtual SrcIter end_src();

    // ----------------------------------------------------------------------------------------
    //  Interface for adding, setting, and getting jet properties
    // ----------------------------------------------------------------------------------------
    // The optional properties in each Jetv2 are stored in a vector of floats,
    // e.g. { Rg_value, mg_value, area, .... }
    // The JetContainer generates the jet and keeps a map of which properties in which location.
    // Default values in the Jet Properties are set to NAN
    // ----------------------------------------------------------------------------------------

    // Get and queary the map of indices of the vector<properties> in the jets
    virtual std::map<Jet::PROPERTY, unsigned int> prop_indices_map() const { return {}; };
    virtual std::vector<Jet::PROPERTY> vec_jet_properties() const { return {}; }; // values from the prop_map
    virtual void print_property_types(std::ostream& /*os*/) const {}; // print the order of properties in jet
    virtual bool has_property(Jet::PROPERTY /*-*/) const { return false; }; 
    virtual size_t n_properties() const { return UINT_MAX; }; // number of properties in the jet
    virtual void print_missing_prop(Jet::PROPERTY/**/) const {}; // print a warning if the property is not in the jet
  
    virtual void print_jets(std::ostream& os=std::cout) { os<<""; return; }; // print the order of properties in jet

    // Add properties to the jets. 
    // If it is a new property(/ies) it will expand all jet prop vectors accordingly
    virtual size_t add_property(Jet::PROPERTY/**/) { return 0; } ; 
    virtual size_t add_property(std::set<Jet::PROPERTY> /**/) { return 0; }; 

    // Get and set values of properties by index (always on current_jet)
    virtual unsigned int find_prop_index (Jet::PROPERTY /*-*/) { return UINT_MAX; };
    virtual float get_prop_by_index (unsigned int /*index*/) const { return NAN; };
    virtual void  set_prop_by_index (unsigned int /*index*/, float /*value*/) {};

    // Get and set values of properties by selections (always on current_jet)
    virtual void select_property(Jet::PROPERTY/**/) {}; // select this property for the current jet
    virtual void select_property(std::vector<Jet::PROPERTY>/**/) {}; // order matters
    virtual float get_selected_property(unsigned int/*selection index=0*/) { return FLT_MAX;};
    virtual void  set_selected_property(float /*set-value*/, unsigned int/*=0*/) {};
                                               
    // ---------------------------------------------------------------------------------------
    //  Add ability for: ```for (auto jet : jet_container->iter_jets()) { ... }```
    // ---------------------------------------------------------------------------------------
    // Use:
    // ```
    //     for (jet : jet_containter->iter_jets()) { 
    //     // jet will increment as Jetv2* through all jets
    //     ... }
    // ```
    // In the loop, current_jet in JetContainer will be updated, too. So, you can reference 
    // values like get_property.
    virtual IterJetTCA iter_jets() ;
    virtual IterJetTCA begin()     ;
    virtual IterJetTCA end()       ;
    // ---------------------------------------------------------------------------------------

    // Methods to sort the TClonesArray of jets
    virtual void sort_jets(Jet::SORT /**/,     bool /*descending */) {};
    virtual void sort_jets(Jet::PROPERTY /**/, bool /*descending */) {};
    virtual bool is_sorted() const {return false;} ;
    virtual void set_is_sorted(bool /*-*/) {};

    virtual void print_sorted_by(std::ostream&/**/) {};
    virtual void set_sorted_by(Jet::SORT /**/, bool /*is_inverse*/, Jet::PROPERTY=Jet::PROPERTY::no_property) {};

    virtual unsigned int get_index_single() const { return UINT_MAX; };
    virtual std::vector<unsigned int> get_index_vec()    const { return {}; };

    virtual void set_rho_median(float /**/) {};
    virtual float get_rho_median() const { return NAN; };
  ClassDefOverride(JetContainer, 1);
};


#endif /* JETBASE_JETCONTAINER_H */
