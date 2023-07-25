#ifndef JETBASE_JET_H
#define JETBASE_JET_H

// Updates in 07/2023 by D. Stewart
// Moving the jets into TClonesArrays with Jetv2 implementation. Major changes:
// 1. Overload IsSortable, IsEqual, Compare for sake of sorting through ROOT::TClonesArrays
// 2. Internally, store jet components in vector<Jet::SRC,unsigned int>, and optional
//    jet properties in  vector<float>; instead of maps as used in Jetv1. The map of
//    jet properties (i.e. linking which property is stored in which vector location)
//    is moved into the JetContainer.
//
//    Therefore: 
//     - map accessors are deprecated, and the vector accessors are added.
//     - some other functions which are 'const' with map cannot be const, 
//       and are therefore also updated.
// 3. Teh enumerators have also been updated.

#include <phool/PHObject.h>

#include <cmath>
#include <cstddef>  // for size_t
#include <iostream>
#include <map>

class TClonesArray;

class Jet : public PHObject
{
 public:
  // enums can be extended with new values, but values not altered

  enum ALGO
  {
    NONE = 0,
    ANTIKT = 1,
    KT = 2,
    CAMBRIDGE = 3
  };

  enum SRC
  {
    VOID = 0,
    PARTICLE = 1,
    TRACK = 2,
    CEMC_TOWER = 3,
    CEMC_CLUSTER = 4,
    HCALIN_TOWER = 5,
    HCALIN_CLUSTER = 6,
    HCALOUT_TOWER = 7,
    HCALOUT_CLUSTER = 8,
    FEMC_TOWER = 9,
    FEMC_CLUSTER = 10,
    FHCAL_TOWER = 11,
    FHCAL_CLUSTER = 12,
    CEMC_TOWER_RETOWER = 13, /* needed for HI jet reco */
    CEMC_TOWER_SUB1 = 14,
    HCALIN_TOWER_SUB1 = 15,
    HCALOUT_TOWER_SUB1 = 16, /* needed for HI jet reco */
    CEMC_TOWER_SUB1CS = 17,
    HCALIN_TOWER_SUB1CS = 18,
    HCALOUT_TOWER_SUB1CS = 19, /* needed for CS subtraction w/ HI jet reco */
    HEPMC_IMPORT = 20,         /*Direct import HEPMC containers, such as sHijing HIJFRG truth jets loaded by JetHepMCLoader*/
    HCAL_TOPO_CLUSTER = 21,    /* I+HOCal 3-D topoCluster input */
    ECAL_TOPO_CLUSTER = 22,    /* EMCal 3-D topoCluster input */
    EEMC_TOWER = 23,
    EEMC_CLUSTER = 24,
    CEMC_TOWERINFO = 25,
    HCALIN_TOWERINFO = 26,
    HCALOUT_TOWERINFO = 27,
    CEMC_TOWERINFO_RETOWER = 28, /* needed for HI jet reco */
    CEMC_TOWERINFO_SUB1 = 29,
    HCALIN_TOWERINFO_SUB1 = 30,
    HCALOUT_TOWERINFO_SUB1 = 31, /* needed for HI jet reco */
  };

  enum class SORT  // used as criteria for sorting output in JetMap
  {
    NO_SORT = 0, // a blank input to not sort input
    PT   = 1, // PT descending order
    E    = 2, // E  descending order
    P    = 3, // P descending order
    MASS = 4, // Mass descending order
    AREA = 5, // AREA descending order --> maybe used in future, as jets don't have area for now...
    PROPERTY = 6, // doesn't say which property; a second parameter of
                  // Jet::PROPERTY will be needed the sorting algorithm
    MASS2 = 7, // Mass descending order
    ETA   = 8, // Mass descending order
  };

  enum class SORT_ORDER
  {
    DESCENDING = 0,
    ASCENDING  = 1,
  };

  enum PROPERTY
  {
    //! jet charge
    prop_JetCharge = 1,

    //! b-jet fraction
    prop_BFrac = 2,

    //! discriminator D = max tower E / average E , used to identify
    //! seeds in 1st iteration of UE determination
    prop_SeedD = 3,

    //! used to tag as seed jet in 1st or 2nd iteration of UE
    //! determination
    prop_SeedItr = 4,

    //! SoftDrop quantities
    prop_zg = 5,
    prop_Rg = 6,
    prop_mu = 7,
    //! photon tag property
    prop_gamma = 8,
    prop_JetHadronFlavor = 9,
    prop_JetHadronZT = 10,

    //! jet area
    prop_area = 11,

    no_property = 12,
  };

  Jet() {}
  ~Jet() override {}

  void identify(std::ostream& os = std::cout) const override;
  int isValid() const override { return 0; }
  PHObject* CloneMe() const override { return nullptr; }

  // jet info ------------------------------------------------------------------

  virtual unsigned int get_id() const { return 0xFFFFFFFF; }
  virtual void set_id(unsigned int)   { return; }

  virtual float get_px() const        { return NAN; }
  virtual void set_px(float)          { return; }

  virtual float get_py() const        { return NAN; }
  virtual void set_py(float)          { return; }

  virtual float get_pz() const        { return NAN; }
  virtual void set_pz(float)          { return; }

  virtual float get_e() const         { return NAN; }
  virtual void set_e(float)           { return; }

  virtual float get_p() const         { return NAN; }
  virtual float get_pt() const        { return NAN; }
  virtual float get_et() const        { return NAN; }
  virtual float get_eta() const       { return NAN; }
  virtual float get_phi() const       { return NAN; }
  virtual float get_mass() const      { return NAN; }
  virtual float get_mass2() const     { return NAN; }

  // messages for detaulf function implementations
  static void mv1(const std::string& method_name, std::ostream& os=std::cout); //msg for Jetv1 imp. only
  static void mv2(const std::string& method_name, std::ostream& os=std::cout); //msg for >Jetv2 imp. only

  // --------------------------------------------------------------------------
  // Functions for jet properties (always float values)
  // -- new with Jetv2 --------------------------------------------------------
  virtual void  resize_properties(size_t /**/)  { mv2("resize_properties()"); }; 
  virtual std::vector<float>& get_vec_properties();
  virtual size_t       n_properties()                  { mv2("n_properties()");      return 0; };
  virtual inline float get_prop_by_index(unsigned int /*index*/) const  { mv2("get_prop_by_index()"); return NAN; }
  virtual inline void  set_prop_by_index(unsigned int /*index*/, float /*value*/) { mv2("set_prop_by_index()"); return; }
  //-- deprecated with Jetv2 ---------------------------------------------------------
  virtual bool  has_property(Jet::PROPERTY /*prop_id*/) const { mv1("has_property()"); return false; }      //
  virtual float get_property(Jet::PROPERTY /*prop_id*/) const { mv1("get_property()"); return NAN; }       //
  virtual void  set_property(Jet::PROPERTY /*prop_id*/, float /*value*/) { mv1("set_property()"); return; } //
  virtual void  print_property(std::ostream& /*os*/) const { mv1("print_property()"); return; } //
  //----------------------------------------------------------------------------------

  
  // --------------------------------------------------------------------------
  // Functions for jet components
  // -- in all Jet versions --------------------------------------------------
  virtual void clear_comp() { }
  virtual void insert_comp(Jet::SRC, unsigned int) { }
  // -- new with Jetv2 --------------------------------------------------------
  virtual size_t num_comp(SRC=Jet::SRC::VOID /**/) { mv2("num_comp()"); return 0; }
  virtual void   print_comp(std::ostream& /**/, bool /**/) { mv2("print_comp()"); }
  virtual std::vector<Jet::SRC> comp_src_vec() { mv2("comp_src_vec()"); return {}; }
  virtual std::map<Jet::SRC,size_t> comp_src_sizemap() { mv2("comp_src_sizemap()"); return {}; }

  typedef std::vector<std::pair<Jet::SRC, unsigned int>> TYPE_comp_vec;
  typedef TYPE_comp_vec::iterator ITER_comp_vec;

  virtual ITER_comp_vec  comp_begin() ;
  virtual ITER_comp_vec  comp_begin(Jet::SRC);
  virtual ITER_comp_vec  comp_end()   ;
  virtual ITER_comp_vec  comp_end(Jet::SRC); 
  virtual TYPE_comp_vec& get_comp_vec() ;
  virtual void sort_comp_ids() { mv2("sort_comp_ids()"); return; }
  //-- deprecated with Jetv2 ---------------------------------------------------------
  virtual bool empty_comp() const  { mv1("empty_comp()"); return true; }
  virtual size_t size_comp() const { mv1("size_comp()"); return 0; }
  virtual size_t count_comp(Jet::SRC /*source*/) const  { mv1("count_comp()"); return 0; }

  typedef std::multimap<Jet::SRC, unsigned int> typ_comp_ids;
  typedef typ_comp_ids::const_iterator ConstIter;
  typedef typ_comp_ids::iterator Iter;

  virtual ConstIter begin_comp() const;
  virtual ConstIter lower_bound_comp(Jet::SRC source) const;
  virtual ConstIter upper_bound_comp(Jet::SRC source) const;
  virtual ConstIter find(Jet::SRC source) const;
  virtual ConstIter end_comp() const;

  virtual Iter begin_comp();
  virtual Iter lower_bound_comp(Jet::SRC source);
  virtual Iter upper_bound_comp(Jet::SRC source);
  virtual Iter find(Jet::SRC source);
  virtual Iter end_comp();
  virtual size_t erase_comp(Jet::SRC)    { mv1("erase_comp"); return 0; } // could implement for v2, but would be expensive
  virtual void erase_comp(Iter /*iter*/) { mv1("erase_comp"); return; }
  virtual void erase_comp(Iter /*first*/, Iter /*last*/) { mv1("erase_comp"); return; }
  //----------------------------------------------------------------------------------
 
  //-- new with Jetv2: Used for allowing TClonesArray to sort jets. Each jet will
  // point to the JetSortSels (Jet Sort Selections) object in JetContainer and then
  // the TClonesArray will sort the jets based on the criteria in the JetSortSels.
  struct SortSelections {
      Jet::SORT       criteria   { Jet::SORT::PT };
      Jet::PROPERTY   property   { Jet::PROPERTY::no_property }; // when sorted by property
      Jet::SORT_ORDER order      { Jet::SORT_ORDER::DESCENDING }; // 
      unsigned int    prop_index { 0 }; // for use when sorting by criteria
  };

  virtual Bool_t IsSortable() const override { mv2("IsSortable()"); return false; }
  virtual void set_sort_selptr (Jet::SortSelections* /*sortopt*/) { mv2("set_sort_selptr()"); }
  virtual Jet::SortSelections* get_sort_selptr() { mv2("get_sort_selptr"); return {}; }
  protected:
  virtual bool  IsEqual(const TObject* /**/) const override {mv2("IsEqual()"); return false;}
  virtual Int_t Compare(const TObject* /**/) const override {mv2("Compare()"); return 0;}

  public:
  // structure to iterate over ther jets in a TClonesArray in the JetContainer
  struct IterJetTCA {
    TClonesArray* tca   { nullptr };
    Jet*& current_jet ; // note this is areference to the current_jet pointer in JetContainer
    int index { 0 };
    int size;

    // build Iterator -- capture reference to current_jet pointer from JetContainer
    IterJetTCA  (TClonesArray* _tca, Jet*& _in_jet) ;
    void operator++();
    Jet* operator*() { return current_jet; };
    bool operator!=(const IterJetTCA& rhs);
  };

  ClassDefOverride(Jet, 1);
};

#endif
