#ifndef JETBASE_JET_H
#define JETBASE_JET_H

// The constituents between Jetv1 and Jetv2 are updated from
//    typedef multimap<Jet::SRC, unsigned int> typ_comp_ids;
//  to
//    typedef vector<pair<Jet::SRC, unsigned int>> TYPE_comp_vec;
//
//
// Functions deprecated between v1 and above
//   v1: 
//  virtual bool   has_property(Jet::PROPERTY  /*prop_id*/) const
//  virtual float  get_property(Jet::PROPERTY  /*prop_id*/) const
//  virtual void   set_property(Jet::PROPERTY  /*prop_id*/, float   /*value*/)
//  virtual void   print_property(std::ostream & /*os*/) const
//  virtual bool   empty_comp() const
//  virtual size_t size_comp()  const
//  virtual size_t count_comp(Jet::SRC /*source*/) const
//  virtual size_t erase_comp(Jet::SRC)
//  virtual void   erase_comp(Iter /*iter*/)
//  virtual void   erase_comp(Iter /*first*/, Iter /*last*/)
//
//  v2:
//  virtual size_t        n_properties()
//  virtual inline float  get_prop_by_index(unsigned int /*index*/) const
//  virtual inline void   set_prop_by_index(unsigned int /*index*/, float /*value*/)
//  virtual size_t        num_comp(SRC = Jet::SRC::VOID /**/)
//  virtual std::map<Jet::SRC, size_t> comp_src_sizemap()

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
  virtual void set_id(unsigned int) { return; }

  virtual float get_px() const { return NAN; }
  virtual void set_px(float) { return; }

  virtual float get_py() const { return NAN; }
  virtual void set_py(float) { return; }

  virtual float get_pz() const { return NAN; }
  virtual void set_pz(float) { return; }

  virtual float get_e() const { return NAN; }
  virtual void set_e(float) { return; }

  virtual float get_p() const { return NAN; }
  virtual float get_pt() const { return NAN; }
  virtual float get_et() const { return NAN; }
  virtual float get_eta() const { return NAN; }
  virtual float get_phi() const { return NAN; }
  virtual float get_mass() const { return NAN; }
  virtual float get_mass2() const { return NAN; }

  // --------------------------------------------------------------------------
  // Functions for jet properties (always float values)
  // --------------------------------------------------------------------------
  //    new with Jetv2
  virtual void resize_properties(size_t /**/) { };
  virtual std::vector<float>& get_vec_properties();
  virtual size_t n_properties() { return 0; };
  virtual inline float get_prop_by_index(unsigned int /*index*/) const { return NAN; }
  virtual inline void set_prop_by_index(unsigned int /*index*/, float /*value*/) { return; }
  //   deprecated by Jetv2
  virtual bool has_property(Jet::PROPERTY /*prop_id*/) const { return false; };
  virtual float get_property(Jet::PROPERTY /*prop_id*/) const { return NAN; };
  virtual void set_property(Jet::PROPERTY /*prop_id*/, float /*value*/) {};
  virtual void print_property(std::ostream& /*os*/) const {};
  //----------------------------------------------------------------------------------

  // --------------------------------------------------------------------------
  // Functions for jet components
  //    in all Jet versions 
  virtual void clear_comp() {}
  virtual void insert_comp(Jet::SRC, unsigned int) {}
  virtual size_t size_comp() const { return 0; };
  //    new with Jetv2
  virtual size_t num_comp(SRC = Jet::SRC::VOID /**/) { return 0; };
  virtual void print_comp(std::ostream& /**/, bool /**/) {};
  virtual std::vector<Jet::SRC> comp_src_vec() { return {}; };
  virtual std::map<Jet::SRC, size_t> comp_src_sizemap() { return {}; };
//
  typedef std::pair<Jet::SRC, unsigned int> TYPE_comp;
  typedef std::vector<TYPE_comp> TYPE_comp_vec;
  typedef TYPE_comp_vec::iterator ITER_comp_vec;
//
  virtual ITER_comp_vec comp_begin();
  virtual ITER_comp_vec comp_begin(Jet::SRC);
  virtual ITER_comp_vec comp_end();
  virtual ITER_comp_vec comp_end(Jet::SRC);
//
  virtual TYPE_comp_vec& get_comp_vec();
  //-- deprecated with Jetv2 ---------------------------------------------------------
  virtual bool empty_comp() const { return true; }
  virtual size_t count_comp(Jet::SRC /*source*/) const { return 0; };
//
  typedef std::multimap<Jet::SRC, unsigned int> typ_comp_ids;
  typedef typ_comp_ids::const_iterator ConstIter;
  typedef typ_comp_ids::iterator Iter;
//
  virtual ConstIter begin_comp() const;
  virtual ConstIter lower_bound_comp(Jet::SRC source) const;
  virtual ConstIter upper_bound_comp(Jet::SRC source) const;
  virtual ConstIter find(Jet::SRC source) const;
  virtual ConstIter end_comp() const;
//
  virtual Iter begin_comp();
  virtual Iter lower_bound_comp(Jet::SRC source);
  virtual Iter upper_bound_comp(Jet::SRC source);
  virtual Iter find(Jet::SRC source);
  virtual Iter end_comp();

  virtual size_t erase_comp(Jet::SRC) { return 0; } 
  virtual void erase_comp(Iter /*iter*/) { return; }
  virtual void erase_comp(Iter /*first*/, Iter /*last*/) { return; }
  //----------------------------------------------------------------------------------

  Bool_t IsSortable() const override { return false; }

  // structure to iterate over the jets in a TClonesArray in the JetContainer
  struct IterJetTCA
  {
    TClonesArray* tca{nullptr};
    Jet*& current_jet;  // note this is a reference to the current_jet pointer in JetContainer
    int index{0};
    int size;

    // build Iterator -- capture reference to current_jet pointer from JetContainer
    IterJetTCA(TClonesArray* _tca, Jet*& _in_jet);
    void operator++();
    Jet* operator*() { return current_jet; };
    bool operator!=(const IterJetTCA& rhs);
  };

  ClassDefOverride(Jet, 1);
};

#endif
