#ifndef G4JET_JET_H
#define G4JET_JET_H

#include <phool/PHObject.h>

#include <cmath>
#include <cstddef>  // for size_t
#include <iostream>
#include <map>

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
  };


  enum SORT  // used as criteria for sorting output in JetMap
  {
    NO_SORT = 0, // a blank input to not sort input
    PT   = 1, // PT descending order
    E    = 2, // E  descending order
    P    = 3, // P descending order
    MASS = 4, // Mass descending order
    AREA = 5, // AREA descending order --> maybe used in future, as jets don't have area for now...
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

  // extended jet info ---------------------------------------------------------

  virtual bool has_property(Jet::PROPERTY /*prop_id*/) const { return false; }
  virtual float get_property(Jet::PROPERTY /*prop_id*/) const { return NAN; }
  virtual void set_property(Jet::PROPERTY /*prop_id*/, float /*value*/) { return; }
  virtual void print_property(std::ostream& /*os*/) const { return; }

  // component id storage ------------------------------------------------------

  /*! \addtogroup clustered component
   * clustered component methods (multimap interface based)
   * source type id --> unique id within that storage
   *  @{
   */

  typedef std::multimap<Jet::SRC, unsigned int> typ_comp_ids;
  typedef typ_comp_ids::const_iterator ConstIter;
  typedef typ_comp_ids::iterator Iter;

  virtual bool empty_comp() const { return true; }
  virtual size_t size_comp() const { return 0; }
  virtual size_t count_comp(Jet::SRC /*source*/) const { return 0; }

  virtual void clear_comp() { return; }
  virtual void insert_comp(Jet::SRC /*source*/, unsigned int /*compid*/) { return; }
  virtual size_t erase_comp(Jet::SRC /*source*/) { return 0; }
  virtual void erase_comp(Iter /*iter*/) { return; }
  virtual void erase_comp(Iter /*first*/, Iter /*last*/) { return; }

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

  ClassDefOverride(Jet, 1);
};

#endif
