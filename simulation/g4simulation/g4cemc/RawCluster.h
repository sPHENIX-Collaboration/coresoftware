#ifndef RAWCLUSTER_H__
#define RAWCLUSTER_H__

#include "RawClusterDefs.h"
#include "RawTowerDefs.h"

#include <phool/PHObject.h>
#include <phool/phool.h>
#include <climits>
#include <cmath>  // def of NAN
#include <map>

class RawTower;

class RawCluster : public PHObject
{
 public:
  typedef std::map<RawTowerDefs::keytype, float> TowerMap;
  typedef TowerMap::iterator TowerIterator;
  typedef TowerMap::const_iterator TowerConstIterator;
  typedef std::pair<TowerIterator, TowerIterator> TowerRange;
  typedef std::pair<TowerConstIterator, TowerConstIterator> TowerConstRange;

  virtual ~RawCluster() {}
  virtual void Reset() { PHOOL_VIRTUAL_WARNING; }
  virtual int isValid() const
  {
    PHOOL_VIRTUAL_WARNING;
    return 0;
  }
  virtual void identify(std::ostream& os = std::cout) const { PHOOL_VIRTUAL_WARNING; }
  virtual RawClusterDefs::keytype get_id() const
  {
    PHOOL_VIRTUAL_WARN("get_id()");
    return 0;
  }

  //! convert cluster location to psuedo-rapidity with a user chosen z-location
  virtual float get_eta(const float z) const
  {
    PHOOL_VIRTUAL_WARN("get_eta()");
    return NAN;
  }

  //! location of cluster in cylindrical coordinate

  virtual float get_phi() const
  {
    PHOOL_VIRTUAL_WARN("get_phi()");
    return NAN;
  }
  virtual float get_energy() const
  {
    PHOOL_VIRTUAL_WARN("get_energy()");
    return NAN;
  }
  virtual float get_ecore() const
  {
    PHOOL_VIRTUAL_WARN("get_ecore()");
    return NAN;
  }
  virtual float get_prob() const
  {
    PHOOL_VIRTUAL_WARN("get_prob()");
    return NAN;
  }
  virtual float get_chi2() const
  {
    PHOOL_VIRTUAL_WARN("get_chi2()");
    return NAN;
  }

  virtual void set_id(const RawClusterDefs::keytype id) { PHOOL_VIRTUAL_WARN("set_id(const unsigned int id)"); }
  //  virtual void set_eta(const float eta) { PHOOL_VIRTUAL_WARN("set_eta(const float eta)");}
  virtual void set_phi(const float phi) { PHOOL_VIRTUAL_WARN("set_phi(const float phi)"); }
  virtual void set_energy(const float energy) { PHOOL_VIRTUAL_WARN("set_energy(const float energy)"); }
  virtual void set_ecore(const float ecore) { PHOOL_VIRTUAL_WARN("set_ecore(const float ecore)"); }
  virtual void set_prob(const float prob) { PHOOL_VIRTUAL_WARN("set_prob(const float prob)"); }
  virtual void set_chi2(const float chi2) { PHOOL_VIRTUAL_WARN("set_chi2(const float chi2)"); }
  virtual void addTower(const RawClusterDefs::keytype twrid, const float etower) { PHOOL_VIRTUAL_WARNING; }
  virtual size_t getNTowers() const
  {
    PHOOL_VIRTUAL_WARNING;
    return 0;
  }
  virtual TowerConstRange get_towers()
  {
    PHOOL_VIRTUAL_WARN("get_towers()");
    static TowerMap dummy;
    return make_pair(dummy.begin(), dummy.end());
  }
  virtual const TowerMap& get_towermap() const
  {
    PHOOL_VIRTUAL_WARN("get_towers()");
    static TowerMap dummy;
    return dummy;
  }

  /** @defgroup property_map property map definitions
   *  @{
   */

 public:
  //! Procedure to add a new PROPERTY tag:
  //! 1.add new tag below with unique value,
  //! 2.add a short name to RawCluster::get_property_info()
  enum PROPERTY
  {  //
    // ----- additional cluster properties -----
    //! cluster core energy for EM shower
    prop_ecore = 0,
    //! cluster template probability for EM shower
    prop_prob,
    //! reduced chi2 for EM shower
    prop_chi2,

    // ----- analysis specific quantities -----
    //! isolation ET
    prop_et_iso = 20,

    // ----- truth cluster quantities -----
    //! truth cluster's PHG4Particle ID
    prop_truth_track_ID = 100,
    //! truth cluster's PHG4Particle flavor
    prop_truth_flavor,

    //! max limit in order to fit into 8 bit unsigned number
    prop_MAX_NUMBER = UCHAR_MAX
  };

  enum PROPERTY_TYPE
  {  //
    type_int = 1,
    type_uint = 2,
    type_float = 3,
    type_unknown = -1
  };

  //! getters
  virtual bool has_property(const PROPERTY prop_id) const { return false; }
  virtual float get_property_float(const PROPERTY prop_id) const { return NAN; }
  virtual int get_property_int(const PROPERTY prop_id) const { return INT_MIN; }
  virtual unsigned int get_property_uint(const PROPERTY prop_id) const { return UINT_MAX; }
  //! setters
  virtual void set_property(const PROPERTY prop_id, const float value) { return; }
  virtual void set_property(const PROPERTY prop_id, const int value) { return; }
  virtual void set_property(const PROPERTY prop_id, const unsigned int value) { return; }
  //! type management
  static std::pair<const std::string, PROPERTY_TYPE> get_property_info(PROPERTY prop_id);
  static bool check_property(const PROPERTY prop_id, const PROPERTY_TYPE prop_type);
  static std::string get_property_type(const PROPERTY_TYPE prop_type);

  /** @} */  // end of property map definitions

 protected:
  RawCluster() {}  // make sure nobody calls ctor of base class
  ClassDef(RawCluster, 1)
};

#endif /*RAWCLUSTER_H__ */
