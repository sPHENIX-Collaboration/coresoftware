#ifndef RAWCLUSTERV1_H__
#define RAWCLUSTERV1_H__

#include "RawCluster.h"

#include <cmath>
#include <map>
#include <vector>

#ifdef __CINT__
#include <stdint.h>
#else
#include <cstdint>
#endif

class RawClusterv1 : public RawCluster
{
 public:
  RawClusterv1();
  virtual ~RawClusterv1() {}
  virtual void Reset();
//  virtual PHObject* clone() const;
  virtual int isValid() const { return towermap.size() > 0; }
  virtual void identify(std::ostream& os = std::cout) const;

  /** @defgroup getters
   *  @{
   */
  //! cluster ID
  RawClusterDefs::keytype get_id() const { return clusterid; }
  //! total energy
  float get_energy() const { return _energy; }
  //! Tower operations
  size_t getNTowers() const { return towermap.size(); }
  RawCluster::TowerConstRange get_towers() { return make_pair(towermap.begin(), towermap.end()); }
  //! return tower map for c++11 range-based for-loop
  const TowerMap& get_towermap() const { return towermap; }
  //
  //! cluster position in 3D
  virtual CLHEP::Hep3Vector get_position() const
  {
    return CLHEP::Hep3Vector(get_x(),get_y(),get_z());
  }
  //!  access to intrinsic cylindrical coordinate system
  float get_phi() const { return _phi; }
  float get_r() const { return _r; }
  float get_z() const { return _z; }
  //
//  //! convert cluster location to psuedo-rapidity given a user chosen z-location
//  virtual float get_eta(const float z) const;
//  //! convert cluster E_T given a user chosen z-location
//  virtual float get_et(const float z) const;
  //
  //! access Cartesian coordinate system
  virtual float get_x() const { return get_r() * std::cos(get_phi()); }
  virtual float get_y() const { return get_r() * std::sin(get_phi()); }
  //
  //! access additional optional properties
  //! cluster core energy for EM shower
  virtual float get_ecore() const { return get_property_float(prop_ecore); }
  //! reduced chi2 for EM shower
  virtual float get_chi2() const { return get_property_float(prop_chi2); }
  //! cluster template probability for EM shower
  virtual float get_prob() const { return get_property_float(prop_prob); }
  //! isolation ET
  virtual float get_et_iso() const { return get_property_float(prop_et_iso); }
//  //! truth cluster's PHG4Particle ID
//  virtual int get_truth_track_ID() const { return get_property_int(prop_truth_track_ID); }
//  //! truth cluster's PHG4Particle flavor
//  virtual int get_truth_flavor() const { return get_property_int(prop_truth_flavor); }
  //
  /** @} */  // end of getters

  /** @defgroup setters
   *  @{
   */
  //! cluster ID
  void set_id(const RawClusterDefs::keytype id) { clusterid = id; }
  //! Tower operations
  void addTower(const RawClusterDefs::keytype twrid, const float etower);
  //! total energy
  void set_energy(const float energy) { _energy = energy; }
  //!  access to intrinsic cylindrical coordinate system
  void set_phi(const float phi) { _phi = phi; }
  void set_z(const float z) { _z = z; }
  void set_r(const float r) { _r = r; }
  //
  //! access additional optional properties
  //! cluster core energy for EM shower
  virtual void set_ecore(const float ecore) { set_property(prop_ecore, ecore); }
  //! reduced chi2 for EM shower
  virtual void set_chi2(const float chi2) { set_property(prop_chi2, chi2); }
  //! cluster template probability for EM shower
  virtual void set_prob(const float prob) { set_property(prop_prob, prob); }
  //! isolation ET
  virtual void set_et_iso(const float e) { set_property(prop_et_iso, e); }
//  //! truth cluster's PHG4Particle ID
//  virtual void set_truth_track_ID(const int i) { set_property(prop_truth_track_ID, i); }
//  //! truth cluster's PHG4Particle flavor
//  virtual void set_truth_flavor(const int f) { set_property(prop_truth_flavor, f); }
  //
  /*
   *
   * @} */  // end of setters

  /** @defgroup property_map property map definitions
   *  @{
   */
 public:
  bool has_property(const PROPERTY prop_id) const;
  float get_property_float(const PROPERTY prop_id) const;
  int get_property_int(const PROPERTY prop_id) const;
  unsigned int get_property_uint(const PROPERTY prop_id) const;
  void set_property(const PROPERTY prop_id, const float value);
  void set_property(const PROPERTY prop_id, const int value);
  void set_property(const PROPERTY prop_id, const unsigned int value);

 protected:
  unsigned int get_property_nocheck(const PROPERTY prop_id) const;
  void set_property_nocheck(const PROPERTY prop_id, const unsigned int ui) { prop_map[prop_id] = ui; }
  //! storage types for additional property
  typedef uint8_t prop_id_t;
  typedef uint32_t prop_storage_t;
  typedef std::map<prop_id_t, prop_storage_t> prop_map_t;

  //! convert between 32bit inputs and storage type prop_storage_t
  union u_property {
    float fdata;
    int32_t idata;
    uint32_t uidata;

    u_property(int32_t in)
      : idata(in)
    {
    }
    u_property(uint32_t in)
      : uidata(in)
    {
    }
    u_property(float in)
      : fdata(in)
    {
    }
    u_property()
      : uidata(0)
    {
    }

    prop_storage_t get_storage() const { return uidata; }
  };

  //! container for additional property
  prop_map_t prop_map;

  /** @} */  // end of property map definitions

  //
 protected:
  //! cluster ID
  RawClusterDefs::keytype clusterid;
  //! total energy
  float _energy;
  //! Tower operations
  TowerMap towermap;

  //! location of cluster in cylindrical coordinate
  float _r;
  float _phi;
  float _z;

  ClassDef(RawClusterv1, 3)
};

#endif /*RAWCLUSTERV1_H__ */
