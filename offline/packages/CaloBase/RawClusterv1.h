#ifndef CALOBASE_RAWCLUSTERV1_H
#define CALOBASE_RAWCLUSTERV1_H

#include "RawCluster.h"
#include "RawClusterDefs.h"

#include <CLHEP/Vector/ThreeVector.h>

#include <cstddef>
#include <cstdint>
#include <iostream>
#include <map>
#include <utility>

class PHObject;

class RawClusterv1 : public RawCluster
{
 public:
  RawClusterv1();
  ~RawClusterv1() override {}

  void Reset() override;
  PHObject* CloneMe() const override { return new RawClusterv1(*this); }
  int isValid() const override { return towermap.size() > 0; }
  void identify(std::ostream& os = std::cout) const override;

  /** @defgroup getters
   *  @{
   */
  //! cluster ID
  RawClusterDefs::keytype get_id() const override { return clusterid; }
  //! total energy
  float get_energy() const override { return _energy; }
  //! Tower operations
  size_t getNTowers() const override { return towermap.size(); }
  RawCluster::TowerConstRange get_towers() const override { return make_pair(towermap.begin(), towermap.end()); }
  //! return tower map for c++11 range-based for-loop
  const TowerMap& get_towermap() const override { return towermap; }
  //
  //! cluster position in 3D
  CLHEP::Hep3Vector get_position() const override
  {
    return CLHEP::Hep3Vector(get_x(), get_y(), get_z());
  }
  //!  access to intrinsic cylindrical coordinate system
  float get_phi() const override { return _phi; }
  float get_r() const override { return _r; }
  float get_z() const override { return _z; }
  //
  //  //! convert cluster location to psuedo-rapidity given a user chosen z-location
  //  virtual float get_eta(const float z) const;
  //  //! convert cluster E_T given a user chosen z-location
  //  virtual float get_et(const float z) const;
  //
  //! access Cartesian coordinate system
  float get_x() const override { return get_r() * std::cos(get_phi()); }
  float get_y() const override { return get_r() * std::sin(get_phi()); }
  //
  //! access additional optional properties
  //! cluster core energy for EM shower
  float get_ecore() const override { return get_property_float(prop_ecore); }
  //! reduced chi2 for EM shower
  float get_chi2() const override { return get_property_float(prop_chi2); }
  //! cluster template probability for EM shower
  float get_prob() const override { return get_property_float(prop_prob); }
  //! isolation ET default
  float get_et_iso() const override { return get_property_float(prop_et_iso_calotower_R03); }
  //! isolation ET the radius and hueristic can be specified
  float get_et_iso(const int radiusx10, bool subtracted, bool clusterTower) const override;
  //  //! truth cluster's PHG4Particle ID
  //  virtual int get_truth_track_ID() const override { return get_property_int(prop_truth_track_ID); }
  //  //! truth cluster's PHG4Particle flavor
  //  virtual int get_truth_flavor() const override { return get_property_int(prop_truth_flavor); }
  //
  /** @} */  // end of getters

  /** @defgroup setters
   *  @{
   */
  //! cluster ID
  void set_id(const RawClusterDefs::keytype id) override { clusterid = id; }
  //! Tower operations
  void addTower(const RawClusterDefs::keytype twrid, const float etower) override;
  //! total energy
  void set_energy(const float energy) override { _energy = energy; }
  //!  access to intrinsic cylindrical coordinate system
  void set_phi(const float phi) override { _phi = phi; }
  void set_z(const float z) override { _z = z; }
  void set_r(const float r) override { _r = r; }
  //
  //! access additional optional properties
  //! cluster core energy for EM shower
  void set_ecore(const float ecore) override { set_property(prop_ecore, ecore); }
  //! reduced chi2 for EM shower
  void set_chi2(const float chi2) override { set_property(prop_chi2, chi2); }
  //! cluster template probability for EM shower
  void set_prob(const float prob) override { set_property(prop_prob, prob); }
  //! isolation ET default
  void set_et_iso(const float e) override { set_property(prop_et_iso_calotower_R03, e); }
  //! isolation ET the radius and hueristic can be specified
  void set_et_iso(const float et_iso, const int radiusx10, bool subtracted, bool clusterTower) override;
  //  //! truth cluster's PHG4Particle ID
  //  virtual void set_truth_track_ID(const int i) override { set_property(prop_truth_track_ID, i); }
  //  //! truth cluster's PHG4Particle flavor
  //  virtual void set_truth_flavor(const int f) override { set_property(prop_truth_flavor, f); }
  //
  /*
   *
   * @} */
  // end of setters

  /** @defgroup property_map property map definitions
   *  @{
   */
 public:
  bool has_property(const PROPERTY prop_id) const override;
  float get_property_float(const PROPERTY prop_id) const override;
  int get_property_int(const PROPERTY prop_id) const override;
  unsigned int get_property_uint(const PROPERTY prop_id) const override;
  void set_property(const PROPERTY prop_id, const float value) override;
  void set_property(const PROPERTY prop_id, const int value) override;
  void set_property(const PROPERTY prop_id, const unsigned int value) override;

 protected:  // protected is declared twice !?
  unsigned int get_property_nocheck(const PROPERTY prop_id) const;
  void set_property_nocheck(const PROPERTY prop_id, const unsigned int ui) { prop_map[prop_id] = ui; }
  //! storage types for additional property
  typedef uint8_t prop_id_t;
  typedef uint32_t prop_storage_t;
  typedef std::map<prop_id_t, prop_storage_t> prop_map_t;

  //! convert between 32bit inputs and storage type prop_storage_t
  union u_property
  {
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

  ClassDefOverride(RawClusterv1, 3)
};

#endif
