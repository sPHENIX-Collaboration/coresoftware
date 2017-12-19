#ifndef RAWCLUSTERV1_H__
#define RAWCLUSTERV1_H__

#include <cmath>
#include <map>
#include <vector>
#include "RawCluster.h"

class RawClusterv1 : public RawCluster
{
 public:
  RawClusterv1();
  virtual ~RawClusterv1() {}
  void Reset();
  int isValid() const { return towermap.size() > 0; }
  void identify(std::ostream& os = std::cout) const;

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
  //!  access to intrinsic cylindrical coordinate system
  float get_phi() const { return _phi; }
  float get_r() const { return _r; }
  float get_z() const { return _z; }
  //! convert cluster location to psuedo-rapidity given a user chosen z-location
  virtual float get_eta(const float z) const;
  //! convert cluster E_T given a user chosen z-location
  virtual float get_et(const float z) const;
  //
  //! access Cartesian coordinate system
  virtual float get_x() const { return get_r() * std::cos(get_phi()); }
  virtual float get_y() const { return get_r() * std::sin(get_phi()); }
  //
  //! access additional optional properties
  float get_ecore() const { return _ecore; }
  float get_chi2() const { return _chi2; }
  float get_prob() const { return _prob; }
  /** @} */ // end of getters

  /** @defgroup setters
   *  @{
   */
  //! cluster ID
  void set_id(const RawClusterDefs::keytype id) { clusterid = id; }
  //! Tower operations
  void addTower(const RawClusterDefs::keytype twrid, const float etower);
  //! total energy
  void set_energy(const float energy) { _energy = energy; }
  //
  //!  access to intrinsic cylindrical coordinate system
  void set_phi(const float phi) { _phi = phi; }
  void set_z(const float z) { _z = z; }
  void set_r(const float r) { _r = r; }
  //
  //! access additional optional properties
  void set_ecore(const float ecore) { _ecore = ecore; }
  void set_chi2(const float chi2) { _chi2 = chi2; }
  void set_prob(const float prob) { _prob = prob; }


  /** @} */ // end of setters

  //
 private:
  //! cluster ID
  RawClusterDefs::keytype clusterid;
  //! total energy
  float _energy;
  //! Tower operations
  TowerMap towermap;

  //! location of cluster in cylindrical coordinate
  float _z;
  float _phi;
  float _r;


  /** @defgroup property_map property map definitions
   *  @{
   */


  /** @} */ // end of property map definitions


  ClassDef(RawClusterv1, 3)
};

#endif /*RAWCLUSTERV1_H__ */
