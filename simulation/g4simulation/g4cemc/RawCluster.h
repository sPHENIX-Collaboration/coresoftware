#ifndef RAWCLUSTER_H__
#define RAWCLUSTER_H__

#include "RawClusterDefs.h"
#include "RawTowerDefs.h"

#include <phool/PHObject.h>
#include <phool/phool.h>
#include <cmath> // def of NAN
#include <vector>
#include <map>

class RawTower;

class RawCluster : public PHObject {

 public:

  typedef std::map<RawTowerDefs::keytype, float> TowerMap;
  typedef TowerMap::iterator TowerIterator;
  typedef TowerMap::const_iterator TowerConstIterator;
  typedef std::pair<TowerIterator, TowerIterator> TowerRange;
  typedef std::pair<TowerConstIterator, TowerConstIterator> TowerConstRange;


  virtual ~RawCluster() {}

  virtual void Reset() { PHOOL_VIRTUAL_WARNING; }
  virtual int isValid() const { PHOOL_VIRTUAL_WARNING; return 0; }
  virtual void identify(std::ostream& os=std::cout) const { PHOOL_VIRTUAL_WARNING; }

  virtual RawClusterDefs::keytype get_id() const { PHOOL_VIRTUAL_WARN("get_id()"); return 0; }
  virtual float get_eta() const { PHOOL_VIRTUAL_WARN("get_eta()"); return NAN; }
  virtual float get_phi() const { PHOOL_VIRTUAL_WARN("get_phi()"); return NAN; }
  virtual float get_energy() const { PHOOL_VIRTUAL_WARN("get_energy()"); return NAN; }
  virtual float get_prob() const { PHOOL_VIRTUAL_WARN("get_prob()"); return NAN; }
  virtual float get_chi2() const { PHOOL_VIRTUAL_WARN("get_chi2()"); return NAN; }

  virtual void set_id(const RawClusterDefs::keytype id) { PHOOL_VIRTUAL_WARN("set_id(const unsigned int id)");}
  virtual void set_eta(const float eta) { PHOOL_VIRTUAL_WARN("set_eta(const float eta)");}
  virtual void set_phi(const float phi) { PHOOL_VIRTUAL_WARN("set_phi(const float phi)");}
  virtual void set_energy(const float energy) { PHOOL_VIRTUAL_WARN("set_energy(const float energy)");}
  virtual void set_prob(const float prob) { PHOOL_VIRTUAL_WARN("set_prob(const float prob)");}
  virtual void set_chi2(const float chi2) { PHOOL_VIRTUAL_WARN("set_chi2(const float chi2)");}

  virtual void addTower(const RawClusterDefs::keytype twrid, const float etower) { PHOOL_VIRTUAL_WARNING; }
  virtual size_t getNTowers() const { PHOOL_VIRTUAL_WARNING; return 0;}
  virtual TowerConstRange get_towers()
  {
    PHOOL_VIRTUAL_WARN("get_towers()");
    TowerMap dummy;
    return make_pair(dummy.begin(), dummy.end());
  }

 protected:
  RawCluster() {} // make sure nobody calls ctor of base class

  ClassDef(RawCluster,1)

};

#endif /*RAWCLUSTER_H__ */
