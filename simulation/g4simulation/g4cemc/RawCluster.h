#ifndef RAWCLUSTER_H__
#define RAWCLUSTER_H__

#include <phool/PHObject.h>
#include <phool/phool.h>
#include <cmath> // def of NAN
#include <vector>
#include <map>

class RawTower;

class RawCluster : public PHObject {

 public:
  RawCluster() {}
  virtual ~RawCluster() {}

  virtual void Reset() { PHOOL_VIRTUAL_WARNING; }
  virtual int isValid() const { PHOOL_VIRTUAL_WARNING; return 0; }
  virtual void identify(std::ostream& os=std::cout) const { PHOOL_VIRTUAL_WARNING; }

  virtual float get_eta() const { PHOOL_VIRTUAL_WARN("get_eta()"); return NAN; }
  virtual float get_phi() const { PHOOL_VIRTUAL_WARN("get_phi()"); return NAN; }
  virtual float get_energy() const { PHOOL_VIRTUAL_WARN("get_energy()"); return NAN; }

  virtual void set_eta(const float eta) { PHOOL_VIRTUAL_WARN("set_eta(const float eta)");}
  virtual void set_phi(const float phi) { PHOOL_VIRTUAL_WARN("set_phi(const float phi)");}
  virtual void set_energy(const float energy) { PHOOL_VIRTUAL_WARN("set_energy(const float energy)");}

  virtual void addTower(const int ieta, const int iphi) { PHOOL_VIRTUAL_WARNING; }
  virtual size_t getNTowers() const { PHOOL_VIRTUAL_WARNING; return 0;}
  virtual std::pair<int,int> getTowerBin(const unsigned int itower) const { PHOOL_VIRTUAL_WARNING; return std::make_pair(0,0);}

  ClassDef(RawCluster,1)

};

#endif /*RAWCLUSTER_H__ */
