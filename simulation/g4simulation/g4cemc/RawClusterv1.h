#ifndef RAWCLUSTERV1_H__
#define RAWCLUSTERV1_H__

#include "RawCluster.h"
#include <vector>
#include <map>

class RawClusterv1 : public RawCluster {

 public:
  RawClusterv1();
  virtual ~RawClusterv1() {}

  void Reset();
  int isValid() const { return towermap.size() > 0; }
  void identify(std::ostream& os=std::cout) const {
    os << "This is the RawClusterv1 object" << std::endl;
  }

  RawClusterDefs::keytype get_id() const { return clusterid; }
  float get_eta() const { return _eta; }
  float get_phi() const { return _phi; }
  float get_energy() const { return _energy; }
  float get_chi2() const { return _chi2; }
  float get_prob() const { return _prob; }

  void set_id(const RawClusterDefs::keytype id) {clusterid = id;}
  void set_eta(const float eta) { _eta = eta; }
  void set_phi(const float phi) { _phi = phi; }
  void set_energy(const float energy) { _energy = energy; }
  void set_chi2(const float chi2) { _chi2 = chi2; }
  void set_prob(const float prob) { _prob = prob; }

  void addTower(const RawClusterDefs::keytype twrid, const float etower);
  size_t getNTowers() const { return towermap.size();}

  RawCluster::TowerConstRange get_towers()
    {
      return make_pair(towermap.begin(),towermap.end());
    }

 private:
  RawClusterDefs::keytype clusterid;
  float _eta;
  float _phi;
  float _energy;
  float _chi2;
  float _prob;
  std::vector<std::pair<int,int> > _towers;

  TowerMap towermap;

  ClassDef(RawClusterv1,1)

};

#endif /*RAWCLUSTERV1_H__ */
