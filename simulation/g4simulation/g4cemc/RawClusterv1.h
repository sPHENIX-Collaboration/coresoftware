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
  int isValid() const { return _towers.size() > 0; }
  void identify(std::ostream& os=std::cout) const {
    os << "This is the RawClusterv1 object" << std::endl;
  }

  float get_eta() const { return _eta; }
  float get_phi() const { return _phi; }
  float get_energy() const { return _energy; }

  void set_eta(const float eta) { _eta = eta; }
  void set_phi(const float phi) { _phi = phi; }
  void set_energy(const float energy) { _energy = energy; }

  void addTower(const int ieta, const int iphi)
  {
    _towers.push_back(std::make_pair(ieta,iphi));
  }
  size_t getNTowers() const { return _towers.size();}
  std::pair<int,int> getTowerBin(const unsigned int itower) const;

 private:
  float _eta;
  float _phi;
  float _energy;

  std::vector<std::pair<int,int> > _towers;

  ClassDef(RawClusterv1,1)

};

#endif /*RAWCLUSTERV1_H__ */
