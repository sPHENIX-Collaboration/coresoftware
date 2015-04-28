#ifndef RAWTOWERV1_H_
#define RAWTOWERV1_H_

#include "RawTower.h"

#include <map>

class RawTowerv1 : public RawTower {

 public:
  RawTowerv1();
  RawTowerv1(const int ieta, const int iphi);
  virtual ~RawTowerv1();

  void Reset();
  int isValid() const;
  void identify(std::ostream& os=std::cout) const;

  int get_bineta() const { return bineta; }
  int get_binphi() const { return binphi; }
  float get_energy() const;

  bool is_adjacent(RawTower& tower);

  std::pair< std::map<unsigned int,float>::const_iterator, std::map<unsigned int,float>::const_iterator > get_g4cells()
  {return make_pair(ecells.begin(), ecells.end());}
  void add_ecell(const unsigned int g4cellid, const float ecell);

 protected:
  int bineta;
  int binphi;

  std::map<unsigned int, float> ecells;

  ClassDef(RawTowerv1,1)
};

#endif /* RAWTOWERV1_H_ */
