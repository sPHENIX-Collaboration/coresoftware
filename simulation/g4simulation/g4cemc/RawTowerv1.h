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

  RawTowerDefs::keytype get_id() const { return bineta << 12 | binphi;}
  int get_bineta() const { return bineta; }
  int get_binphi() const { return binphi; }
  float get_energy() const;

  void set_light_yield(float l)  { light_yield = l; }
  float get_light_yield() const { return light_yield; };

  bool is_adjacent(RawTower& tower);

  CellConstRange get_g4cells()
  {return make_pair(ecells.begin(), ecells.end());}
  void add_ecell(const PHG4CylinderCellDefs::keytype g4cellid, const float ecell);

 protected:
  int bineta;
  int binphi;
  float light_yield;

  CellMap ecells;

  ClassDef(RawTowerv1,2)
};

#endif /* RAWTOWERV1_H_ */
