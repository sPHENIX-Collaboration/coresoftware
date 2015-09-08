#ifndef RAWTOWERV1_H_
#define RAWTOWERV1_H_

#include "RawTower.h"

#include "RawTowerDefs.h"

#include <map>

class RawTowerv1 : public RawTower {

 public:
  RawTowerv1();
  RawTowerv1(RawTowerDefs::keytype id);
  RawTowerv1(const unsigned int ieta, const unsigned int iphi);
  virtual ~RawTowerv1();

  void Reset();
  int isValid() const;
  void identify(std::ostream& os=std::cout) const;

  RawTowerDefs::keytype get_id() const { return towerid;}
  int get_bineta() const { return (towerid >> RawTowerDefs::eta_idbits)&0xFFF ; }
  int get_binphi() const { return towerid&0xFFF; }
  double get_energy() const;

  void set_light_yield(const float l)  { light_yield = l; }
  float get_light_yield() const { return light_yield; };

  RawTower::CellConstRange get_g4cells()
  {return make_pair(ecells.begin(), ecells.end());}
  void add_ecell(const PHG4CylinderCellDefs::keytype g4cellid, const float ecell);

 protected:
  RawTowerDefs::keytype towerid;
  float light_yield;

  CellMap ecells;

  ClassDef(RawTowerv1,2)
};

#endif /* RAWTOWERV1_H_ */
