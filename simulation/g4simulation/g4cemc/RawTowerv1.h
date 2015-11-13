#ifndef RAWTOWERV1_H_
#define RAWTOWERV1_H_

#include "RawTower.h"

#include "RawTowerDefs.h"

#include <map>

class RawTowerv1 : public RawTower {

 public:
  RawTowerv1();
  RawTowerv1(const RawTower & tower);
  RawTowerv1(RawTowerDefs::keytype id);
  RawTowerv1(const unsigned int ieta, const unsigned int iphi);
  virtual ~RawTowerv1();

  void Reset();
  int isValid() const;
  void identify(std::ostream& os=std::cout) const;

  RawTowerDefs::keytype get_id() const { return towerid;}
  int get_bineta() const { return RawTowerDefs::decode_index1(towerid); }
  int get_binphi() const { return RawTowerDefs::decode_index2(towerid); }
  double get_energy() const {return energy;}
  void set_energy(const double e) {energy = e;}

  RawTower::CellConstRange get_g4cells() const
  {return make_pair(ecells.begin(), ecells.end());}
  void add_ecell(const PHG4CylinderCellDefs::keytype g4cellid, const float ecell);

 protected:
  RawTowerDefs::keytype towerid;

  //! energy assigned to the tower. Depending on stage of process and DST node name, it could be energy deposition, light yield or calibrated energies
  double energy;

  CellMap ecells;

  ClassDef(RawTowerv1,3)
};

#endif /* RAWTOWERV1_H_ */
