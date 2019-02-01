#ifndef CALOBASE_RAWTOWERV1_H
#define CALOBASE_RAWTOWERV1_H

#include "RawTower.h"

#include "RawTowerDefs.h"

#include <map>

class RawTowerv1 : public RawTower
{
 public:
  RawTowerv1();
  RawTowerv1(const RawTower& tower);
  RawTowerv1(RawTowerDefs::keytype id);
  RawTowerv1(const unsigned int ieta, const unsigned int iphi);
  RawTowerv1(const RawTowerDefs::CalorimeterId caloid, const unsigned int ieta,
             const unsigned int iphi);
  virtual ~RawTowerv1() {}

  void Reset();
  int isValid() const;
  void identify(std::ostream& os = std::cout) const;

  void set_id(RawTowerDefs::keytype id) { towerid = id; }
  RawTowerDefs::keytype get_id() const { return towerid; }
  int get_bineta() const { return RawTowerDefs::decode_index1(towerid); }
  int get_binphi() const { return RawTowerDefs::decode_index2(towerid); }
  double get_energy() const { return energy; }
  void set_energy(const double e) { energy = e; }
  float get_time() const { return time; }
  void set_time(const float t) { time = t; }

  //---cell access--------------------------------------------------------------

  bool empty_g4cells() const { return ecells.empty(); }
  size_t size_g4cells() const { return ecells.size(); }
  RawTower::CellConstRange get_g4cells() const
  {
    return make_pair(ecells.begin(), ecells.end());
  }
  RawTower::CellIterator find_g4cell(int id) { return ecells.find(id); }
  RawTower::CellConstIterator find_g4cell(int id) const { return ecells.find(id); }
  void add_ecell(const CellKeyType g4cellid,
                 const float ecell);
  void clear_g4cells() { ecells.clear(); }

  //---shower access------------------------------------------------------------

  bool empty_g4showers() const { return eshowers.empty(); }
  size_t size_g4showers() const { return eshowers.size(); }
  RawTower::ShowerConstRange get_g4showers() const
  {
    return make_pair(eshowers.begin(), eshowers.end());
  }
  RawTower::ShowerIterator find_g4shower(int id) { return eshowers.find(id); }
  RawTower::ShowerConstIterator find_g4shower(int id) const { return eshowers.find(id); }
  void add_eshower(const int g4showerid, const float eshower);
  void clear_g4showers() { eshowers.clear(); }

 protected:
  RawTowerDefs::keytype towerid;

  //! energy assigned to the tower. Depending on stage of process and DST node
  //! name, it could be energy deposition, light yield or calibrated energies
  double energy;
  //! Time stamp assigned to the tower. Depending on the tower maker, it could
  //! be rise time or peak time.
  float time;

  CellMap ecells;      //< default truth storage
  ShowerMap eshowers;  //< alternate truth storage for smaller filesizes

  ClassDef(RawTowerv1, 5)
};

#endif
