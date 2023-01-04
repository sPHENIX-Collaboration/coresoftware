#ifndef CALOBASE_RAWTOWERV1_H
#define CALOBASE_RAWTOWERV1_H

#include "RawTower.h"

#include "RawTowerDefs.h"

#include <cstddef>
#include <iostream>
#include <map>
#include <utility>

class RawTowerv1 : public RawTower
{
 public:
  RawTowerv1() {}
  RawTowerv1(const RawTower& tower);
  RawTowerv1(RawTowerDefs::keytype id);
  RawTowerv1(const unsigned int ieta, const unsigned int iphi);
  RawTowerv1(const RawTowerDefs::CalorimeterId caloid, const unsigned int ieta,
             const unsigned int iphi);
  ~RawTowerv1() override {}

  void Reset() override;
  int isValid() const override;
  void identify(std::ostream& os = std::cout) const override;

  void set_id(RawTowerDefs::keytype id) override { towerid = id; }
  RawTowerDefs::keytype get_id() const override { return towerid; }
  int get_bineta() const override;
  int get_binphi() const override;
  int get_binl() const override { return RawTowerDefs::decode_index3v2(towerid); }
  double get_energy() const override { return energy; }
  void set_energy(const double e) override { energy = e; }
  float get_time() const override { return time; }
  void set_time(const float t) override { time = t; }

  //---cell access--------------------------------------------------------------

  bool empty_g4cells() const override { return ecells.empty(); }
  size_t size_g4cells() const override { return ecells.size(); }
  RawTower::CellConstRange get_g4cells() const override
  {
    return make_pair(ecells.begin(), ecells.end());
  }
  RawTower::CellIterator find_g4cell(CellKeyType id) override { return ecells.find(id); }
  RawTower::CellConstIterator find_g4cell(CellKeyType id) const override { return ecells.find(id); }
  void add_ecell(const CellKeyType g4cellid,
                 const float ecell) override;
  void clear_g4cells() override { ecells.clear(); }

  //---shower access------------------------------------------------------------

  bool empty_g4showers() const override { return eshowers.empty(); }
  size_t size_g4showers() const override { return eshowers.size(); }
  RawTower::ShowerConstRange get_g4showers() const override
  {
    return make_pair(eshowers.begin(), eshowers.end());
  }
  RawTower::ShowerIterator find_g4shower(int id) override { return eshowers.find(id); }
  RawTower::ShowerConstIterator find_g4shower(int id) const override { return eshowers.find(id); }
  void add_eshower(const int g4showerid, const float eshower) override;
  void clear_g4showers() override { eshowers.clear(); }

 protected:
  RawTowerDefs::keytype towerid = ~0;

  //! energy assigned to the tower. Depending on stage of process and DST node
  //! name, it could be energy deposition, light yield or calibrated energies
  double energy = 0.;
  //! Time stamp assigned to the tower. Depending on the tower maker, it could
  //! be rise time or peak time.
  float time = NAN;

  CellMap ecells;      //< default truth storage
  ShowerMap eshowers;  //< alternate truth storage for smaller filesizes

  ClassDefOverride(RawTowerv1, 5)
};

#endif
