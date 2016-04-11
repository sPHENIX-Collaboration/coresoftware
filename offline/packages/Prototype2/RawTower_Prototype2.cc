#include "RawTower_Prototype2.h"
#include <g4cemc/RawTowerDefs.h>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <map>

using namespace std;

ClassImp(RawTower_Prototype2)

RawTower_Prototype2::RawTower_Prototype2() :
    towerid(~0), // initialize all bits on
    energy(0), time(NAN)
{
}

RawTower_Prototype2::RawTower_Prototype2(const RawTower & tower)
{
  towerid = (tower.get_id());
  energy = (tower.get_energy());
  time = (tower.get_time());

  CellConstRange cell_range = tower.get_g4cells();

  for (CellConstIterator cell_iter = cell_range.first;
      cell_iter != cell_range.second; ++cell_iter)
    {
      add_ecell(cell_iter->first, cell_iter->second);
    }

  ShowerConstRange shower_range = tower.get_g4showers();

  for (ShowerConstIterator shower_iter = shower_range.first;
      shower_iter != shower_range.second; ++shower_iter)
    {
      add_eshower(shower_iter->first, shower_iter->second);
    }
}

RawTower_Prototype2::RawTower_Prototype2(RawTowerDefs::keytype id) :
    towerid(id), energy(0), time(NAN)
{
}

RawTower_Prototype2::RawTower_Prototype2(const unsigned int ieta, const unsigned int iphi) :
    towerid(0), energy(0)
{
  towerid = RawTowerDefs::encode_towerid(RawTowerDefs::NONE, ieta, iphi);
}

RawTower_Prototype2::RawTower_Prototype2(const RawTowerDefs::CalorimeterId caloid,
    const unsigned int ieta, const unsigned int iphi) :
    towerid(0), energy(0), time(NAN)
{
  towerid = RawTowerDefs::encode_towerid(caloid, ieta, iphi);
}

RawTower_Prototype2::~RawTower_Prototype2()
{
}

void
RawTower_Prototype2::Reset()
{
  energy = 0;
  time = NAN;
}

int
RawTower_Prototype2::isValid() const
{
  return get_energy() != 0;
}

void
RawTower_Prototype2::identify(std::ostream& os) const
{
  os << "RawTower_Prototype2: etabin: " << get_bineta() << ", phibin: " << get_binphi()
      << " energy=" << get_energy() << std::endl;
}

