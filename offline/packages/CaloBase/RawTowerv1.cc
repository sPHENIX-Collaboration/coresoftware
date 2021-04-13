#include "RawTowerv1.h"

#include <cmath>
#include <iostream>

using namespace std;

RawTowerv1::RawTowerv1()
  : towerid(~0)
  ,  // initialize all bits on
  energy(0)
  , scint_gammas(0)
  , cerenkov_gammas(0)
  , time(NAN)
{
}

RawTowerv1::RawTowerv1(const RawTower& tower)
{
  towerid = (tower.get_id());
  energy = (tower.get_energy());
  scint_gammas = (tower.get_scint_gammas());
  cerenkov_gammas = (tower.get_cerenkov_gammas());
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

RawTowerv1::RawTowerv1(RawTowerDefs::keytype id)
  : towerid(id)
  , energy(0)
  , scint_gammas(0)
  , cerenkov_gammas(0)
  , time(NAN)
{
}

RawTowerv1::RawTowerv1(const unsigned int ieta, const unsigned int iphi)
  : towerid(0)
  , energy(0)
  , scint_gammas(0)
  , cerenkov_gammas(0)
  , time(NAN)
{
  towerid = RawTowerDefs::encode_towerid(RawTowerDefs::NONE, ieta, iphi);
}

RawTowerv1::RawTowerv1(const RawTowerDefs::CalorimeterId caloid,
                       const unsigned int ieta, const unsigned int iphi)
  : towerid(0)
  , energy(0)
  , scint_gammas(0)
  , cerenkov_gammas(0)
  , time(NAN)
{
  towerid = RawTowerDefs::encode_towerid(caloid, ieta, iphi);
}

void RawTowerv1::Reset()
{
  energy = 0;
  scint_gammas = 0;
  cerenkov_gammas = 0;
  time = NAN;
  ecells.clear();
  eshowers.clear();
}

int RawTowerv1::isValid() const
{
  return get_energy() != 0;
}

void RawTowerv1::identify(std::ostream& os) const
{
  os << "RawTowerv1: etabin: " << get_bineta() << ", phibin: " << get_binphi()
     << " energy=" << get_energy() << std::endl;
}

void RawTowerv1::add_ecell(const CellKeyType g4cellid,
                           const float ecell)
{
  if (ecells.find(g4cellid) == ecells.end())
  {
    ecells[g4cellid] = ecell;
  }
  else
  {
    ecells[g4cellid] += ecell;
  }
}

void RawTowerv1::add_eshower(const int g4showerid, const float eshower)
{
  if (eshowers.find(g4showerid) == eshowers.end())
  {
    eshowers[g4showerid] = eshower;
  }
  else
  {
    eshowers[g4showerid] += eshower;
  }
}
