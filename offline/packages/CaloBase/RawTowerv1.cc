#include "RawTowerv1.h"

#include <cmath>
#include <iostream>

using namespace std;

RawTowerv1::RawTowerv1(const RawTower& tower)
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

RawTowerv1::RawTowerv1(RawTowerDefs::keytype id)
  : towerid(id)
{
}

RawTowerv1::RawTowerv1(const unsigned int ieta, const unsigned int iphi)
{
  towerid = RawTowerDefs::encode_towerid(RawTowerDefs::NONE, ieta, iphi);
}

RawTowerv1::RawTowerv1(const RawTowerDefs::CalorimeterId caloid,
                       const unsigned int ieta, const unsigned int iphi)
{
  towerid = RawTowerDefs::encode_towerid(caloid, ieta, iphi);
}

void RawTowerv1::Reset()
{
  energy = 0;
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
  RawTowerDefs::CalorimeterId caloid = RawTowerDefs::decode_caloid(towerid);
  switch (caloid)
  {
  case RawTowerDefs::LFHCAL:
    os << "RawTowerv1: etabin: " << get_bineta() << ", phibin: " << get_binphi() << ", l-bin: " << get_binl()
       << " energy=" << get_energy() << std::endl;
    return;
  default:
    os << "RawTowerv1: etabin: " << get_bineta() << ", phibin: " << get_binphi()
       << " energy=" << get_energy() << std::endl;
    return;
  }
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

int RawTowerv1::get_bineta() const
{
  RawTowerDefs::CalorimeterId caloid = RawTowerDefs::decode_caloid(towerid);
  switch (caloid)
  {
  case RawTowerDefs::LFHCAL:
    return RawTowerDefs::decode_index1v2(towerid);
  default:
    return RawTowerDefs::decode_index1(towerid);
  }
  return -1;
}

int RawTowerv1::get_binphi() const
{
  RawTowerDefs::CalorimeterId caloid = RawTowerDefs::decode_caloid(towerid);
  switch (caloid)
  {
  case RawTowerDefs::LFHCAL:
    return RawTowerDefs::decode_index2v2(towerid);
  default:
    return RawTowerDefs::decode_index2(towerid);
  }
  return -1;
}
