#include "RawTowerv2.h"

#include "RawTowerDefs.h"

#include <iostream>
#include <algorithm>

#include <cmath>
#include <map>

using namespace std;

ClassImp(RawTowerv2)

RawTowerv2::RawTowerv2() :
    towerid(~0), // initialize all bits on
  energy(0), time(NAN), _tower_type(-1)
{
}

RawTowerv2::RawTowerv2(const RawTower & tower)
{
  towerid = (tower.get_id());
  energy = (tower.get_energy());
  time = (tower.get_time());
  _tower_type = (tower.get_tower_type()); 

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

RawTowerv2::RawTowerv2(RawTowerDefs::keytype id) :
    towerid(id), energy(0), time(NAN)
{
}

RawTowerv2::RawTowerv2(const unsigned int ieta, const unsigned int iphi) :
    towerid(0), energy(0)
{
  towerid = RawTowerDefs::encode_towerid(RawTowerDefs::NONE, ieta, iphi);
}

RawTowerv2::RawTowerv2(const RawTowerDefs::CalorimeterId caloid,
    const unsigned int ieta, const unsigned int iphi) :
    towerid(0), energy(0), time(NAN)
{
  towerid = RawTowerDefs::encode_towerid(caloid, ieta, iphi);
}

RawTowerv2::~RawTowerv2()
{
}

void
RawTowerv2::Reset()
{
  energy = 0;
  time = NAN;
  ecells.clear();
  eshowers.clear();
}

int
RawTowerv2::isValid() const
{
  return get_energy() != 0;
}

void
RawTowerv2::identify(std::ostream& os) const
{
  os << "RawTowerv2: etabin: " << get_bineta() << ", phibin: " << get_binphi()
      << " energy=" << get_energy() << std::endl;
}

void
RawTowerv2::add_ecell(const PHG4CylinderCellDefs::keytype g4cellid,
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

void 
RawTowerv2::add_eshower(const int g4showerid, const float eshower)
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
