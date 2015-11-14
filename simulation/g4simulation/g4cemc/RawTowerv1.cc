#include "RawTowerv1.h"

#include "RawTowerDefs.h"

#include <iostream>
#include <algorithm>

#include <cmath>
#include <map>

using namespace std;

ClassImp(RawTowerv1)

RawTowerv1::RawTowerv1() : 
  towerid(~0), // initialize all bits on
  energy(0)
{}

RawTowerv1::RawTowerv1(const RawTower & tower)
{
  towerid = (tower.get_id());
  energy = (tower.get_energy());

  CellConstRange cell_range = tower.get_g4cells();

  for (CellConstIterator cell_iter = cell_range.first;
      cell_iter != cell_range.second; ++cell_iter)
    {
      add_ecell(cell_iter->first, cell_iter->second);
    }
}

RawTowerv1::RawTowerv1(RawTowerDefs::keytype id) : 
  towerid(id),
  energy(0)
{}

RawTowerv1::RawTowerv1(const unsigned int ieta, const unsigned int iphi) :
  towerid(0),
  energy(0)
{
  towerid = RawTowerDefs::encode_towerid( RawTowerDefs::NONE , ieta , iphi );
}

RawTowerv1::RawTowerv1(const RawTowerDefs::CalorimeterId caloid, const unsigned int ieta, const unsigned int iphi) :
  towerid(0),
  energy(0)
{
  towerid = RawTowerDefs::encode_towerid( caloid , ieta , iphi );
}

RawTowerv1::~RawTowerv1()
{}

void RawTowerv1::Reset()
{
  ecells.clear();
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

void 
RawTowerv1::add_ecell(const PHG4CylinderCellDefs::keytype g4cellid, const float ecell)
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


