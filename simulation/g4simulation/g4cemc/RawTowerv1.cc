#include "RawTowerv1.h"
#include <iostream>
#include <algorithm>

#include <cmath>
#include <map>

using namespace std;

ClassImp(RawTowerv1)

RawTowerv1::RawTowerv1() : 
  towerid(~0), // initialize all bits on
  light_yield(NAN)
{}

RawTowerv1::RawTowerv1(RawTowerDefs::keytype id) : 
  towerid(id),
  light_yield(NAN)
{}

RawTowerv1::RawTowerv1(const unsigned int ieta, const unsigned int iphi) :
  towerid(0),
  light_yield(NAN)
{
  if (ieta < 0xFFF && iphi < 0xFFF)
    {
  towerid = (ieta << 12) + iphi;
    }
  else
    {
      cout << "too large eta or phi bin, eta: " << ieta
	   << ", phi: " << iphi << ", max val: " << 0xFFF << endl;
      exit(1);
    }
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

double 
RawTowerv1::get_energy() const
{
  RawTower::CellConstIterator iter;
  double esum = 0;
  for (iter = ecells.begin(); iter != ecells.end(); ++iter)
    {
      esum += iter->second;
    }
  return esum;
}

