#include "RawTowerContainer.h"

#include "RawTower.h"

#include <cstdlib>
#include <iostream>

using namespace std;

void RawTowerContainer::compress(const double emin)
{
  if (emin <= 0)  // no need to loop through the map if we don't apply a cut
  {
    return;
  }
  Iterator itr = _towers.begin();
  Iterator last = _towers.end();
  for (; itr != last;)
  {
    RawTower *tower = (itr->second);
    if (tower->get_energy() < emin)
    {
      delete tower;
      _towers.erase(itr++);
    }
    else
    {
      ++itr;
    }
  }
}

RawTowerContainer::ConstRange
RawTowerContainer::getTowers() const
{
  return make_pair(_towers.begin(), _towers.end());
}

RawTowerContainer::Range
RawTowerContainer::getTowers()
{
  return make_pair(_towers.begin(), _towers.end());
}

RawTowerContainer::ConstIterator
RawTowerContainer::AddTower(const unsigned int ieta, const int unsigned iphi, RawTower *rawtower)
{
  RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(_caloid, ieta, iphi);
  _towers[key] = rawtower;
  rawtower->set_id(key);  // force tower key to be synced to container key

  return _towers.find(key);
}

RawTowerContainer::ConstIterator
RawTowerContainer::AddTower(RawTowerDefs::keytype key, RawTower *twr)
{
  if (RawTowerDefs::decode_caloid(key) != _caloid)
  {
    cout << "RawTowerContainer::AddTower - Error - adding tower to wrong container! Container CaloID = "
         << _caloid << ", requested CaloID = " << RawTowerDefs::decode_caloid(key) << " based on key " << key << endl;
    exit(2);
  }

  _towers[key] = twr;
  twr->set_id(key);  // force tower key to be synced to container key

  return _towers.find(key);
}

RawTower *
RawTowerContainer::getTower(RawTowerDefs::keytype key)
{
  ConstIterator it = _towers.find(key);
  if (it != _towers.end())
  {
    return it->second;
  }
  return nullptr;
}

const RawTower *
RawTowerContainer::getTower(RawTowerDefs::keytype key) const
{
  ConstIterator it = _towers.find(key);
  if (it != _towers.end())
  {
    return it->second;
  }
  return nullptr;
}

RawTower *
RawTowerContainer::getTower(const unsigned int ieta, const unsigned int iphi)
{
  RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(_caloid, ieta, iphi);
  return getTower(key);
}

const RawTower *
RawTowerContainer::getTower(const unsigned int ieta, const unsigned int iphi) const
{
  RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(_caloid, ieta, iphi);
  return getTower(key);
}

RawTower *
RawTowerContainer::getTower(const unsigned int ieta, const unsigned int iphi, const unsigned int il)
{
  RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(_caloid, ieta, iphi, il);
  return getTower(key);
}

const RawTower *
RawTowerContainer::getTower(const unsigned int ieta, const unsigned int iphi, const unsigned int il) const
{
  RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(_caloid, ieta, iphi, il);
  return getTower(key);
}

int RawTowerContainer::isValid() const
{
  return (!_towers.empty());
}

void RawTowerContainer::Reset()
{
  while (_towers.begin() != _towers.end())
  {
    delete _towers.begin()->second;
    _towers.erase(_towers.begin());
  }
}

void RawTowerContainer::identify(std::ostream &os) const
{
  os << "RawTowerContainer, number of towers: " << size() << std::endl;
}

double
RawTowerContainer::getTotalEdep() const
{
  double totalenergy = 0;
  ConstIterator iter;
  for (iter = _towers.begin(); iter != _towers.end(); ++iter)
  {
    totalenergy += iter->second->get_energy();
  }
  return totalenergy;
}
