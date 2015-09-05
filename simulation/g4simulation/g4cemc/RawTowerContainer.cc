#include "RawTowerContainer.h"
#include "RawTower.h"

#include <cstdlib>
#include <iostream>

ClassImp(RawTowerContainer)

using namespace std;

RawTowerDefs::keytype
RawTowerContainer::genkey(const unsigned int ieta, const unsigned int iphi) const
{
  if (ieta > 0xFFF || iphi > 0xFFF)
    {
      cout << "ieta " << ieta << " or iphi " << iphi 
	   << " exceed max length of " << 0xFFF << endl;
      cout << "reconsider the generation of unique keys" << endl;
      exit(1);
    }
  RawTowerDefs::keytype key = 0;
  key |= (ieta << RawTowerDefs::eta_idbits);
  key |= iphi;
  return key;
}

void 
RawTowerContainer::compress(const double emin)
{
  if (emin <= 0) // no need to loop through the map if we don't apply a cut
    {
      return;
    }
  Iterator itr = _towers.begin();
  Iterator last = _towers.end();
  for (; itr != last; )
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
RawTowerContainer::getTowers( void ) const
{
  return make_pair(_towers.begin(), _towers.end());
}


RawTowerContainer::Range
RawTowerContainer::getTowers( void )
{
  return make_pair(_towers.begin(), _towers.end());
}


RawTowerContainer::ConstIterator
RawTowerContainer::AddTower(const unsigned int ieta, const int unsigned iphi, RawTower *rawtower)
{
  RawTowerDefs::keytype key = genkey(ieta,iphi);
  _towers[key] = rawtower;
  return _towers.find(key);
}

RawTower *
RawTowerContainer::getTower(RawTowerDefs::keytype key)
{
  Iterator it = _towers.find(key);
  if (it != _towers.end())
    {
      return it->second;
    }
  return NULL;
}

RawTower *
RawTowerContainer::getTower(const unsigned int ieta, const unsigned int iphi)
{
  RawTowerDefs::keytype key = genkey(ieta,iphi);
  return getTower(key);
}

int 
RawTowerContainer::isValid() const
{
  return (!_towers.empty());
}

void
RawTowerContainer::Reset()
{
  while (_towers.begin() != _towers.end())
    {
      delete _towers.begin()->second;
      _towers.erase(_towers.begin());
    }
}

void 
RawTowerContainer::identify(std::ostream& os) const
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
