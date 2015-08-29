#include "RawTowerContainer.h"
#include "RawTower.h"

#include <cstdlib>
#include <iostream>

ClassImp(RawTowerContainer)

using namespace std;

unsigned int
RawTowerContainer::genkey(const unsigned int ieta, const unsigned int iphi) const
{
  if (ieta > 0xFFFF || iphi > 0xFFFF)
    {
      cout << "ieta " << ieta << " or iphi " << iphi 
	   << " exceed max length of " << 0xFFFF << endl;
      cout << "reconsider the generation of unique keys" << endl;
      exit(1);
    }
  unsigned int key = 0;
  key |= (ieta << 16);
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
  std::map<unsigned int, RawTower*>::iterator itr = _towers.begin();
  std::map<unsigned int, RawTower*>::iterator last = _towers.end();
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
RawTowerContainer::AddTower(const int ieta, const int iphi, RawTower *rawtower)
{
  unsigned int key = genkey(ieta,iphi);
  _towers[key] = rawtower;
  return _towers.find(key);
}

RawTower *
RawTowerContainer::getTower(const int ieta, const int iphi)
{
  unsigned int key = genkey(ieta,iphi);
  Iterator it = _towers.find(key);
  if (it != _towers.end())
    {
      return it->second;
    }
  return NULL;
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
