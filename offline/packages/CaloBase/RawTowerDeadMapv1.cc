#include "RawTowerDeadMapv1.h"

#include <cstdlib>
#include <iostream>
#include <map>

using namespace std;

const RawTowerDeadMapv1::Map&
RawTowerDeadMapv1::getDeadTowers() const
{
  return m_DeadTowers;
}

RawTowerDeadMapv1::Map&
RawTowerDeadMapv1::getDeadTowers()
{
  return m_DeadTowers;
}

void RawTowerDeadMapv1::addDeadTower(const unsigned int ieta, const int unsigned iphi)
{
  RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(_caloid, ieta, iphi);
  m_DeadTowers.insert(key);
}

void RawTowerDeadMapv1::addDeadTower(RawTowerDefs::keytype key)
{
  if (RawTowerDefs::decode_caloid(key) != _caloid)
  {
    cout << "RawTowerDeadMapv1::addDeadTower - Error - adding tower to wrong container! Container CaloID = "
         << _caloid << ", requested CaloID = " << RawTowerDefs::decode_caloid(key) << " based on key " << key << endl;
    exit(2);
  }
  m_DeadTowers.insert(key);
}

bool RawTowerDeadMapv1::isDeadTower(RawTowerDefs::keytype key)
{
  auto it = m_DeadTowers.find(key);
  if (it != m_DeadTowers.end())
  {
    return true;
  }
  return false;
}

bool RawTowerDeadMapv1::isDeadTower(const unsigned int ieta, const unsigned int iphi)
{
  RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(_caloid, ieta, iphi);
  return isDeadTower(key);
}

int RawTowerDeadMapv1::isValid() const
{
  return size() > 0;
}

void RawTowerDeadMapv1::Reset()
{
  m_DeadTowers.clear();
}

void RawTowerDeadMapv1::identify(std::ostream& os) const
{
  os << "RawTowerDeadMapv1, number of towers: " << size() << std::endl;
}
