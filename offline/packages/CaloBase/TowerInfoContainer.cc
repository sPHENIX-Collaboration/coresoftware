#include "TowerInfoDefs.h"
#include "TowerInfoContainer.h"



unsigned int TowerInfoContainer::encode_epd(unsigned int towerIndex)
{
  unsigned int key = TowerInfoDefs::encode_epd(towerIndex);
  return key;
}

unsigned int TowerInfoContainer::encode_emcal(unsigned int towerIndex)
{
  unsigned int key = TowerInfoDefs::encode_emcal(towerIndex);
  return key;
}

unsigned int TowerInfoContainer::encode_hcal(unsigned int towerIndex)
{
  unsigned int key = TowerInfoDefs::encode_hcal(towerIndex);
  return key;
}

unsigned int TowerInfoContainer::encode_mbd(unsigned int towerIndex)
{
  unsigned int key = TowerInfoDefs::encode_mbd(towerIndex);
  return key;
}
unsigned int TowerInfoContainer::encode_zdc(unsigned int towerIndex)
{
  unsigned int key = TowerInfoDefs::encode_zdc(towerIndex);
  return key;
}

unsigned int TowerInfoContainer::encode_key(unsigned int towerIndex)
{
  int key = 0;
  if (_detector == DETECTOR::EMCAL)
    {
      key = TowerInfoContainer::encode_emcal(towerIndex);
    }
  else if (_detector == DETECTOR::HCAL)
    {
      key = TowerInfoContainer::encode_hcal(towerIndex);
    }
  else if (_detector == DETECTOR::SEPD)
    {
    key = TowerInfoContainer::encode_epd(towerIndex);
    }
  else if (_detector == DETECTOR::MBD)
    {
    key = TowerInfoContainer::encode_mbd(towerIndex);
    }
  else if (_detector == DETECTOR::ZDC)
    {
    key = TowerInfoContainer::encode_zdc(towerIndex);
    }
  return key;
}

unsigned int TowerInfoContainer::decode_epd(unsigned int tower_key)
{
  unsigned int index = TowerInfoDefs::decode_epd(tower_key);  
  return index;
}

unsigned int TowerInfoContainer::decode_emcal(unsigned int tower_key)
{

  unsigned int index = TowerInfoDefs::decode_emcal(tower_key);  
  return index;
}

unsigned int TowerInfoContainer::decode_hcal(unsigned int tower_key)
{
  unsigned int index = TowerInfoDefs::decode_hcal(tower_key);  
  return index;
}

unsigned int TowerInfoContainer::decode_mbd(unsigned int tower_key)
{
  unsigned int index = TowerInfoDefs::decode_mbd(tower_key);  
  return index;
}
unsigned int TowerInfoContainer::decode_zdc(unsigned int tower_key)
{
  unsigned int index = TowerInfoDefs::decode_zdc(tower_key);  
  return index;
}

unsigned int TowerInfoContainer::decode_key(unsigned int tower_key)
{
  int index = 0;

  if (_detector == DETECTOR::EMCAL)
  {
    index = TowerInfoContainer::decode_emcal(tower_key);
  }
  else if (_detector == DETECTOR::HCAL)
  {
    index = TowerInfoContainer::decode_hcal(tower_key);
  }
  else if (_detector == DETECTOR::SEPD)
  {
    index = TowerInfoContainer::decode_epd(tower_key);
  }
  else if (_detector == DETECTOR::MBD)
  {
    index = TowerInfoContainer::decode_mbd(tower_key);
  }
  else if (_detector == DETECTOR::ZDC)
  {
    index = TowerInfoContainer::decode_zdc(tower_key);
  }
  return index;
}

unsigned int TowerInfoContainer::getTowerPhiBin(unsigned int key)
{
  unsigned int phibin = TowerInfoDefs::getCaloTowerPhiBin(key);
  return phibin;
}

unsigned int TowerInfoContainer::getTowerEtaBin(unsigned int key)
{
  unsigned int etabin = TowerInfoDefs::getCaloTowerEtaBin(key);
  return etabin;
}