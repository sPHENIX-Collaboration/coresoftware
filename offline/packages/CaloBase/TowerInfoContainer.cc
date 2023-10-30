#include "TowerInfoDefs.h"
#include "TowerInfoContainer.h"

void TowerInfoContainer::identify(std::ostream& os) const
{
  os << "TowerInfoContainer Base Class " << std::endl;
}


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
