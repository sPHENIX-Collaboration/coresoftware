#include "TowerInfoContainer.h"
#include "TowerInfoDefs.h"

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

unsigned int TowerInfoContainer::encode_key(unsigned int towerIndex)
{
  unsigned int key = 0;
  if (get_detectorid() == DETECTOR::EMCAL)
  {
    key = TowerInfoContainer::encode_emcal(towerIndex);
  }
  else if (get_detectorid() == DETECTOR::HCAL)
  {
    key = TowerInfoContainer::encode_hcal(towerIndex);
  }
  else if (get_detectorid() == DETECTOR::SEPD)
  {
    key = TowerInfoContainer::encode_epd(towerIndex);
  }
  else if (get_detectorid() == DETECTOR::MBD)
  {
    key = TowerInfoContainer::encode_mbd(towerIndex);
  }
  else if (get_detectorid() == DETECTOR::ZDC)
  {
    key = TowerInfoContainer::encode_zdc(towerIndex);
  }
  return key;
}

unsigned int TowerInfoContainer::decode_key(unsigned int tower_key)
{
  unsigned int index = 0;

  if (get_detectorid() == DETECTOR::EMCAL)
  {
    index = TowerInfoContainer::decode_emcal(tower_key);
  }
  else if (get_detectorid() == DETECTOR::HCAL)
  {
    index = TowerInfoContainer::decode_hcal(tower_key);
  }
  else if (get_detectorid() == DETECTOR::SEPD)
  {
    index = TowerInfoContainer::decode_epd(tower_key);
  }
  else if (get_detectorid() == DETECTOR::MBD)
  {
    index = TowerInfoContainer::decode_mbd(tower_key);
  }
  else if (get_detectorid() == DETECTOR::ZDC)
  {
    index = TowerInfoContainer::decode_zdc(tower_key);
  }
  return index;
}

int TowerInfoContainer::get_channels(DETECTOR detec)
{
  int nchannels = 744;
  if (detec == DETECTOR::SEPD)
  {
    nchannels = 744;
  }
  else if (detec == DETECTOR::EMCAL)
  {
    nchannels = 24576;
  }
  else if (detec == DETECTOR::HCAL)
  {
    nchannels = 1536;
  }
  else if (detec == DETECTOR::MBD)
  {
    nchannels = 256;
  }
  else if (detec == DETECTOR::ZDC)
  {
    nchannels = 52;
  }
  return nchannels;
}
