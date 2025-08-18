#ifndef CALOBASE_TOWERINFODEFS_H
#define CALOBASE_TOWERINFODEFS_H

#include "RawTowerDefs.h"

namespace TowerInfoDefs
{

  unsigned int encode_emcal(const unsigned int towerIndex);
  unsigned int encode_emcal(const unsigned int etabin, const unsigned int phibin);
  unsigned int decode_emcal(const unsigned int tower_key);

  unsigned int encode_hcal(const unsigned int towerIndex);
  unsigned int encode_hcal(const unsigned int etabin, const unsigned int phibin);
  unsigned int decode_hcal(const unsigned int tower_key);

  unsigned int encode_epd(const unsigned int towerIndex);
  unsigned int encode_epd(const unsigned int arm, const unsigned int rbin, const unsigned int phibin);
  unsigned int decode_epd(const unsigned int tower_key);

  unsigned int get_epd_arm(unsigned int key);
  unsigned int get_epd_sector(unsigned int key);
  unsigned int get_epd_rbin(unsigned int key);
  unsigned int get_epd_phibin(unsigned int key);
  unsigned int getCaloTowerPhiBin(const unsigned int key);
  unsigned int getCaloTowerEtaBin(const unsigned int key);
  std::pair<int, int> getEMCalSectorIB(unsigned int towerIndex);

  unsigned int get_mbd_arm(const unsigned int key);
  unsigned int get_mbd_side(const unsigned int key);  // side is same as arm
  unsigned int get_mbd_type(const unsigned int key);  // 0 for time 1 for charge
  unsigned int get_mbd_channel(const unsigned int key);
  unsigned int encode_mbd(const unsigned int towerIndex);
  unsigned int decode_mbd(const unsigned int key);

  unsigned int encode_zdc(const unsigned int towerIndex);  // has ZDC,SMD,VETO
  unsigned int decode_zdc(const unsigned int key);
  bool isZDC(const unsigned int towerIndex);
  int get_zdc_side(const unsigned int key);
  bool isSMD(const unsigned int towerIndex);
  int get_smd_side(const unsigned int key);
  bool isVeto(const unsigned int towerIndex);
  int get_veto_side(const unsigned int key);

  RawTowerDefs::keytype get_emcal_geokey_at_channel(const unsigned int towerIndex);
  RawTowerDefs::keytype get_hcalin_geokey_at_channel(const unsigned int towerIndex);
  RawTowerDefs::keytype get_hcalout_geokey_at_channel(const unsigned int towerIndex);

}  // namespace TowerInfoDefs
#endif
