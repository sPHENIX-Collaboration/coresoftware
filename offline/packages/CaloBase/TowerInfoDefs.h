#ifndef CALOBASE_TOWERINFODEFS_H
#define CALOBASE_TOWERINFODEFS_H

#include "RawTowerDefs.h"


namespace TowerInfoDefs 
{

   unsigned int encode_emcal(const unsigned int towerIndex);
   unsigned int encode_emcal (const unsigned int etabin, const unsigned int phibin);
   unsigned int decode_emcal(const unsigned int tower_key);

   unsigned int encode_hcal(const unsigned int towerIndex);
   unsigned int encode_hcal (const unsigned int etabin, const unsigned int phibin);
   unsigned int decode_hcal(const unsigned int tower_key);

   unsigned int encode_epd(const unsigned int towerIndex);
   unsigned int encode_epd (const unsigned int arm, const unsigned int rbin, const unsigned int phibin);
   unsigned int decode_epd(const unsigned int tower_key);




   unsigned int get_epd_arm(unsigned int key);
   unsigned int get_epd_sector(unsigned int key);
   unsigned int get_epd_rbin(unsigned int key);
   unsigned int get_epd_phibin(unsigned int key);
   unsigned int getCaloTowerPhiBin(const unsigned int key);
   unsigned int getCaloTowerEtaBin(const unsigned int key) ;

   int get_zdc_side(const unsigned int key) ;
   unsigned int get_zdc_module_index(const unsigned int key);
   int get_smd_side(const unsigned int key);
   int get_smd_xy(const unsigned int key);
   int get_smd_finger_index(const unsigned int key);



   unsigned int encode_zdc(const unsigned int towerIndex);
   unsigned int encode_smd(const unsigned int towerIndex);
   unsigned int decode_smd(const unsigned int key);
   unsigned int decode_zdc(const unsigned int key);


   RawTowerDefs::keytype get_emcal_geokey_at_channel(const unsigned int towerIndex);
   RawTowerDefs::keytype get_hcalin_geokey_at_channel(const unsigned int towerIndex) ;
   RawTowerDefs::keytype get_hcalout_geokey_at_channel(const unsigned int towerIndex) ;



}
#endif
