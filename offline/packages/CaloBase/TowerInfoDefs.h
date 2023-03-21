#ifndef CALOBASE_TOWERINFODEFS_H
#define CALOBASE_TOWERINFODEFS_H

#include "RawTowerDefs.h"
#include <bitset>
#include <cstdlib>
#include <iostream>
#include <string>
/*! Namespace with functions to encode / decode calo information. 
 */


int emcadc[8][8] = {
    {62, 60, 46, 44, 30, 28, 14, 12},
    {63, 61, 47, 45, 31, 29, 15, 13},
    {58, 56, 42, 40, 26, 24, 10, 8},
    {59, 57, 43, 41, 27, 25, 11, 9},
    {54, 52, 38, 36, 22, 20, 6, 4},
    {55, 53, 39, 37, 23, 21, 7, 5},
    {50, 48, 34, 32, 18, 16, 2, 0},
    {51, 49, 35, 33, 19, 17, 3, 1}};

int hcaladc[8][2] = {
    {0, 1},
    {2, 3},
    {4, 5},
    {6, 7},
    {8, 9},
    {10, 11},
    {12, 13},
    {14, 15}};









namespace TowerInfoDefs
{

  // convert from tower index to key
  inline unsigned int encode_emcal(const unsigned int towerIndex)
  {
    int phimap[64] = {0};
    int etamap[64] = {0};
    for (int j = 0; j < 8; j++)
      {
	for (int k = 0; k < 8; k++)
	  {
	    etamap[emcadc[j][k]] = j;
	    phimap[emcadc[j][k]] = k;
	  }
      }
    int channels_per_sector = 64;
    int supersector = 64 * 12;
    int  nchannelsperpacket = 64 * 3;
    int  maxphibin = 7;
    int  maxetabin = 23;
    int etabinoffset[4] = {24,0,48,72};
    int supersectornumber = towerIndex / supersector;
    int packet = (towerIndex % supersector) / nchannelsperpacket;  // 0 = S small |eta|, 1 == S big |eta|, 2 == N small |eta|, 3 == N big |eta|
    if (packet < 0 || packet > 3 )
      {
	std::cout << "Attempting to access channel with invalid value in EMCal " << packet << std::endl;
	exit(1);
      }
    int interfaceboard = ((towerIndex % supersector) % nchannelsperpacket) / channels_per_sector;
    int interfaceboard_channel = ((towerIndex % supersector) % nchannelsperpacket) % channels_per_sector; 
    int localphibin = phimap[interfaceboard_channel];
    if (packet == 0 || packet == 1)
      {
	localphibin = maxphibin - localphibin;
      } 
    int localetabin = etamap[interfaceboard_channel];
    int packet_etabin = localetabin + 8 * interfaceboard;
    if (packet == 0 || packet == 1)
      {
	packet_etabin = maxetabin - packet_etabin;
      }
    unsigned int globaletabin = packet_etabin + etabinoffset[packet];
    unsigned int globalphibin = localphibin + supersectornumber * 8;
    unsigned int key = globalphibin + (globaletabin << 16U);
    return key;
  } 




  inline unsigned int encode_emcal (const unsigned int etabin, const unsigned int phibin)
  {
    unsigned int key = phibin + (etabin << 16U);
    return key;
  }




  // convert from tower index to key
  inline unsigned int encode_hcal(const unsigned int towerIndex)
  {
 int phimap[64] = {0};
  int etamap[64] = {0};
  for (int j = 0; j < 8; j++)
    {
      for (int k = 0; k < 2; k++)
      {
        etamap[hcaladc[j][k]] = j;
        phimap[hcaladc[j][k]] = k;
      }
    }
  int channels_per_sector = 16;
  int supersector = 16 * 4 * 3;
  int nchannelsperpacket = channels_per_sector * 4;
  int etabinoffset[4] = {0,8,16,0};
  int phibinoffset[4] = {0,2,4,6};
  int supersectornumber = towerIndex / supersector;
  int packet = (towerIndex % supersector) / nchannelsperpacket;  // 0 = S small |eta|, 1 == S big |eta|, 2 == N small |eta|, 3 == N big |eta|
  if (packet < 0 || packet > 3 )
    {
      std::cout << "Attempting to access channel with invalid value ih HCAL " << packet << std::endl;
      exit(1);
    }

  int interfaceboard = ((towerIndex % supersector) % nchannelsperpacket) / channels_per_sector;
  int interfaceboard_channel = ((towerIndex % supersector) % nchannelsperpacket) % channels_per_sector;
  int localphibin = phimap[interfaceboard_channel] + phibinoffset[interfaceboard];
  int localetabin = etamap[interfaceboard_channel];
  int packet_etabin = localetabin;
  unsigned int globaletabin = packet_etabin + etabinoffset[packet];
  unsigned int globalphibin = localphibin + supersectornumber * 8;
  unsigned int key = globalphibin + (globaletabin << 16U);
  return key;

  } 
  inline unsigned int encode_hcal (const unsigned int etabin, const unsigned int phibin)
  {
    unsigned int key = phibin + (etabin << 16U);
    return key;
  }


  inline unsigned int encode_epd(const unsigned int towerIndex)  // convert from tower index to key
  {
    int channels_per_sector = 31;
    int supersector = channels_per_sector * 12;
    unsigned int supersectornumber = towerIndex / supersector;
    int sector = ((towerIndex % supersector)) / channels_per_sector;
    int channel = ((towerIndex % supersector)) % channels_per_sector;
    int rmap[31] = {0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14, 14, 15, 15};
    int phimap_sepd[31] = {0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1};
    unsigned int globalphi = phimap_sepd[channel] + 2 * sector;
    unsigned int r = rmap[channel];
    unsigned int key = globalphi + (r << 10U) + (supersectornumber << 20U);
    return key;
    
  } 

  inline unsigned int encode_epd (const unsigned int arm, const unsigned int rbin, const unsigned int phibin) 
  {
    unsigned int key = phibin + (rbin << 10U) + (arm << 20U);
    return key;
  }


  inline unsigned int decode_epd(const unsigned int tower_key) // decode the tower_key to corresponding index for the emcal
  {
    int channels_per_sector = 31;
    int supersector = channels_per_sector * 12;
    unsigned int ns_sector = tower_key >> 20U;
    unsigned int rbin = (tower_key - (ns_sector << 20U)) >> 10U;
    unsigned int phibin = tower_key - (ns_sector << 20U) - (rbin << 10U);
    int epdchnlmap[16][2] = {{0, 0}, {1, 2}, {3, 4}, {5, 6}, {7, 8}, {9, 10}, {11, 12}, {13, 14}, {15, 16}, {17, 18}, {19, 20}, {21, 22}, {23, 24}, {25, 26}, {27, 28}, {29, 30}};
    int sector = phibin / 2;
    int channel = 0;
    if (rbin > 0)
      {
	channel = epdchnlmap[rbin][phibin - 2 * sector];
      }
    else
      {
	channel = 0;
      }
    unsigned int index = 0;
    index = ns_sector * supersector + sector * channels_per_sector + channel;
    return index;
  }


  inline unsigned int decode_emcal(const unsigned int tower_key)  //convert from calorimeter key to channel index
  {
    int etabinoffset[4] = {0};
    int etabinmap[4] = {0};
    int channels_per_sector = 64;
    int supersector = 64 * 12;
    int nchannelsperpacket = 64 * 3;
    int maxphibin = 7;
    int maxetabin = 23;
    etabinoffset[0] = 24;
    etabinoffset[1] = 0;
    etabinoffset[2] = 48;
    etabinoffset[3] = 72;
    
    etabinmap[0] = 1;
    etabinmap[1] = 0;
    etabinmap[2] = 2;
    etabinmap[3] = 3;
    unsigned int etabin = tower_key >> 16U;
    unsigned int phibin = tower_key - (etabin << 16U);
    int packet = etabinmap[(int) etabin / 24];
    int localetabin = etabin - etabinoffset[packet];
    int localphibin = phibin % 8;
    int supersectornumber = phibin / 8;
    int ib = 0;
    if (packet == 0 || packet == 1)
      {
	localetabin = maxetabin - localetabin;
      }
    ib = localetabin / 8;
    unsigned int index = 0;
    if (packet == 0 || packet == 1)
      {
	localphibin = maxphibin - localphibin;
      }
    localetabin = localetabin % 8;
    unsigned int localindex = emcadc[localetabin][localphibin];
    index = localindex + channels_per_sector * ib + packet * nchannelsperpacket + supersector * supersectornumber;
    return index;
  }
  
  inline unsigned int decode_hcal(const unsigned int tower_key) //convert from calorimeter key to channel index
  {
    int channels_per_sector = 16;
    int supersector = 16 * 4 * 3;
    int nchannelsperpacket = channels_per_sector * 4;
    int etabinoffset[3] = {0,8,16};
    int phibinoffset[4] = {0,2,4,6};
    unsigned int etabin = tower_key >> 16U;
    unsigned int phibin = tower_key - (etabin << 16U);
    int packet = etabin / 8;
    int localetabin = etabin - etabinoffset[packet];
    int localphibin = phibin % 8;
    int supersectornumber = phibin / 8;
    int ib = localphibin / 2;
    unsigned int index = 0;
    localphibin = localphibin - phibinoffset[ib];
    unsigned int localindex = hcaladc[localetabin][localphibin];
    index = localindex + channels_per_sector * ib + packet * nchannelsperpacket + supersector * supersectornumber;
    return index;
  }


  inline unsigned int getCaloTowerPhiBin(const unsigned int key) // convert from calorimeter key to phi bin coordinates
  {
    unsigned int etabin = key >> 16U;
    unsigned int phibin = key - (etabin << 16U);
    return phibin;
  }
  
  inline unsigned int getCaloTowerEtaBin(const unsigned int key) // convert from calorimeter key to eta bin coordinates
  {
    unsigned int etabin = key >> 16U;
    return etabin;
  }
  

  inline RawTowerDefs::keytype get_emcal_geokey_at_channel(const unsigned int towerIndex) // convienent for interface to geometry class
  {
    unsigned int towerkey = encode_emcal(towerIndex);
    unsigned int etabin = getCaloTowerEtaBin(towerkey);
    unsigned int phibin = getCaloTowerPhiBin(towerkey);
    const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::CEMC, etabin, phibin);
    return key;
  }


  inline RawTowerDefs::keytype get_hcalin_geokey_at_channel(const unsigned int towerIndex) // convienent for interface to geometry class
  {
    unsigned int towerkey = encode_hcal(towerIndex);
    unsigned int etabin = getCaloTowerEtaBin(towerkey);
    unsigned int phibin = getCaloTowerPhiBin(towerkey);
    const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALIN, etabin, phibin);
    return key;
  }


  inline RawTowerDefs::keytype get_hcalout_geokey_at_channel(const unsigned int towerIndex) // convienent for interface to geometry class
  {
    unsigned int towerkey = encode_hcal(towerIndex);
    unsigned int etabin = getCaloTowerEtaBin(towerkey);
    unsigned int phibin = getCaloTowerPhiBin(towerkey);
    const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALOUT, etabin, phibin);
    return key;
  }









}
#endif
