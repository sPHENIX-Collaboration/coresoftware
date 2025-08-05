#include "TowerInfoDefs.h"
#include "RawTowerDefs.h"

#include <cstdlib>
#include <iostream>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
static const int emcadc[8][8] = {
    {62, 60, 46, 44, 30, 28, 14, 12},
    {63, 61, 47, 45, 31, 29, 15, 13},
    {58, 56, 42, 40, 26, 24, 10, 8},
    {59, 57, 43, 41, 27, 25, 11, 9},
    {54, 52, 38, 36, 22, 20, 6, 4},
    {55, 53, 39, 37, 23, 21, 7, 5},
    {50, 48, 34, 32, 18, 16, 2, 0},
    {51, 49, 35, 33, 19, 17, 3, 1}};

static const int hcaladc[8][2] = {
    {0, 1},
    {2, 3},
    {4, 5},
    {6, 7},
    {8, 9},
    {10, 11},
    {12, 13},
    {14, 15}};
static const int epdchnlmap[16][2] = {
    {0, 0},
    {1, 2},
    {3, 4},
    {5, 6},
    {7, 8},
    {9, 10},
    {11, 12},
    {13, 14},
    {15, 16},
    {17, 18},
    {19, 20},
    {21, 22},
    {23, 24},
    {25, 26},
    {27, 28},
    {29, 30}};

static const int epd_phimap[31] = {0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1};
static const int epd_rmap[31] = {0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14, 14, 15, 15};

unsigned int TowerInfoDefs::encode_emcal(const unsigned int towerIndex)
{
  static int phimap[64];
  static int etamap[64];
  static int etabinoffset[4];
  static int ifirst = 1;
  if (ifirst == 1)
  {
    for (int j = 0; j < 8; j++)
    {
      for (int k = 0; k < 8; k++)
      {
        etamap[emcadc[j][k]] = j;
        phimap[emcadc[j][k]] = k;
      }
    }
    etabinoffset[0] = 24;
    etabinoffset[1] = 0;
    etabinoffset[2] = 48;
    etabinoffset[3] = 72;
    ifirst = 0;
  }
  const int channels_per_sector = 64;
  const int supersector = 64 * 12;
  const int nchannelsperpacket = 64 * 3;
  const int maxphibin = 7;
  const int maxetabin = 23;
  int supersectornumber = towerIndex / supersector;
  int packet = (towerIndex % supersector) / nchannelsperpacket;  // 0 = S small |eta|, 1 == S big |eta|, 2 == N small |eta|, 3 == N big |eta|
  if (packet < 0 || packet > 3)
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

unsigned int TowerInfoDefs::encode_emcal(const unsigned int etabin, const unsigned int phibin)
{
  unsigned int key = phibin + (etabin << 16U);
  return key;
}

unsigned int TowerInfoDefs::decode_emcal(const unsigned int tower_key)
{
  const int etabinoffset[4] = {24, 0, 48, 72};
  const int etabinmap[4] = {1, 0, 2, 3};
  const int channels_per_sector = 64;
  const int supersector = 64 * 12;
  const int nchannelsperpacket = 64 * 3;
  const int maxphibin = 7;
  const int maxetabin = 23;

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

unsigned int TowerInfoDefs::encode_hcal(const unsigned int towerIndex)
{
  static int phimap[64];
  static int etamap[64];
  static int etabinoffset[4];
  static int phibinoffset[4];
  static int ifirst = 1;
  if (ifirst == 1)
  {
    for (int j = 0; j < 8; j++)
    {
      for (int k = 0; k < 2; k++)
      {
        etamap[hcaladc[j][k]] = j;
        phimap[hcaladc[j][k]] = k;
      }
    }
    etabinoffset[0] = 0;
    etabinoffset[1] = 8;
    etabinoffset[2] = 16;
    etabinoffset[3] = 0;

    phibinoffset[0] = 0;
    phibinoffset[1] = 2;
    phibinoffset[2] = 4;
    phibinoffset[3] = 6;
    ifirst = 0;
  }

  const int channels_per_sector = 16;
  const int supersector = 16 * 4 * 3;
  const int nchannelsperpacket = channels_per_sector * 4;
  // const int etabinoffset[4] = {0,8,16,0};
  // const int phibinoffset[4] = {0,2,4,6};
  int supersectornumber = towerIndex / supersector;
  int packet = (towerIndex % supersector) / nchannelsperpacket;  // 0 = S small |eta|, 1 == S big |eta|, 2 == N small |eta|, 3 == N big |eta|
  if (packet < 0 || packet > 3)
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

// convert from etabin-phibin to key
unsigned int TowerInfoDefs::encode_hcal(const unsigned int etabin, const unsigned int phibin)
{
  unsigned int key = phibin + (etabin << 16U);
  return key;
}

unsigned int TowerInfoDefs::decode_hcal(const unsigned int tower_key)
{
  int channels_per_sector = 16;
  int supersector = 16 * 4 * 3;
  int nchannelsperpacket = channels_per_sector * 4;
  int etabinoffset[3] = {0, 8, 16};
  int phibinoffset[4] = {0, 2, 4, 6};
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

// convert from calorimeter key to phi bin
unsigned int TowerInfoDefs::getCaloTowerPhiBin(const unsigned int key)
{
  unsigned int etabin = key >> 16U;
  unsigned int phibin = key - (etabin << 16U);
  return phibin;
}

// convert from calorimeter key to eta bin
unsigned int TowerInfoDefs::getCaloTowerEtaBin(const unsigned int key)
{
  unsigned int etabin = key >> 16U;
  return etabin;
}

// convert from emcal tower index to (sector, interface board)
std::pair<int, int> TowerInfoDefs::getEMCalSectorIB(const unsigned int towerIndex)
{
  constexpr int channels_per_sector = 384;
  constexpr int channels_per_ib = 64;

  int k = towerIndex / channels_per_sector;

  int sector = (k % 2) ? (k - 1) / 2 : (k / 2) + 32;
  int ib = (towerIndex % channels_per_sector) / channels_per_ib;

  return std::make_pair(sector, ib);
}

unsigned int TowerInfoDefs::encode_epd(const unsigned int towerIndex)  // convert from tower index to key
{
  constexpr unsigned int channels_per_sector = 31;
  constexpr unsigned int supersector = channels_per_sector * 12;
  unsigned int supersectornumber = towerIndex / supersector;
  unsigned int sector = ((towerIndex % supersector)) / channels_per_sector;
  unsigned int channel = ((towerIndex % supersector)) % channels_per_sector;
  unsigned int key = channel + (sector << 5U) + (supersectornumber << 9U);
  return key;
}

// convert from arm-rbin-phibin to key
unsigned int TowerInfoDefs::encode_epd(const unsigned int arm, const unsigned int rbin, const unsigned int phibin)
{
  if (rbin == 0 && phibin > 11)
  {
    std::cout << __PRETTY_FUNCTION__ << " encode_epd invalid phibin value: " << phibin << " where max valid phibin is 11" << std::endl;
    exit(1);
  }

  unsigned int sector = phibin / 2;
  if (rbin == 0)
  {
    sector = phibin;
  }

  int channel = 0;
  if (rbin != 0)
  {
    channel = epdchnlmap[rbin][phibin - 2 * sector];
  }

  unsigned int key = channel + (sector << 5U) + (arm << 9U);
  return key;
}

unsigned int TowerInfoDefs::decode_epd(const unsigned int tower_key)
{
  int channels_per_sector = 31;
  int supersector = channels_per_sector * 12;
  unsigned int ns_sector = tower_key >> 9U;
  unsigned int sector = (tower_key - (ns_sector << 9U)) >> 5U;
  unsigned int channel = tower_key - (ns_sector << 9U) - (sector << 5U);
  unsigned int index = ns_sector * supersector + sector * channels_per_sector + channel;
  return index;
}

// convert from epd key to arm bin
unsigned int TowerInfoDefs::get_epd_arm(unsigned int key)
{
  unsigned int arm = key >> 9U;
  return arm;
}
// convert from epd key to sector number
unsigned int TowerInfoDefs::get_epd_sector(unsigned int key)
{
  unsigned int arm = get_epd_arm(key);
  unsigned int sector = (key - (arm << 9U)) >> 5U;
  return sector;
}
// convert from epd key to r bin
unsigned int TowerInfoDefs::get_epd_rbin(unsigned int key)
{
  unsigned int arm = get_epd_arm(key);
  unsigned int sector = get_epd_sector(key);
  unsigned int channel = key - (sector << 5U) - (arm << 9U);
  unsigned int rbin = epd_rmap[channel];
  return rbin;
}
// convert from epd key to phi bin
unsigned int TowerInfoDefs::get_epd_phibin(unsigned int key)
{
  int flip[24] = {1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1};
  unsigned int arm = get_epd_arm(key);
  unsigned int rbin = get_epd_rbin(key);
  unsigned int sector = get_epd_sector(key);
  unsigned int channel = key - (sector << 5U) - (arm << 9U);
  unsigned int phibin = epd_phimap[channel] + 2 * sector;
    
  if (arm == 1)
  {
    phibin = phibin + flip[phibin];
  }
    
  if (rbin == 0)
  {
    phibin = sector;
  }

  return phibin;
}

unsigned int TowerInfoDefs::encode_zdc(const unsigned int towerIndex)
{
  if (towerIndex > 51)
  {
    std::cout << "Attempting to access zdc channel with invalid number " << towerIndex << std::endl;
    exit(1);
  }
  unsigned int key = towerIndex;
  return key;
}

unsigned int TowerInfoDefs::decode_zdc(const unsigned int key)
{
  unsigned int index = 999;
  for (unsigned int i = 0; i < 52; i++)
  {
    if (encode_zdc(i) == key)
    {
      index = i;
      break;
    }
  }
  return index;
}

bool TowerInfoDefs::isZDC(const unsigned int towerIndex)
{
  bool is_zdc = false;

  if (towerIndex < 16)
  {
    is_zdc = true;
  }
  return is_zdc;
}

// get zdc side, 0 = south, 1 = north
int TowerInfoDefs::get_zdc_side(const unsigned int key)
{
  if (key & 8U)
  {
    return 1;
  }
  return 0;
}

bool TowerInfoDefs::isSMD(const unsigned int towerIndex)
{
  bool is_smd = false;

  if ((towerIndex > 17 && towerIndex < 34) || (towerIndex > 35 && towerIndex < 52))
  {
    is_smd = true;
  }
  return is_smd;
}

// get smd side, 0 = south, 1 = north
int TowerInfoDefs::get_smd_side(const unsigned int key)
{
  if (key < 34)
  {
    return 1;
  }
  return 0;
}

bool TowerInfoDefs::isVeto(const unsigned int towerIndex)
{
  bool is_veto = false;

  if ((towerIndex > 15 && towerIndex < 18) || (towerIndex > 33 && towerIndex < 36))
  {
    is_veto = true;
  }
  return is_veto;
}
// get veto side, 0 = south, 1 = north
int TowerInfoDefs::get_veto_side(const unsigned int key)
{
  if (key & 2U)
  {
    return 0;
  }
  return 1;
}

// 128 channels per side, goes 8 times and 8 charges and so on
unsigned int TowerInfoDefs::encode_mbd(const unsigned int pmtIndex)
{
  unsigned int arm = pmtIndex / 128;
  unsigned int type = (pmtIndex % 16) / 8;
  unsigned int channel = (pmtIndex % 8) + ((pmtIndex / 16) * 8);
  if (channel > 63)
  {
    channel -= 64;
  }

  unsigned int key = (arm << 7U) | (type << 6U) | channel;

  return key;
}

unsigned int TowerInfoDefs::decode_mbd(const unsigned int key)
{
  unsigned int arm = (key >> 7U) & 0x1U;
  unsigned int type = (key >> 6U) & 0x1U;
  unsigned int channel = key & 0x3fU;

  unsigned int index = (arm * 128) + (type * 8) + (channel % 8) + (channel / 8) * 16;

  return index;
}

unsigned int TowerInfoDefs::get_mbd_arm(const unsigned int key)
{
  return (key >> 7U) & 0x1U;
}

unsigned int TowerInfoDefs::get_mbd_side(const unsigned int key)
{
  return get_mbd_arm(key);
}

unsigned int TowerInfoDefs::get_mbd_type(const unsigned int key)
{
  return (key >> 6U) & 0x1U;
}

unsigned int TowerInfoDefs::get_mbd_channel(const unsigned int key)
{
  return key & 0x3fU;
}

// convienent for interface to geometry class
RawTowerDefs::keytype TowerInfoDefs::get_emcal_geokey_at_channel(const unsigned int towerIndex)
{
  unsigned int towerkey = encode_emcal(towerIndex);
  unsigned int etabin = getCaloTowerEtaBin(towerkey);
  unsigned int phibin = getCaloTowerPhiBin(towerkey);
  const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::CEMC, etabin, phibin);
  return key;
}

// convienent for interface to geometry class
RawTowerDefs::keytype TowerInfoDefs::get_hcalin_geokey_at_channel(const unsigned int towerIndex)
{
  unsigned int towerkey = encode_hcal(towerIndex);
  unsigned int etabin = getCaloTowerEtaBin(towerkey);
  unsigned int phibin = getCaloTowerPhiBin(towerkey);
  const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALIN, etabin, phibin);
  return key;
}

// convienent for interface to geometry class
RawTowerDefs::keytype TowerInfoDefs::get_hcalout_geokey_at_channel(const unsigned int towerIndex)
{
  unsigned int towerkey = encode_hcal(towerIndex);
  unsigned int etabin = getCaloTowerEtaBin(towerkey);
  unsigned int phibin = getCaloTowerPhiBin(towerkey);
  const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALOUT, etabin, phibin);
  return key;
}

#pragma GCC diagnostic pop
