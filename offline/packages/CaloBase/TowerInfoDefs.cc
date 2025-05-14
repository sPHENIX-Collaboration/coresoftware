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

  static constexpr int epdchnlmap_data[2][16][24] = {
    { // arm 0
      {  618,  432,  401,  556,  525,  711,  463,  494,  649,  680,  742,  587,  999,  999,  999,  999,  999,  999,  999,  999,  999,  999,  999,  999 },
      {  603,  619,  417,  433,  386,  402,  541,  557,  510,  526,  696,  712,  448,  464,  479,  495,  634,  650,  665,  681,  727,  743,  572,  588 },
      {  601,  616,  415,  430,  384,  399,  539,  554,  508,  523,  694,  709,  446,  461,  477,  492,  632,  647,  663,  678,  725,  740,  570,  585 },
      {  602,  617,  416,  431,  385,  400,  540,  555,  509,  524,  695,  710,  447,  462,  478,  493,  633,  648,  664,  679,  726,  741,  571,  586 },
      {  599,  614,  413,  428,  382,  397,  537,  552,  506,  521,  692,  707,  444,  459,  475,  490,  630,  645,  661,  676,  723,  738,  568,  583 },
      {  600,  615,  414,  429,  383,  398,  538,  553,  507,  522,  693,  708,  445,  460,  476,  491,  631,  646,  662,  677,  724,  739,  569,  584 },
      {  597,  612,  411,  426,  380,  395,  535,  550,  504,  519,  690,  705,  442,  457,  473,  488,  628,  643,  659,  674,  721,  736,  566,  581 },
      {  598,  613,  412,  427,  381,  396,  536,  551,  505,  520,  691,  706,  443,  458,  474,  489,  629,  644,  660,  675,  722,  737,  567,  582 },
      {  595,  610,  409,  424,  378,  393,  533,  548,  502,  517,  688,  703,  440,  455,  471,  486,  626,  641,  657,  672,  719,  734,  564,  579 },
      {  596,  611,  410,  425,  379,  394,  534,  549,  503,  518,  689,  704,  441,  456,  472,  487,  627,  642,  658,  673,  720,  735,  565,  580 },
      {  593,  608,  407,  422,  376,  391,  531,  546,  500,  515,  686,  701,  438,  453,  469,  484,  624,  639,  655,  670,  717,  732,  562,  577 },
      {  594,  609,  408,  423,  377,  392,  532,  547,  501,  516,  687,  702,  439,  454,  470,  485,  625,  640,  656,  671,  718,  733,  563,  578 },
      {  591,  606,  405,  420,  374,  389,  529,  544,  498,  513,  684,  699,  436,  451,  467,  482,  622,  637,  653,  668,  715,  730,  560,  575 },
      {  592,  607,  406,  421,  375,  390,  530,  545,  499,  514,  685,  700,  437,  452,  468,  483,  623,  638,  654,  669,  716,  731,  561,  576 },
      {  589,  604,  403,  418,  372,  387,  527,  542,  496,  511,  682,  697,  434,  449,  465,  480,  620,  635,  651,  666,  713,  728,  558,  573 },
      {  590,  605,  404,  419,  373,  388,  528,  543,  497,  512,  683,  698,  435,  450,  466,  481,  621,  636,  652,  667,  714,  729,  559,  574 }
    },
    { // arm 1
      {  277,  370,  339,  246,  215,  122,   60,   29,   91,  153,  184,  308,  999,  999,  999,  999,  999,  999,  999,  999,  999,  999,  999,  999 },
      {  278,  262,  371,  355,  340,  324,  247,  231,  216,  200,  123,  107,   61,   45,   30,    0,   92,   76,  154,  138,  185,  169,  309,  293 },
      {  275,  260,  368,  353,  337,  322,  244,  229,  213,  198,  120,  105,   58,   43,   27,   13,   89,   74,  151,  136,  182,  167,  306,  291 },
      {  276,  261,  369,  354,  338,  323,  245,  230,  214,  199,  121,  106,   59,   44,   28,   14,   90,   75,  152,  137,  183,  168,  307,  292 },
      {  273,  258,  366,  351,  335,  320,  242,  227,  211,  196,  118,  103,   56,   41,   25,   11,   87,   72,  149,  134,  180,  165,  304,  289 },
      {  274,  259,  367,  352,  336,  321,  243,  228,  212,  197,  119,  104,   57,   42,   26,   12,   88,   73,  150,  135,  181,  166,  305,  290 },
      {  271,  256,  364,  349,  333,  318,  240,  225,  209,  194,  116,  101,   54,   39,   23,    9,   85,   70,  147,  132,  178,  163,  302,  287 },
      {  272,  257,  365,  350,  334,  319,  241,  226,  210,  195,  117,  102,   55,   40,   24,   10,   86,   71,  148,  133,  179,  164,  303,  288 },
      {  269,  254,  362,  347,  331,  316,  238,  223,  207,  192,  114,   99,   52,   37,   21,    7,   83,   68,  145,  130,  176,  161,  300,  285 },
      {  270,  255,  363,  348,  332,  317,  239,  224,  208,  193,  115,  100,   53,   38,   22,    8,   84,   69,  146,  131,  177,  162,  301,  286 },
      {  267,  252,  360,  345,  329,  314,  236,  221,  205,  190,  112,   97,   50,   35,   19,    5,   81,   66,  143,  128,  174,  159,  298,  283 },
      {  268,  253,  361,  346,  330,  315,  237,  222,  206,  191,  113,   98,   51,   36,   20,    6,   82,   67,  144,  129,  175,  160,  299,  284 },
      {  265,  250,  358,  343,  327,  312,  234,  219,  203,  188,  110,   95,   48,   33,   17,    3,   79,   64,  141,  126,  172,  157,  296,  281 },
      {  266,  251,  359,  344,  328,  313,  235,  220,  204,  189,  111,   96,   49,   34,   18,    4,   80,   65,  142,  127,  173,  158,  297,  282 },
      {  263,  248,  356,  341,  325,  310,  232,  217,  201,  186,  108,   93,   46,   31,   15,    1,   77,   62,  139,  124,  170,  155,  294,  279 },
      {  264,  249,  357,  342,  326,  311,  233,  218,  202,  187,  109,   94,   47,   32,   16,    2,   78,   63,  140,  125,  171,  156,  295,  280 }
    }
  };
  

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

//for simulation only
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

unsigned int TowerInfoDefs::encode_epd(const unsigned int arm, const unsigned int rbin, const unsigned int phibin, bool isData) {
  int idx = epdchnlmap_data[arm][rbin][phibin];

  if(!isData) {
    return TowerInfoDefs::encode_epd(arm,rbin,phibin);
  }

  return TowerInfoDefs::encode_epd(static_cast<unsigned int>(idx));

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
