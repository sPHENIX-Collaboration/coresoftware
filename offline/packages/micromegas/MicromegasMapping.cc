/*!
 * \file MicromegasMapping.cc
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "MicromegasMapping.h"
#include "MicromegasDefs.h"

#include <algorithm>
#include <cassert>

namespace
{

  class mec8_channel_id
  {
    public:

    /*
    * 0 or 1, corresponding to cable1 and cable 2 as defined by Takao in
    * https://indico.bnl.gov/event/18458/contributions/73400/attachments/46043/77969/FEE_to_MEC_map_Feb17_2023.xlsx
    */
// NOLINTNEXTLINE(misc-non-private-member-variables-in-classes)
    int m_cable_id = 0;

    /*
    * 0 or 1, corresponding to j2 or j3 in transition board drawing in
    * https://wiki.sphenix.bnl.gov/index.php/File:HDR-225938-XX.PNG.png
    * https://wiki.sphenix.bnl.gov/index.php/File:HDR-225940-XX.PNG.png
    * and as defined by Takao in
    * https://indico.bnl.gov/event/18458/contributions/73400/attachments/46043/77969/FEE_to_MEC_map_Feb17_2023.xlsx
    */
// NOLINTNEXTLINE(misc-non-private-member-variables-in-classes)
    int m_connector_id = 0;

    /*
    * 1 to 70 as defined in
    * https://wiki.sphenix.bnl.gov/index.php/File:HDR-225938-XX.PNG.png
    * https://wiki.sphenix.bnl.gov/index.php/File:HDR-225940-XX.PNG.png
    */
// NOLINTNEXTLINE(misc-non-private-member-variables-in-classes)
    int m_channel_id = 0;

    // constructor
    mec8_channel_id( int cable_id, int connector_id, int channel_id ):
      m_cable_id( cable_id ),
      m_connector_id( connector_id ),
      m_channel_id( channel_id )
    {}

  };

  // less than operator
  inline bool operator < (const mec8_channel_id& lhs, const mec8_channel_id& rhs )
  {
    if( lhs.m_cable_id != rhs.m_cable_id )
    {
      return lhs.m_cable_id < rhs.m_cable_id;
    }
    if( lhs.m_connector_id != rhs.m_connector_id )
    {
      return lhs.m_connector_id < rhs.m_connector_id;
    }
    return lhs.m_channel_id < rhs.m_channel_id;
  }

  // print mapping
  [[maybe_unused]] void print_mapping( const std::string& name, const std::array<int, MicromegasDefs::m_nchannels_fee>& array )
  {
    int count = 0;
    std::cout << "int " << name << "[" << array.size() << "] = {" << std::endl << "  ";
    for( size_t i =0; i < array.size(); ++i )
    {
      if( i > 0 )
      {
	std::cout << ", ";
      }
      if( count == 32 )
      {
        std::cout << std::endl << "  ";
        count = 0;
      }
      std::cout << array[i];
      ++count;
    }

    std::cout << std::endl << "};" << std::endl;
  }

}

//____________________________________________________________________________________________________
MicromegasMapping::MicromegasMapping():
m_detectors( {
  /*
   * see https://wiki.sphenix.bnl.gov/index.php/TPc_Outer_Tracker_(TPOT)#Fiber_mapping.2C_as_of_2023.2F3.2F22
   */
  // south side
  {5,  MicromegasDefs::genHitSetKey(55, MicromegasDefs::SegmentationType::SEGMENTATION_PHI, 4 ), "sec20.0", "R3.1", "M9P",  "SEIP" },
  {7,  MicromegasDefs::genHitSetKey(56, MicromegasDefs::SegmentationType::SEGMENTATION_Z,   4 ), "sec20.1", "R3.2", "M9Z",  "SEIZ" },
  {6,  MicromegasDefs::genHitSetKey(55, MicromegasDefs::SegmentationType::SEGMENTATION_PHI, 0 ), "sec21.0", "R3.3", "M5P",  "SCOP" },
  {8,  MicromegasDefs::genHitSetKey(56, MicromegasDefs::SegmentationType::SEGMENTATION_Z,   0 ), "sec21.1", "R3.4", "M5Z",  "SCOZ" },
  {9,  MicromegasDefs::genHitSetKey(55, MicromegasDefs::SegmentationType::SEGMENTATION_PHI, 1 ), "sec21.2", "R3.5", "M8P",  "SCIP" },
  /* {10, MicromegasDefs::genHitSetKey(56, MicromegasDefs::SegmentationType::SEGMENTATION_Z,   1 ), "sec21.3", "R3.6", "M8Z",  "SCIZ" }, */
  // updated after fiber swapping on May 23, to fix flaky initialization of the FEE
  {23, MicromegasDefs::genHitSetKey(56, MicromegasDefs::SegmentationType::SEGMENTATION_Z,   1 ), "sec21.3", "R3.9", "M8Z",  "SCIZ" },
  {24, MicromegasDefs::genHitSetKey(55, MicromegasDefs::SegmentationType::SEGMENTATION_PHI, 6 ), "sec22.0", "R3.8", "M6P",  "SWIP" },
  {25, MicromegasDefs::genHitSetKey(56, MicromegasDefs::SegmentationType::SEGMENTATION_Z,   6 ), "sec22.1", "R3.8", "M6Z",  "SWIZ" },

  {11, MicromegasDefs::genHitSetKey(55, MicromegasDefs::SegmentationType::SEGMENTATION_PHI, 5 ), "sec8.0",  "R2.1", "M2P",  "NEIP" },
  {12, MicromegasDefs::genHitSetKey(56, MicromegasDefs::SegmentationType::SEGMENTATION_Z,   5 ), "sec8.1",  "R2.2", "M2Z",  "NEIZ" },
  {19, MicromegasDefs::genHitSetKey(55, MicromegasDefs::SegmentationType::SEGMENTATION_PHI, 3 ), "sec9.0",  "R2.3", "M10P", "NCOP" },
  {18, MicromegasDefs::genHitSetKey(56, MicromegasDefs::SegmentationType::SEGMENTATION_Z,   3 ), "sec9.1",  "R2.4", "M10Z", "NCOZ" },
  {0,  MicromegasDefs::genHitSetKey(55, MicromegasDefs::SegmentationType::SEGMENTATION_PHI, 2 ), "sec9.2",  "R2.5", "M4P",  "NCIP" },
  {1,  MicromegasDefs::genHitSetKey(56, MicromegasDefs::SegmentationType::SEGMENTATION_Z,   2 ), "sec9.3",  "R2.6", "M4Z",  "NCIZ" },
  {15, MicromegasDefs::genHitSetKey(55, MicromegasDefs::SegmentationType::SEGMENTATION_PHI, 7 ), "sec10.0", "R2.7", "M7P",  "NWIP" },
  {14, MicromegasDefs::genHitSetKey(56, MicromegasDefs::SegmentationType::SEGMENTATION_Z,   7 ), "sec10.1", "R2.8", "M7Z",  "NWIZ" }
} )
{
  std::cout << "MicromegasMapping::MicromegasMapping." << std::endl;

  // sort vector based on layer/tile
  std::sort( m_detectors.begin(), m_detectors.end(), []( const DetectorId& lhs, const DetectorId& rhs )
     {
       if(TrkrDefs::getLayer( lhs.m_hitsetkey ) != TrkrDefs::getLayer( rhs.m_hitsetkey ) )
       {
	 return TrkrDefs::getLayer( lhs.m_hitsetkey ) < TrkrDefs::getLayer( rhs.m_hitsetkey );
       }
       else
       {
	 return MicromegasDefs::getTileId( lhs.m_hitsetkey ) < MicromegasDefs::getTileId( rhs.m_hitsetkey );
       }
     } );

  // fill detector map from vector
  for( const auto& detector_id:m_detectors )
  {
    m_detector_map.emplace( detector_id.m_fee_id, detector_id );
  }

  // construct channel mapping
  construct_channel_mapping();
}

//____________________________________________________________________________________________________
std::vector<int> MicromegasMapping::get_fee_id_list() const
{
  std::vector<int> out;
  std::transform( m_detectors.begin(), m_detectors.end(), std::back_inserter( out ), []( const DetectorId& det_id ){ return det_id.m_fee_id; } );
  return out;
}

//____________________________________________________________________________________________________
TrkrDefs::hitsetkey MicromegasMapping::get_hitsetkey( int fee_id ) const
{
  const auto iter = m_detector_map.find( fee_id );
  if( iter == m_detector_map.end() )
  {
    std::cout << "MicromegasMapping::get_hitsetkey - invalid fee_id: " << fee_id << std::endl;
    return 0;
  }
  else
  {
    return iter->second.m_hitsetkey;
  }
}

//____________________________________________________________________________________________________
std::string MicromegasMapping::get_detname_saclay( int fee_id ) const
{
  const auto iter = m_detector_map.find( fee_id );
  if( iter == m_detector_map.end() )
  {
    std::cout << "MicromegasMapping::get_detname_saclay - invalid fee_id: " << fee_id << std::endl;
    return std::string();
  }
  else
  {
    return iter->second.m_detname_saclay;
  }
}

//____________________________________________________________________________________________________
std::string MicromegasMapping::get_detname_sphenix( int fee_id ) const
{
  const auto iter = m_detector_map.find( fee_id );
  if( iter == m_detector_map.end() )
  {
    std::cout << "MicromegasMapping::get_detname_sphenix - invalid fee_id: " << fee_id << std::endl;
    return std::string();
  }
  else
  {
    return iter->second.m_detname_sphenix;
  }
}

//____________________________________________________________________________________________________
int MicromegasMapping::get_physical_strip( int fee_id, int channel_id) const
{
  // bound check
  if( channel_id < 0 || channel_id >= MicromegasDefs::m_nchannels_fee )
  {
    std::cout << "MicromegasMapping::get_physical_strip - invalid channel: " << channel_id << std::endl;
    return -1;
  }

  // get hitsetkey and orientation
  const auto hitsetkey = get_hitsetkey(fee_id);
  const auto segmentation_type = MicromegasDefs::getSegmentationType(hitsetkey);
  switch (segmentation_type)
  {
    case MicromegasDefs::SegmentationType::SEGMENTATION_Z:
    return m_fee_to_strip_mapping_z[channel_id];

    case MicromegasDefs::SegmentationType::SEGMENTATION_PHI:
    return m_fee_to_strip_mapping_phi[channel_id];
  }

  // never reached
  return -1;
}

//____________________________________________________________________________________________________
std::string MicromegasMapping::get_detname_saclay_from_hitsetkey( TrkrDefs::hitsetkey key ) const
{
  const auto iter = std::find_if( m_detectors.begin(), m_detectors.end(), [key](const DetectorId& detector ) { return detector.m_hitsetkey == key; } );
  if( iter == m_detectors.end() )
  {
    std::cout << "MicromegasMapping::get_detname_saclay_from_hitsetkey - invalid key: " << key << std::endl;
    return std::string();
  }
  else
  {
    return iter->m_detname_saclay;
  }
}

//____________________________________________________________________________________________________
std::string MicromegasMapping::get_detname_sphenix_from_hitsetkey( TrkrDefs::hitsetkey key ) const
{
  const auto iter = std::find_if( m_detectors.begin(), m_detectors.end(), [key](const DetectorId& detector ) { return detector.m_hitsetkey == key; } );
  if( iter == m_detectors.end() )
  {
    std::cout << "MicromegasMapping::get_detname_sphenix_from_hitsetkey - invalid key: " << key << std::endl;
    return std::string();
  }
  else
  {
    return iter->m_detname_sphenix;
  }
}

//____________________________________________________________________________________________________
int MicromegasMapping::get_fee_id_from_hitsetkey( TrkrDefs::hitsetkey key ) const
{
  const auto iter = std::find_if( m_detectors.begin(), m_detectors.end(), [key](const DetectorId& detector ) { return detector.m_hitsetkey == key; } );
  if( iter == m_detectors.end() )
  {
    std::cout << "MicromegasMapping::get_fee_id_from_hitsetkey - invalid key: " << key << std::endl;
    return -1;
  }
  else
  {
    return iter->m_fee_id;
  }
}

//____________________________________________________________________________________________________
void MicromegasMapping::construct_channel_mapping()
{
  // procedure details: https://indico.bnl.gov/event/19349/contributions/75908/attachments/47219/80098/talk.pdf

  /*
   * map channel id on FEE board (0-255) to mec8 connector and channel
   * there are 2 cables (left and right) per FEE board
   * each cable has two MEC8 connector
   * on each MEC8 connector, channel 1,2, 35, 36, 69 and 70 are connected to the ground
   * the other 64 channels are signals
   */
  // source: https://indico.bnl.gov/event/18458/contributions/73400/attachments/46043/77969/FEE_to_MEC_map_Feb17_2023.xlsx
  std::array<mec8_channel_id, MicromegasDefs::m_nchannels_fee> fee_to_mec8_mapping =
  {{
    {0,0,34}, {0,0,33}, {0,0,32}, {0,0,31}, {0,0,30}, {0,0,29}, {0,0,28}, {0,0,27},
    {0,0,26}, {0,0,25}, {0,0,24}, {0,0,23}, {0,0,22}, {0,0,21}, {0,0,20}, {0,0,19},
    {0,0,18}, {0,0,17}, {0,0,16}, {0,0,15}, {0,0,14}, {0,0,13}, {0,0,12}, {0,0,11},
    {0,0,10}, {0,0,9},  {0,0,8},  {0,0,7},  {0,0,6},  {0,0,5},  {0,0,4},  {0,0,3},

    {0,0,68}, {0,0,67}, {0,0,66}, {0,0,65}, {0,0,64}, {0,0,63}, {0,0,62}, {0,0,61},
    {0,0,60}, {0,0,59}, {0,0,58}, {0,0,57}, {0,0,56}, {0,0,55}, {0,0,54}, {0,0,53},
    {0,0,52}, {0,0,51}, {0,0,50}, {0,0,49}, {0,0,48}, {0,0,47}, {0,0,46}, {0,0,45},
    {0,0,44}, {0,0,43}, {0,0,42}, {0,0,41}, {0,0,40}, {0,0,39}, {0,0,38}, {0,0,37},

    {0,1,34}, {0,1,33}, {0,1,32}, {0,1,31}, {0,1,30}, {0,1,29}, {0,1,28}, {0,1,27},
    {0,1,26}, {0,1,25}, {0,1,24}, {0,1,23}, {0,1,22}, {0,1,21}, {0,1,20}, {0,1,19},
    {0,1,18}, {0,1,17}, {0,1,16}, {0,1,15}, {0,1,14}, {0,1,13}, {0,1,12}, {0,1,11},
    {0,1,10}, {0,1,9},  {0,1,8},  {0,1,7},  {0,1,6},  {0,1,5},  {0,1,4},  {0,1,3},

    {0,1,68}, {0,1,67}, {0,1,66}, {0,1,65}, {0,1,64}, {0,1,63}, {0,1,62}, {0,1,61},
    {0,1,60}, {0,1,59}, {0,1,58}, {0,1,57}, {0,1,56}, {0,1,55}, {0,1,54}, {0,1,53},
    {0,1,52}, {0,1,51}, {0,1,50}, {0,1,49}, {0,1,48}, {0,1,47}, {0,1,46}, {0,1,45},
    {0,1,44}, {0,1,43}, {0,1,42}, {0,1,41}, {0,1,40}, {0,1,39}, {0,1,38}, {0,1,37},

    {1,0,34}, {1,0,33}, {1,0,32}, {1,0,31}, {1,0,30}, {1,0,29}, {1,0,28}, {1,0,27},
    {1,0,26}, {1,0,25}, {1,0,24}, {1,0,23}, {1,0,22}, {1,0,21}, {1,0,20}, {1,0,19},
    {1,0,18}, {1,0,17}, {1,0,16}, {1,0,15}, {1,0,14}, {1,0,13}, {1,0,12}, {1,0,11},
    {1,0,10}, {1,0,9},  {1,0,8},  {1,0,7},  {1,0,6},  {1,0,5},  {1,0,4},  {1,0,3},

    {1,0,68}, {1,0,67}, {1,0,66}, {1,0,65}, {1,0,64}, {1,0,63}, {1,0,62}, {1,0,61},
    {1,0,60}, {1,0,59}, {1,0,58}, {1,0,57}, {1,0,56}, {1,0,55}, {1,0,54}, {1,0,53},
    {1,0,52}, {1,0,51}, {1,0,50}, {1,0,49}, {1,0,48}, {1,0,47}, {1,0,46}, {1,0,45},
    {1,0,44}, {1,0,43}, {1,0,42}, {1,0,41}, {1,0,40}, {1,0,39}, {1,0,38}, {1,0,37},

    {1,1,34}, {1,1,33}, {1,1,32}, {1,1,31}, {1,1,30}, {1,1,29}, {1,1,28}, {1,1,27},
    {1,1,26}, {1,1,25}, {1,1,24}, {1,1,23}, {1,1,22}, {1,1,21}, {1,1,20}, {1,1,19},
    {1,1,18}, {1,1,17}, {1,1,16}, {1,1,15}, {1,1,14}, {1,1,13}, {1,1,12}, {1,1,11},
    {1,1,10}, {1,1,9},  {1,1,8},  {1,1,7},  {1,1,6},  {1,1,5},  {1,1,4},  {1,1,3},

    {1,1,68}, {1,1,67}, {1,1,66}, {1,1,65}, {1,1,64}, {1,1,63}, {1,1,62}, {1,1,61},
    {1,1,60}, {1,1,59}, {1,1,58}, {1,1,57}, {1,1,56}, {1,1,55}, {1,1,54}, {1,1,53},
    {1,1,52}, {1,1,51}, {1,1,50}, {1,1,49}, {1,1,48}, {1,1,47}, {1,1,46}, {1,1,45},
    {1,1,44}, {1,1,43}, {1,1,42}, {1,1,41}, {1,1,40}, {1,1,39}, {1,1,38}, {1,1,37}
  }};

  // map mec8 channel id (1-70) to mec8 signal id on detector as defined by audrey in z views (1-64)
  /* sources:
   * https://wiki.sphenix.bnl.gov/index.php/File:HDR-225938-XX.PNG.png
   * https://wiki.sphenix.bnl.gov/index.php/File:HDR-225940-XX.PNG.png
   * https://indico.bnl.gov/event/19038/contributions/75495/attachments/47029/79750/MappingGerber.pdf
   * https://indico.bnl.gov/event/17391/contributions/68961/attachments/47061/79816/TpotProgress.pdf, slide 2
   */
  std::map<int,int> mec8_to_signal_mapping_z = {
    /* bottom row on mec8 (female) connector corresponds to bottom row on detector as defined by Audrey */
    {1,-1},
    {3,1},  {5,2},   {7,3},   {9,4},   {11,5},  {13,6},  {15,7},  {17,8},
    {19,9}, {21,10}, {23,11}, {25,12}, {27,13}, {29,14}, {31,15}, {33,16},
    {35,-1},
    {37,17}, {39,18}, {41,19}, {43,20}, {45,21}, {47,22}, {49,23}, {51,24},
    {53,25}, {55,26}, {57,27}, {59,28}, {61,29}, {63,30}, {65,31}, {67,32},
    {69,-1},
    /* top row on mec8 (female) connector corresponds to top row on detector as defined by Audrey */
    {2,-1},
    {4,33},  {6,34},  {8,35},  {10,36}, {12,37}, {14,38}, {16,39}, {18,40},
    {20,41}, {22,42}, {24,43}, {26,44}, {28,45}, {30,46}, {32,47}, {34,48},
    {36,-1},
    {38,49}, {40,50}, {42,51}, {44,52}, {46,53}, {48,54}, {50,55}, {52,56},
    {54,57}, {56,58}, {58,59}, {60,60}, {62,61}, {64,62}, {66,63}, {68,64},
    {70,-1}
  };

  // map all mec8 channel id to signal id on detector as defined by audrey in z views (1-256)
  /* sources:
   * https://indico.bnl.gov/event/19038/contributions/75495/attachments/47029/79750/MappingGerber.pdf
   * https://indico.bnl.gov/event/17391/contributions/68961/attachments/47061/79816/TpotProgress.pdf, slide 2
   */
  std::map<mec8_channel_id,int> mec8_to_signal_mapping_z_all;
  for( const auto& [mec8_channel,signal_id]:mec8_to_signal_mapping_z )
  {

    // ignore ground channels
    if( signal_id == -1 ) { continue;
}

    // cable 0, connector 0 corresponds to signal ids 1 to 64
    mec8_to_signal_mapping_z_all.insert(std::pair<mec8_channel_id,int>({0,0,mec8_channel}, signal_id));

    // cable 0, connector 1 corresponds to signal ids 65 to 128
    mec8_to_signal_mapping_z_all.insert(std::pair<mec8_channel_id,int>({0,1,mec8_channel}, signal_id+64));

    // cable 1, connector 0 corresponds to signal ids 129 to 192
    mec8_to_signal_mapping_z_all.insert(std::pair<mec8_channel_id,int>({1,0,mec8_channel}, signal_id+128));

    // cable 1, connector 1 corresponds to signal ids 193 to 256
    mec8_to_signal_mapping_z_all.insert(std::pair<mec8_channel_id,int>({1,1,mec8_channel}, signal_id+192));
  }

  // map phisical strips in the detector to MEC8 signal id on detector as defined by Audrey
  // source: https://indico.bnl.gov/event/19038/contributions/75495/attachments/47029/79824/MappingTPOT.xlsx
  std::map<int,int> strip_to_signal_id_mapping_all =
  {
    {1,1}, {2,33}, {3,2}, {4,34}, {5,3}, {6,35}, {7,4}, {8,36},
    {9,5}, {10,37}, {11,6}, {12,38}, {13,7}, {14,39}, {15,8}, {16,40},
    {17,9}, {18,41}, {19,10}, {20,42}, {21,11}, {22,43}, {23,12}, {24,44},
    {25,13}, {26,45}, {27,14}, {28,46}, {29,15}, {30,47}, {31,16}, {32,48},
    {33,17}, {34,49}, {35,18}, {36,50}, {37,19}, {38,51}, {39,20}, {40,52},
    {41,21}, {42,53}, {43,22}, {44,54}, {45,23}, {46,55}, {47,24}, {48,56},
    {49,25}, {50,57}, {51,26}, {52,58}, {53,27}, {54,59}, {55,28}, {56,60},
    {57,29}, {58,61}, {59,30}, {60,62}, {61,31}, {62,63}, {63,32}, {64,64},

    {65,65}, {66,97}, {67,66}, {68,98}, {69,67}, {70,99}, {71,68}, {72,100},
    {73,69}, {74,101}, {75,70}, {76,102}, {77,71}, {78,103}, {79,72}, {80,104},
    {81,73}, {82,105}, {83,74}, {84,106}, {85,75}, {86,107}, {87,76}, {88,108},
    {89,77}, {90,109}, {91,78}, {92,110}, {93,79}, {94,111}, {95,80}, {96,112},
    {97,81}, {98,113}, {99,82}, {100,114}, {101,83}, {102,115}, {103,84}, {104,116},
    {105,85}, {106,117}, {107,86}, {108,118}, {109,87}, {110,119}, {111,88}, {112,120},
    {113,89}, {114,121}, {115,90}, {116,122}, {117,91}, {118,123}, {119,92}, {120,124},
    {121,93}, {122,125}, {123,94}, {124,126}, {125,95}, {126,127}, {127,96}, {128,128},

    {129,129}, {130,161}, {131,130}, {132,162}, {133,131}, {134,163}, {135,132}, {136,164},
    {137,133}, {138,165}, {139,134}, {140,166}, {141,135}, {142,167}, {143,136}, {144,168},
    {145,137}, {146,169}, {147,138}, {148,170}, {149,139}, {150,171}, {151,140}, {152,172},
    {153,141}, {154,173}, {155,142}, {156,174}, {157,143}, {158,175}, {159,144}, {160,176},
    {161,145}, {162,177}, {163,146}, {164,178}, {165,147}, {166,179}, {167,148}, {168,180},
    {169,149}, {170,181}, {171,150}, {172,182}, {173,151}, {174,183}, {175,152}, {176,184},
    {177,153}, {178,185}, {179,154}, {180,186}, {181,155}, {182,187}, {183,156}, {184,188},
    {185,157}, {186,189}, {187,158}, {188,190}, {189,159}, {190,191}, {191,160}, {192,192},

    {193,193}, {194,225}, {195,194}, {196,226}, {197,195}, {198,227}, {199,196}, {200,228},
    {201,197}, {202,229}, {203,198}, {204,230}, {205,199}, {206,231}, {207,200}, {208,232},
    {209,201}, {210,233}, {211,202}, {212,234}, {213,203}, {214,235}, {215,204}, {216,236},
    {217,205}, {218,237}, {219,206}, {220,238}, {221,207}, {222,239}, {223,208}, {224,240},
    {225,209}, {226,241}, {227,210}, {228,242}, {229,211}, {230,243}, {231,212}, {232,244},
    {233,213}, {234,245}, {235,214}, {236,246}, {237,215}, {238,247}, {239,216}, {240,248},
    {241,217}, {242,249}, {243,218}, {244,250}, {245,219}, {246,251}, {247,220}, {248,252},
    {249,221}, {250,253}, {251,222}, {252,254}, {253,223}, {254,255}, {255,224}, {256,256}
  };

  // add mapping from strip number as defined by Audrey and geant convention
  auto get_strip_geant_z = []( int strip_audrey )
  {
    /*
     * for z views, audrey and geant strips are numbered in oposite directions
     * geant strips start from zero
    */
    return MicromegasDefs::m_nchannels_fee-strip_audrey;
  };

  // construct fee channel id to strip
  /*
   * fee channel id is from 0 to 255
   * strip id from 0 to 255 (-1 with respect to Audrey's convention
   */
  for( int channel_id = 0; channel_id < MicromegasDefs::m_nchannels_fee; ++channel_id )
  {

    // mec8 channel id
    auto mec8_channel_id = fee_to_mec8_mapping[channel_id];

    // mec8 signal id as defined by Audrey
    int mec8_signal_id = mec8_to_signal_mapping_z_all.at(mec8_channel_id);

    // find mec8_signal_id in detector mapping
    const auto iter = std::find_if( strip_to_signal_id_mapping_all.begin(), strip_to_signal_id_mapping_all.end(), [mec8_signal_id]( const std::pair<int,int>& pair ) { return pair.second == mec8_signal_id; } );
    assert( iter != strip_to_signal_id_mapping_all.end() );
    const int strip_audrey = iter->first;

    // convert to geant convention
    const int strip_geant = get_strip_geant_z( strip_audrey );

    // store in array
    m_fee_to_strip_mapping_z[channel_id] = strip_geant;
  }

  // print_mapping( "m_fee_to_strip_mapping_z", m_fee_to_strip_mapping_z );

  // map mec8 channel id (1-70) to mec8 signal id on detector as defined by audrey in phi views
  /* sources:
   * https://wiki.sphenix.bnl.gov/index.php/File:HDR-225938-XX.PNG.png
   * https://wiki.sphenix.bnl.gov/index.php/File:HDR-225940-XX.PNG.png
   * https://indico.bnl.gov/event/19038/contributions/75495/attachments/47029/79823/MappingGerber.pdf
   * https://indico.bnl.gov/event/17391/contributions/68961/attachments/47061/79816/TpotProgress.pdf, slide 2
   */
  std::map<int,int> mec8_to_signal_mapping_phi = {
    /* bottom row on mec8 (female) connector corresponds to top row on detector as defined by Audrey */
    {1,-1},
    {3,64},  {5,63},  {7,62},  {9,61},  {11,60}, {13,59}, {15,58}, {17,57},
    {19,56}, {21,55}, {23,54}, {25,53}, {27,52}, {29,51}, {31,50}, {33,49},
    {35,-1},
    {37,48}, {39,47}, {41,46}, {43,45}, {45,44}, {47,43}, {49,42}, {51,41},
    {53,40}, {55,39}, {57,38}, {59,37}, {61,36}, {63,35}, {65,34}, {67,33},
    {69,-1},
    /* top row on mec8 (female) connector corresponds to bottom row on detector as defined by Audrey */
    {2,-1},
    {4,32},  {6,31},  {8,30},  {10,29}, {12,28}, {14,27}, {16,26}, {18,25},
    {20,24}, {22,23}, {24,22}, {26,21}, {28,20}, {30,19}, {32,18}, {34,17},
    {36,-1},
    {38,16}, {40,15}, {42,14}, {44,13}, {46,12}, {48,11}, {50,10}, {52,9},
    {54,8},  {56,7},  {58,6},  {60,5},  {62,4},  {64,3},  {66,2},  {68,1},
    {70,-1}
  };

  // map all mec8 channel id to signal id on detector as defined by audrey in phi views
  /* sources:
   * https://indico.bnl.gov/event/19038/contributions/75495/attachments/47029/79823/MappingGerber.pdf
   * https://indico.bnl.gov/event/17391/contributions/68961/attachments/47061/79816/TpotProgress.pdf, slide 2
   */
  std::map<mec8_channel_id,int> mec8_to_signal_mapping_phi_all;
  for( const auto& [mec8_channel,signal_id]:mec8_to_signal_mapping_phi )
  {

    // ignore ground channels
    if( signal_id == -1 ) { continue;
}

    // cable 0, connector 0 corresponds to signal ids 193 to 256
    mec8_to_signal_mapping_phi_all.insert(std::pair<mec8_channel_id,int>({0,0,mec8_channel}, signal_id+192));

    // cable 0, connector 1 corresponds to signal ids 129 to 160
    mec8_to_signal_mapping_phi_all.insert(std::pair<mec8_channel_id,int>({0,1,mec8_channel}, signal_id+128));

    // cable 1, connector 0 corresponds to signal ids 65 to 128
    mec8_to_signal_mapping_phi_all.insert(std::pair<mec8_channel_id,int>({1,0,mec8_channel}, signal_id+64));

    // cable 1, connector 1 corresponds to signal ids 1 to 64
    mec8_to_signal_mapping_phi_all.insert(std::pair<mec8_channel_id,int>({1,1,mec8_channel}, signal_id));
  }

  // add mapping from strip number as defined by Audrey and geant convention
  auto get_strip_geant_phi = []( int strip_audrey )
  {
    /*
     * for phi views, audrey and geant strips are numbered in oposite directions
     * geant strips start from zero
    */
    return MicromegasDefs::m_nchannels_fee-strip_audrey;
  };

  // construct fee channel id to strip
  /*
   * fee channel id is from 0 to 255
   * strip id from 0 to 255 (-1 with respect to Audrey's convention
   */
  for( int channel_id = 0; channel_id < MicromegasDefs::m_nchannels_fee; ++channel_id )
  {

    // mec8 channel id
    auto mec8_channel_id = fee_to_mec8_mapping[channel_id];

    // mec8 signal id as defined by Audrey
    int mec8_signal_id = mec8_to_signal_mapping_phi_all.at(mec8_channel_id);

    // find mec8_signal_id in detector mapping
    const auto iter = std::find_if( strip_to_signal_id_mapping_all.begin(), strip_to_signal_id_mapping_all.end(), [mec8_signal_id]( const std::pair<int,int>& pair ) { return pair.second == mec8_signal_id; } );
    assert( iter != strip_to_signal_id_mapping_all.end() );
    const int strip_audrey = iter->first;

    // convert to geant convention
    const int strip_geant = get_strip_geant_phi( strip_audrey );

    // store in array
    m_fee_to_strip_mapping_phi[channel_id] = strip_geant;
  }

  // print_mapping( "m_fee_to_strip_mapping_phi", m_fee_to_strip_mapping_phi );

}
