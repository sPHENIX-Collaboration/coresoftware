/*!
 * \file MicromegasMapping.cc
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "MicromegasMapping.h"
#include "MicromegasDefs.h"

#include <algorithm>

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
  {10, MicromegasDefs::genHitSetKey(56, MicromegasDefs::SegmentationType::SEGMENTATION_Z,   1 ), "sec21.3", "R3.6", "M8Z",  "SCIZ" },
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

  /// fill detector map from vector
  for( const auto& detector_id:m_detectors )
  { m_detector_map.emplace( detector_id.m_fee_id, detector_id ); }
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
  } else return iter->second.m_hitsetkey;
}

//____________________________________________________________________________________________________
std::string MicromegasMapping::get_detname_saclay( int fee_id ) const
{
  const auto iter = m_detector_map.find( fee_id );
  if( iter == m_detector_map.end() )
  {
    std::cout << "MicromegasMapping::get_detname_saclay - invalid fee_id: " << fee_id << std::endl;
    return std::string();
  } else return iter->second.m_detname_saclay;
}

//____________________________________________________________________________________________________
std::string MicromegasMapping::get_detname_sphenix( int fee_id ) const
{
  const auto iter = m_detector_map.find( fee_id );
  if( iter == m_detector_map.end() )
  {
    std::cout << "MicromegasMapping::get_detname_sphenix - invalid fee_id: " << fee_id << std::endl;
    return std::string();
  } else return iter->second.m_detname_sphenix;
}

//____________________________________________________________________________________________________
int MicromegasMapping::get_physical_strip( int /*fee_id*/, int channel_id) const
{
  /*
   * this maps channel id (0-255) on a given FEE board, to physical strip number in the detector
   * Just a placeholder for now. We just return the same index.
   */
  return channel_id;
}

