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
   * fee_id (=fiber numbers) for south side correspond to that used during test noise run on March 2 2023
   * for now, same fee_id incremented by 10, have been assigned to south side detectors, but this is arbitrary
   * all this will be updated with final TPOT cabling
   */
  // Phi layer
  {3, MicromegasDefs::genHitSetKey(55, MicromegasDefs::SegmentationType::SEGMENTATION_PHI, 0 ), "M5P",  "SCOP" },
  {8, MicromegasDefs::genHitSetKey(55, MicromegasDefs::SegmentationType::SEGMENTATION_PHI, 1 ), "M8P",  "SCIP" },
  {2, MicromegasDefs::genHitSetKey(55, MicromegasDefs::SegmentationType::SEGMENTATION_PHI, 4 ), "M9P",  "SEIP" },
  {7, MicromegasDefs::genHitSetKey(55, MicromegasDefs::SegmentationType::SEGMENTATION_PHI, 6 ), "M6P",  "SWIP" },

  {13, MicromegasDefs::genHitSetKey(55, MicromegasDefs::SegmentationType::SEGMENTATION_PHI, 3 ), "M10P", "NCOP" },
  {18, MicromegasDefs::genHitSetKey(55, MicromegasDefs::SegmentationType::SEGMENTATION_PHI, 2 ), "M4P",  "NCIP" },
  {12, MicromegasDefs::genHitSetKey(55, MicromegasDefs::SegmentationType::SEGMENTATION_PHI, 5 ), "M2P",  "NEIP" },
  {17, MicromegasDefs::genHitSetKey(55, MicromegasDefs::SegmentationType::SEGMENTATION_PHI, 7 ), "M7P",  "NWIP" },

  // z layer
  {4, MicromegasDefs::genHitSetKey(56, MicromegasDefs::SegmentationType::SEGMENTATION_Z, 0 ), "M5Z",  "SCOZ" },
  {9, MicromegasDefs::genHitSetKey(56, MicromegasDefs::SegmentationType::SEGMENTATION_Z, 1 ), "M8Z",  "SCIZ" },
  {1, MicromegasDefs::genHitSetKey(56, MicromegasDefs::SegmentationType::SEGMENTATION_Z, 4 ), "M9Z",  "SEIZ" },
  {5, MicromegasDefs::genHitSetKey(56, MicromegasDefs::SegmentationType::SEGMENTATION_Z, 6 ), "M6Z",  "SWIZ" },

  {14, MicromegasDefs::genHitSetKey(56, MicromegasDefs::SegmentationType::SEGMENTATION_Z, 3 ), "M10Z", "NCOZ" },
  {19, MicromegasDefs::genHitSetKey(56, MicromegasDefs::SegmentationType::SEGMENTATION_Z, 2 ), "M4Z",  "NCIZ" },
  {11, MicromegasDefs::genHitSetKey(56, MicromegasDefs::SegmentationType::SEGMENTATION_Z, 5 ), "M2Z",  "NEIZ" },
  {15, MicromegasDefs::genHitSetKey(56, MicromegasDefs::SegmentationType::SEGMENTATION_Z, 7 ), "M7Z",  "NWIZ" }
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

