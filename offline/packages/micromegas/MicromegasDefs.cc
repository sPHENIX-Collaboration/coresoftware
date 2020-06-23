/*!
 * \file MicromegasDefs.cc
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "MicromegasDefs.h"


namespace MicromegasDefs
{

  //________________________________________________________________
  TrkrDefs::hitsetkey genHitSetKey(uint8_t layer, SegmentationType type, uint8_t tile )
  {
    TrkrDefs::hitsetkey key = TrkrDefs::genHitSetKey(TrkrDefs::TrkrId::micromegasId, layer);
    
    TrkrDefs::hitsetkey tmp = type;
    key |= (tmp << kBitShiftSegmentation);

    tmp = tile;
    key |= (tmp << kBitShiftTileId);
    
//     std::cout << "MicromegasDefs::genHitSetKey -"
//       << " layer: " << int(layer)
//       << " segmentation: " << type 
//       << " tile: " << int(tile)
//       << std::endl;
    
    return key;
  }

  //________________________________________________________________
  SegmentationType getSegmentationType(TrkrDefs::hitsetkey key)
  {
    TrkrDefs::hitsetkey tmp = (key >> kBitShiftSegmentation);
    return (SegmentationType)tmp;
  }

  //________________________________________________________________
  uint8_t getTileId(TrkrDefs::hitsetkey key)
  {
    TrkrDefs::hitsetkey tmp = (key >> kBitShiftTileId);
    return tmp;
  }

  //________________________________________________________________
  TrkrDefs::hitkey genHitKey(uint16_t strip)
  {
    TrkrDefs::hitkey key = strip << kBitShiftStrip;
    return key;
  }

  //________________________________________________________________
  uint16_t getStrip( TrkrDefs::hitkey key )
  {
    TrkrDefs::hitkey tmp = (key >> kBitShiftStrip);
    return tmp;
  }

  //________________________________________________________________
  TrkrDefs::cluskey genClusterKey(TrkrDefs::hitsetkey hitsetkey, uint32_t clusid)
  {
    TrkrDefs::cluskey tmp = hitsetkey;
    TrkrDefs::cluskey key = (tmp << TrkrDefs::kBitShiftClusId);
    key |= clusid;
    return key;
  }

  //________________________________________________________________
  SegmentationType getSegmentationType(TrkrDefs::cluskey key)
  {
    TrkrDefs::hitsetkey tmp = (key >> TrkrDefs::kBitShiftClusId);
    return getSegmentationType( tmp );
  }

  //________________________________________________________________
  uint8_t getTileId(TrkrDefs::cluskey key)
  {
    TrkrDefs::hitsetkey tmp = (key >> TrkrDefs::kBitShiftClusId);
    return getTileId( tmp );
  }

}
