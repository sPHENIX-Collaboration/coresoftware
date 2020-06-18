/*!
 * \file MicromegasDefs.cc
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "MicromegasDefs.h"


namespace MicromegasDefs
{

  //________________________________________________________________
  TrkrDefs::hitsetkey genHitSetKey(uint8_t layer, uint8_t tile )
  {
    TrkrDefs::hitsetkey key = TrkrDefs::genHitSetKey(TrkrDefs::TrkrId::micromegasId, layer);
    TrkrDefs::hitsetkey tmp = tile;
    key |= (tmp << kBitShiftTileId);
    return key;
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

}
