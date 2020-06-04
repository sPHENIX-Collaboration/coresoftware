/*!
 * \file MicromegasDefs.cc
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "MicromegasDefs.h"


namespace MicromegasDefs
{

  //________________________________________________________________
  uint8_t getTileId(TrkrDefs::hitsetkey key)
  {
    TrkrDefs::hitsetkey tmp = (key >> kBitShiftTileId);
    return tmp;
  }

  //________________________________________________________________
  TrkrDefs::hitsetkey genHitSetKey(uint8_t layer, uint8_t tile )
  {
    TrkrDefs::hitsetkey key = TrkrDefs::genHitSetKey(TrkrDefs::TrkrId::micromegasId, layer);
    TrkrDefs::hitsetkey tmp = tile;
    key |= (tmp << kBitShiftTileId);
    return key;
  }

    //________________________________________________________________
  TrkrDefs::hitkey genHitKey(uint16_t strip)
  {
    TrkrDefs::hitkey key = strip << kBitShiftStrip;
    return key;
  }

}
