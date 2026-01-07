/*!
 * \file MicromegasDefs.cc
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "MicromegasDefs.h"

namespace
{
   //* converninece trait for underlying type
  template<class T>
    using underlying_type_t = typename std::underlying_type<T>::type;

  //* convert an strong type enum to integral type
  template<class T>
    constexpr underlying_type_t<T>
    to_underlying_type(T value) noexcept
  { return static_cast<underlying_type_t<T>>(value);}

  /*!
   * hitsetkey layout:
   * Micromegas specific lower 16 bits
   * 24 - 32  tracker id
   * 16 - 24  layer
   * 8 - 16 segmentation type
   * 0 - 8 tile id
   */
  constexpr unsigned int kBitShiftSegmentation = 8;
  constexpr unsigned int kBitShiftTileId = 0;

  //! bit shift for hit key
  constexpr unsigned int kBitShiftStrip = 0;
  constexpr unsigned int kBitShiftSample = 8;

}

namespace MicromegasDefs
{

  //________________________________________________________________
  TrkrDefs::hitsetkey genHitSetKey(uint8_t layer, SegmentationType type, uint8_t tile )
  {
    TrkrDefs::hitsetkey key = TrkrDefs::genHitSetKey(TrkrDefs::TrkrId::micromegasId, layer);

    TrkrDefs::hitsetkey tmp = to_underlying_type(type)&0x1U;
    key |= (tmp << kBitShiftSegmentation);

    tmp = tile&0xFFU;
    key |= (tmp << kBitShiftTileId);

    return key;
  }

  //________________________________________________________________
  SegmentationType getSegmentationType(TrkrDefs::hitsetkey key)
  {
    TrkrDefs::hitsetkey tmp = (key >> kBitShiftSegmentation);
    return static_cast<SegmentationType>(tmp&0x1U);
  }

  //________________________________________________________________
  uint8_t getTileId(TrkrDefs::hitsetkey key)
  {
    TrkrDefs::hitsetkey tmp = (key >> kBitShiftTileId);
    return tmp&0xFFU;
  }

  //________________________________________________________________
  TrkrDefs::hitkey genHitKey(uint16_t strip, uint16_t sample)
  {
    const TrkrDefs::hitkey key = (strip&0xFFU) << kBitShiftStrip;
    const TrkrDefs::hitkey tmp = (sample&0xFFFFU) << kBitShiftSample;
    return key|tmp;
  }

  //________________________________________________________________
  uint8_t getStrip( TrkrDefs::hitkey key )
  {
    TrkrDefs::hitkey tmp = (key >> kBitShiftStrip);
    return tmp & 0xFFU;
  }

  //________________________________________________________________
  uint16_t getSample( TrkrDefs::hitkey key )
  {
    TrkrDefs::hitkey tmp = (key >> kBitShiftSample);
    return tmp & 0xFFFFU;
  }

  //________________________________________________________________
  SegmentationType getSegmentationType(TrkrDefs::cluskey key)
  {
    const TrkrDefs::hitsetkey tmp = TrkrDefs::getHitSetKeyFromClusKey(key);
    return getSegmentationType( tmp );
  }

  //________________________________________________________________
  uint8_t getTileId(TrkrDefs::cluskey key)
  {
    const TrkrDefs::hitsetkey tmp = TrkrDefs::getHitSetKeyFromClusKey(key);
    return getTileId( tmp );
  }

}
