#ifndef MVTX_MVTXDEFUTIL_H 
#define MVTX_MVTXDEFUTIL_H 

#include <trackbase/TrkrDefs.h>

namespace MvtxDefs
{

#ifndef __CINT__
  // hitsetkey layout:
  //  Mvtx specific lower 16 bits
  //   24 - 32  tracker id
  //   16 - 24  layer
  //   8  - 16  stave id
  //   0  -  8  chip id
  static const unsigned int kBitShiftStaveId __attribute__((unused)) = 8;
  static const unsigned int kBitShiftChipId __attribute__((unused)) = 0;

  // max values for col and row index in chip
  static const uint16_t MAXCOL __attribute__((unused)) = 1024;
  static const uint16_t MAXROW __attribute__((unused)) = 512;

#endif // __CINT__

  /// get the stave id
  uint8_t getStaveId(TrkrDefs::hitsetkey key);
  uint8_t getStaveId(TrkrDefs::cluskey key);

  /// get the chip id
  uint8_t getChipId(TrkrDefs::hitsetkey key);
  uint8_t getChipId(TrkrDefs::cluskey key);

  /// generate hitsetkey (This is mvtx so tracker id is already known)
  TrkrDefs::hitsetkey genHitSetKey(const char lyr, const uint8_t stave, const uint8_t chip);

  /// generate cluskey
  TrkrDefs::cluskey genClusKey(const char lyr, const uint8_t stave, const uint8_t chip, const uint32_t clusid);
  TrkrDefs::cluskey genClusKey(const TrkrDefs::hitsetkey hskey, const uint32_t clusid);


}

#endif  //MVTX_MVTXDEFUTIL_H 
