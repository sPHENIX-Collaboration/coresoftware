#ifndef __MvtxDefUtil_H__
#define __MvtxDefUtil_H__

#include <trackbase/TrkrDefUtil.h>

class MvtxDefUtil : public TrkrDefUtil
{
 public:
  /// ctor
  MvtxDefUtil(){};

  /// dtor
  ~MvtxDefUtil(){};

  /// get the stave id
  uint8_t GetStaveId(TrkrDefs::hitsetkey key);
  uint8_t GetStaveId(TrkrDefs::cluskey key);

  /// get the chip id
  uint8_t GetChipId(TrkrDefs::hitsetkey key);
  uint8_t GetChipId(TrkrDefs::cluskey key);

  /// generate hitsetkey (This is mvtx so tracker id is already known)
  TrkrDefs::hitsetkey GenHitSetKey(const char lyr, const uint8_t stave, const uint8_t chip);

  /// generate cluskey
  TrkrDefs::cluskey GenClusKey(const char lyr, const uint8_t stave, const uint8_t chip, const uint32_t clusid);
  TrkrDefs::cluskey GenClusKey(const TrkrDefs::hitsetkey hskey, const uint32_t clusid);

 private:
  // hitsetkey layout:
  //  Mvtx specific lower 16 bits
  //   24 - 32  tracker id
  //   16 - 24  layer
  //   8  - 16  stave id
  //   0  -  8  chip id
  static const unsigned int kBitShiftStaveId = 8;
  static const unsigned int kBitShiftChipId = 0;
};

#endif  //__MvtxDefUtil_H__
