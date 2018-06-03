#include "MvtxDefUtil.h"

uint8_t
MvtxDefUtil::getStaveId(TrkrDefs::hitsetkey key)
{
  TrkrDefs::hitsetkey tmp = (key >> kBitShiftStaveId);
  return tmp;
}

uint8_t
MvtxDefUtil::getStaveId(TrkrDefs::cluskey key)
{
  TrkrDefs::hitsetkey tmp = (key >> kBitShiftClusId);
  return getStaveId(tmp);
}

uint8_t
MvtxDefUtil::getChipId(TrkrDefs::hitsetkey key)
{
  TrkrDefs::hitsetkey tmp = (key >> kBitShiftChipId);
  return tmp;
}

uint8_t
MvtxDefUtil::getChipId(TrkrDefs::cluskey key)
{
  TrkrDefs::hitsetkey tmp = (key >> kBitShiftClusId);
  return getChipId(tmp);
}

TrkrDefs::hitsetkey
MvtxDefUtil::genHitSetKey(const char lyr, const uint8_t stave, const uint8_t chip)
{
  TrkrDefs::hitsetkey key = TrkrDefUtil::genHitSetKey(TrkrDefs::TrkrId::mvtxId, lyr);
  TrkrDefs::hitsetkey tmp = stave;
  key |= (tmp << kBitShiftStaveId);
  tmp = chip;
  key |= (tmp << kBitShiftChipId);
  return key;
}

TrkrDefs::cluskey
MvtxDefUtil::genClusKey(const char lyr, const uint8_t stave, const uint8_t chip, const uint32_t clusid)
{
  TrkrDefs::cluskey tmp = genHitSetKey(lyr, stave, chip);
  TrkrDefs::cluskey key = (tmp << kBitShiftClusId);
  key |= clusid;
  return key;
}

TrkrDefs::cluskey
MvtxDefUtil::GenClusKey(const TrkrDefs::hitsetkey hskey, const uint32_t clusid)
{
  TrkrDefs::cluskey tmp = hskey;
  TrkrDefs::cluskey key = (tmp << kBitShiftClusId);
  key |= clusid;
  return key;
}
