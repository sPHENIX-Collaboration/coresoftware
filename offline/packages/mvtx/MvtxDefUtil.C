#include "MvtxDefUtil.h"

uint8_t
MvtxDefUtil::GetStaveId(TrkrDefs::hitsetkey key)
{
  TrkrDefs::hitsetkey tmp = (key >> kBitShiftStaveId);
  return tmp;
}

uint8_t
MvtxDefUtil::GetStaveId(TrkrDefs::cluskey key)
{
  TrkrDefs::hitsetkey tmp = (key >> kBitShiftClusId);
  return GetStaveId(tmp);
}

uint8_t
MvtxDefUtil::GetChipId(TrkrDefs::hitsetkey key)
{
  TrkrDefs::hitsetkey tmp = (key >> kBitShiftChipId);
  return tmp;
}

uint8_t
MvtxDefUtil::GetChipId(TrkrDefs::cluskey key)
{
  TrkrDefs::hitsetkey tmp = (key >> kBitShiftClusId);
  return GetChipId(tmp);
}

TrkrDefs::hitsetkey
MvtxDefUtil::GenHitSetKey(const char lyr, const uint8_t stave, const uint8_t chip)
{
  TrkrDefs::hitsetkey key = TrkrDefUtil::GenHitSetKey(TrkrDefs::TRKRID::mvtx_id, lyr);
  TrkrDefs::hitsetkey tmp = stave;
  key |= (tmp << kBitShiftStaveId);
  tmp = chip;
  key |= (tmp << kBitShiftChipId);
  return key;
}

TrkrDefs::cluskey
MvtxDefUtil::GenClusKey(const char lyr, const uint8_t stave, const uint8_t chip, const uint32_t clusid)
{
  TrkrDefs::cluskey tmp = GenHitSetKey(lyr, stave, chip);
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
