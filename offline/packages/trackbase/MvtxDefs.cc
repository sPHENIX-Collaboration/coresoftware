#include "MvtxDefs.h"

namespace
{
  // hitsetkey layout:
  //  Mvtx specific lower 16 bits
  //   24 - 31  tracker id  // 8 bits
  //   16 - 23  layer   // 8 bits
  //     9 - 15  stave id // 7 bits
  //     5 - 8     chip  id // 4 bits
  //     0-4   strobe  id // 5 bits

  static constexpr unsigned int kBitShiftStaveIdOffset = 9;
  static constexpr unsigned int kBitShiftStaveIdWidth = 7;
  static constexpr unsigned int kBitShiftChipIdOffset = 5;
  static constexpr unsigned int kBitShiftChipIdWidth = 4;
  static constexpr unsigned int kBitShiftStrobeIdOffset = 0;
  static constexpr unsigned int kBitShiftStrobeIdWidth = 5;
  static constexpr int strobeOffset = 16;

  // bit shift for hitkey
  static const unsigned int kBitShiftCol __attribute__((unused)) = 16;
  static const unsigned int kBitShiftRow __attribute__((unused)) = 0;
}

uint8_t
MvtxDefs::getStaveId(TrkrDefs::hitsetkey key)
{
  TrkrDefs::hitsetkey tmp = (key >> kBitShiftStaveIdOffset);
  // zero the bits not in the stave id field
  uint8_t tmp1 = tmp;
  tmp1 = (tmp1 << (8 - kBitShiftStaveIdWidth));
  tmp1 = (tmp1 >> (8 - kBitShiftStaveIdWidth));
  return tmp1;
}

uint8_t
MvtxDefs::getStaveId(TrkrDefs::cluskey key)
{
  TrkrDefs::hitsetkey tmp = (key >> TrkrDefs::kBitShiftClusId);
  return getStaveId(tmp);
}

uint8_t
MvtxDefs::getChipId(TrkrDefs::hitsetkey key)
{
  TrkrDefs::hitsetkey tmp = (key >> kBitShiftChipIdOffset);
  uint8_t tmp1 = tmp;
  tmp1 = (tmp1 << (8 - kBitShiftChipIdWidth));
  tmp1 = (tmp1 >> (8 - kBitShiftChipIdWidth));
  return tmp1;
}

uint8_t
MvtxDefs::getChipId(TrkrDefs::cluskey key)
{
  TrkrDefs::hitsetkey tmp = (key >> TrkrDefs::kBitShiftClusId);
  return getChipId(tmp);
}

int MvtxDefs::getStrobeId(TrkrDefs::hitsetkey key)
{
  TrkrDefs::hitsetkey tmp = (key >> kBitShiftStrobeIdOffset);
  uint8_t tmp1 = tmp;
  tmp1 = (tmp1 << (8 - kBitShiftStrobeIdWidth));
  tmp1 = (tmp1 >> (8 - kBitShiftStrobeIdWidth));

  int tmp2 = (int) tmp1 - strobeOffset;  // get back to the signed strobe

  return tmp2;
}

int MvtxDefs::getStrobeId(TrkrDefs::cluskey key)
{
  TrkrDefs::hitsetkey tmp = (key >> TrkrDefs::kBitShiftClusId);
  return getStrobeId(tmp);
}

uint16_t
MvtxDefs::getCol(TrkrDefs::hitkey key)
{
  TrkrDefs::hitkey tmp = (key >> kBitShiftCol);
  return tmp;
}

uint16_t
MvtxDefs::getRow(TrkrDefs::hitkey key)
{
  TrkrDefs::hitkey tmp = (key >> kBitShiftRow);
  return tmp;
}

TrkrDefs::hitkey
MvtxDefs::genHitKey(const uint16_t col, const uint16_t row)
{
  TrkrDefs::hitkey key = (col << kBitShiftCol);
  TrkrDefs::hitkey tmp = (row << kBitShiftRow);
  key |= tmp;
  return key;
}

TrkrDefs::hitsetkey
MvtxDefs::genHitSetKey(const uint8_t lyr, const uint8_t stave, const uint8_t chip, const int strobe_in)
{
  TrkrDefs::hitsetkey key = TrkrDefs::genHitSetKey(TrkrDefs::TrkrId::mvtxId, lyr);

  // offset strobe to make it positive, fit inside 5 bits
  int strobe = strobe_in + strobeOffset;
  if (strobe < 0)
  {
    strobe = 0;
  }
  if (strobe > 31)
  {
    strobe = 31;
  }
  unsigned int ustrobe = (unsigned int) strobe;

  TrkrDefs::hitsetkey tmp = stave;
  key |= (tmp << kBitShiftStaveIdOffset);
  tmp = chip;
  key |= (tmp << kBitShiftChipIdOffset);
  tmp = ustrobe;
  key |= (tmp << kBitShiftStrobeIdOffset);
  return key;
}

TrkrDefs::cluskey
MvtxDefs::genClusKey(const uint8_t lyr, const uint8_t stave, const uint8_t chip, const int strobe, const uint32_t clusid)
{
  TrkrDefs::hitsetkey key = genHitSetKey(lyr, stave, chip, strobe);
  // return TrkrDefs::genClusKey( key, clusid );
  return TrkrDefs::genClusKey(key, clusid);
}

TrkrDefs::hitsetkey
MvtxDefs::resetStrobeHitSetKey(const TrkrDefs::hitsetkey hitsetkey)
{
  // Note: this method uses the fact that the crossing is in the first 5 bits
  TrkrDefs::hitsetkey tmp = hitsetkey;
  // zero the crossing bits by shifting them out of the word, then shift back
  tmp = (tmp >> kBitShiftStrobeIdWidth);
  tmp = (tmp << kBitShiftStrobeIdWidth);
  unsigned int zero_strobe = strobeOffset;
  tmp |= (zero_strobe << kBitShiftStrobeIdOffset);

  return tmp;
}
