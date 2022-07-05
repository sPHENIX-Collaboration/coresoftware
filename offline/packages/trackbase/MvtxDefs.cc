#include "MvtxDefs.h"

uint8_t
MvtxDefs::getStaveId(TrkrDefs::hitsetkey key)
{
  TrkrDefs::hitsetkey tmp = (key >> MvtxDefs::kBitShiftStaveIdOffset);
  // zero the bits not in the stave id field
  uint8_t tmp1 = tmp;
  tmp1 = (tmp1 <<  (8 - MvtxDefs::kBitShiftStaveIdWidth));
  tmp1 = (tmp1 >>  (8 - MvtxDefs::kBitShiftStaveIdWidth));
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
  TrkrDefs::hitsetkey tmp = (key >> MvtxDefs::kBitShiftChipIdOffset);
  uint8_t tmp1 = tmp;
  tmp1 = (tmp1 <<  (8 - MvtxDefs::kBitShiftChipIdWidth));
  tmp1 = (tmp1 >>  (8 - MvtxDefs::kBitShiftChipIdWidth));
  return tmp1;
}

uint8_t
MvtxDefs::getChipId(TrkrDefs::cluskey key)
{
  TrkrDefs::hitsetkey tmp = (key >> TrkrDefs::kBitShiftClusId);
  return getChipId(tmp);
}

int
MvtxDefs::getStrobeId(TrkrDefs::hitsetkey key)
{
  TrkrDefs::hitsetkey tmp = (key >> MvtxDefs::kBitShiftStrobeIdOffset);
  uint8_t tmp1 = tmp;
  tmp1 = (tmp1 <<  (8 - MvtxDefs::kBitShiftStrobeIdWidth));
  tmp1 = (tmp1 >>  (8 - MvtxDefs::kBitShiftStrobeIdWidth));

  int tmp2 = (int) tmp1 - strobeOffset; // get back to the signed strobe

  return tmp2;
}

int
MvtxDefs::getStrobeId(TrkrDefs::cluskey key)
{
  TrkrDefs::hitsetkey tmp = (key >> TrkrDefs::kBitShiftClusId);
  return getStrobeId(tmp);
}

uint16_t
MvtxDefs::getCol(TrkrDefs::hitkey key)
{
  TrkrDefs::hitkey tmp = (key >> MvtxDefs::kBitShiftCol);
  return tmp;
}

uint16_t
MvtxDefs::getRow(TrkrDefs::hitkey key)
{
  TrkrDefs::hitkey tmp = (key >> MvtxDefs::kBitShiftRow);
  return tmp;
}

TrkrDefs::hitkey
MvtxDefs::genHitKey(const uint16_t col, const uint16_t row)
{
  TrkrDefs::hitkey key = (col << MvtxDefs::kBitShiftCol);
  TrkrDefs::hitkey tmp = (row << MvtxDefs::kBitShiftRow);
  key |= tmp;
  return key;
}

TrkrDefs::hitsetkey
MvtxDefs::genHitSetKey(const uint8_t lyr, const uint8_t stave, const uint8_t chip, const int strobe_in)
{
  TrkrDefs::hitsetkey key = TrkrDefs::genHitSetKey(TrkrDefs::TrkrId::mvtxId, lyr);

 // offset strobe to make it positive, fit inside 5 bits
  int strobe = strobe_in + strobeOffset;
  if(strobe < 0) strobe = 0;  
  if(strobe > 31) strobe = 31;
  unsigned int ustrobe = (unsigned int) strobe;

  TrkrDefs::hitsetkey tmp = stave;
  key |= (tmp << MvtxDefs::kBitShiftStaveIdOffset);
  tmp = chip;
  key |= (tmp << MvtxDefs::kBitShiftChipIdOffset);
  tmp = ustrobe;
  key |= (tmp << MvtxDefs::kBitShiftStrobeIdOffset);
  return key;
}

TrkrDefs::cluskey 
MvtxDefs::genClusKey(const uint8_t lyr, const uint8_t stave, const uint8_t chip, const int strobe, const uint32_t clusid)
{
  TrkrDefs::hitsetkey key = genHitSetKey(lyr, stave, chip, strobe);
  //return TrkrDefs::genClusKey( key, clusid );
  return TrkrDefs::genClusKey(key,clusid);
}

TrkrDefs::hitsetkey
MvtxDefs::resetStrobeHitSetKey(const TrkrDefs::hitsetkey hitsetkey)
{
  // Note: this method uses the fact that the crossing is in the first 5 bits
   TrkrDefs::hitsetkey tmp = hitsetkey;
   // zero the crossing bits by shifting them out of the word, then shift back
   tmp = (tmp >>  MvtxDefs::kBitShiftStrobeIdWidth);
   tmp = (tmp << MvtxDefs::kBitShiftStrobeIdWidth);
   unsigned int zero_strobe = strobeOffset;
   tmp |= (zero_strobe <<  MvtxDefs::kBitShiftStrobeIdOffset);

  return tmp;
}
