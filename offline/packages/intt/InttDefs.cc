/**
 * @file intt/InttDefs.cc
 * @author D. McGlinchey
 * @date June 2018
 * @brief Implementation file for InttDefs.h
 */
#include "InttDefs.h"

uint8_t
InttDefs::getStaveId(TrkrDefs::hitsetkey key)
{
  TrkrDefs::hitsetkey tmp = (key >> InttDefs::kBitShiftStaveId);
  return tmp;
}

uint8_t
InttDefs::getStaveId(TrkrDefs::cluskey key)
{
  TrkrDefs::hitsetkey tmp = (key >> TrkrDefs::kBitShiftClusId);
  return getStaveId(tmp);
}

uint8_t
InttDefs::getChipId(TrkrDefs::hitsetkey key)
{
  TrkrDefs::hitsetkey tmp = (key >> InttDefs::kBitShiftChipId);
  return tmp;
}

uint8_t
InttDefs::getChipId(TrkrDefs::cluskey key)
{
  TrkrDefs::hitsetkey tmp = (key >> TrkrDefs::kBitShiftClusId);
  return getChipId(tmp);
}

uint16_t
InttDefs::getCol(TrkrDefs::hitkey key)
{
  TrkrDefs::hitkey tmp = (key >> InttDefs::kBitShiftCol);
  return tmp;
}

uint16_t
InttDefs::getRow(TrkrDefs::hitkey key)
{
  TrkrDefs::hitkey tmp = (key >> InttDefs::kBitShiftRow);
  return tmp;
}

TrkrDefs::hitkey 
InttDefs::genHitKey(const uint16_t col, const uint16_t row)
{
  TrkrDefs::hitkey key = (col << InttDefs::kBitShiftCol);
  TrkrDefs::hitkey tmp = (row << InttDefs::kBitShiftRow);
  key |= tmp;
  return key;
}

TrkrDefs::hitsetkey
InttDefs::genHitSetKey(const char lyr, const uint8_t stave, const uint8_t chip)
{
  TrkrDefs::hitsetkey key = TrkrDefs::genHitSetKey(TrkrDefs::TrkrId::mvtxId, lyr);
  TrkrDefs::hitsetkey tmp = stave;
  key |= (tmp << InttDefs::kBitShiftStaveId);
  tmp = chip;
  key |= (tmp << InttDefs::kBitShiftChipId);
  return key;
}

TrkrDefs::cluskey
InttDefs::genClusKey(const char lyr, const uint8_t stave, const uint8_t chip, const uint32_t clusid)
{
  TrkrDefs::cluskey tmp = genHitSetKey(lyr, stave, chip);
  TrkrDefs::cluskey key = (tmp << TrkrDefs::kBitShiftClusId);
  key |= clusid;
  return key;
}

TrkrDefs::cluskey
InttDefs::genClusKey(const TrkrDefs::hitsetkey hskey, const uint32_t clusid)
{
  TrkrDefs::cluskey tmp = hskey;
  TrkrDefs::cluskey key = (tmp << TrkrDefs::kBitShiftClusId);
  key |= clusid;
  return key;
}
