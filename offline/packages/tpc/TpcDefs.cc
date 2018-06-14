/**
 * @file tpc/TpcDefs.cc
 * @author D. McGlinchey
 * @date June 2018
 * @brief Implementation file for TpcDefs.h
 */
#include "TpcDefs.h"

uint8_t
TpcDefs::getStaveId(TrkrDefs::hitsetkey key)
{
  TrkrDefs::hitsetkey tmp = (key >> TpcDefs::kBitShiftStaveId);
  return tmp;
}

uint8_t
TpcDefs::getStaveId(TrkrDefs::cluskey key)
{
  TrkrDefs::hitsetkey tmp = (key >> TrkrDefs::kBitShiftClusId);
  return getStaveId(tmp);
}

uint8_t
TpcDefs::getChipId(TrkrDefs::hitsetkey key)
{
  TrkrDefs::hitsetkey tmp = (key >> TpcDefs::kBitShiftChipId);
  return tmp;
}

uint8_t
TpcDefs::getChipId(TrkrDefs::cluskey key)
{
  TrkrDefs::hitsetkey tmp = (key >> TrkrDefs::kBitShiftClusId);
  return getChipId(tmp);
}

uint16_t
TpcDefs::getCol(TrkrDefs::hitkey key)
{
  TrkrDefs::hitkey tmp = (key >> TpcDefs::kBitShiftCol);
  return tmp;
}

uint16_t
TpcDefs::getRow(TrkrDefs::hitkey key)
{
  TrkrDefs::hitkey tmp = (key >> TpcDefs::kBitShiftRow);
  return tmp;
}

TrkrDefs::hitkey 
TpcDefs::genHitKey(const uint16_t col, const uint16_t row)
{
  TrkrDefs::hitkey key = (col << TpcDefs::kBitShiftCol);
  TrkrDefs::hitkey tmp = (row << TpcDefs::kBitShiftRow);
  key |= tmp;
  return key;
}

TrkrDefs::hitsetkey
TpcDefs::genHitSetKey(const char lyr, const uint8_t stave, const uint8_t chip)
{
  TrkrDefs::hitsetkey key = TrkrDefs::genHitSetKey(TrkrDefs::TrkrId::mvtxId, lyr);
  TrkrDefs::hitsetkey tmp = stave;
  key |= (tmp << TpcDefs::kBitShiftStaveId);
  tmp = chip;
  key |= (tmp << TpcDefs::kBitShiftChipId);
  return key;
}

TrkrDefs::cluskey
TpcDefs::genClusKey(const char lyr, const uint8_t stave, const uint8_t chip, const uint32_t clusid)
{
  TrkrDefs::cluskey tmp = genHitSetKey(lyr, stave, chip);
  TrkrDefs::cluskey key = (tmp << TrkrDefs::kBitShiftClusId);
  key |= clusid;
  return key;
}

TrkrDefs::cluskey
TpcDefs::genClusKey(const TrkrDefs::hitsetkey hskey, const uint32_t clusid)
{
  TrkrDefs::cluskey tmp = hskey;
  TrkrDefs::cluskey key = (tmp << TrkrDefs::kBitShiftClusId);
  key |= clusid;
  return key;
}
