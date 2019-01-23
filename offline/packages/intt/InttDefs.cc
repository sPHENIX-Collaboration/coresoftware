/**
 * @file intt/InttDefs.cc
 * @author D. McGlinchey
 * @date June 2018
 * @brief Implementation file for InttDefs.h
 */
#include "InttDefs.h"

uint8_t
InttDefs::getLadderZId(TrkrDefs::hitsetkey key)
{
  TrkrDefs::hitsetkey tmp = (key >> InttDefs::kBitShiftLadderZId);
  return tmp;
}

uint8_t
InttDefs::getLadderZId(TrkrDefs::cluskey key)
{
  TrkrDefs::hitsetkey tmp = (key >> TrkrDefs::kBitShiftClusId);
  return getLadderZId(tmp);
}

uint8_t
InttDefs::getLadderPhiId(TrkrDefs::hitsetkey key)
{
  TrkrDefs::hitsetkey tmp = (key >> InttDefs::kBitShiftLadderPhiId);
  return tmp;
}

uint8_t
InttDefs::getLadderPhiId(TrkrDefs::cluskey key)
{
  TrkrDefs::hitsetkey tmp = (key >> TrkrDefs::kBitShiftClusId);
  return getLadderPhiId(tmp);
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
InttDefs::genHitSetKey(const uint8_t lyr, const uint8_t ladder_z_index, uint8_t ladder_phi_index)
{
  TrkrDefs::hitsetkey key = TrkrDefs::genHitSetKey(TrkrDefs::TrkrId::inttId, lyr);
  TrkrDefs::hitsetkey tmp = ladder_z_index;
  key |= (tmp << InttDefs::kBitShiftLadderZId);
  tmp = ladder_phi_index;
  key |= (tmp << InttDefs::kBitShiftLadderPhiId);
  return key;
}

TrkrDefs::cluskey
InttDefs::genClusKey(const uint8_t lyr, const uint8_t ladder_z_index, const uint8_t ladder_phi_index, const uint32_t clusid)
{
  TrkrDefs::cluskey tmp = genHitSetKey(lyr, ladder_z_index, ladder_phi_index);
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
