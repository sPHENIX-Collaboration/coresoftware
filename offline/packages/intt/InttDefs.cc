/**
 * @file intt/InttDefs.cc
 * @author D. McGlinchey
 * @date June 2018
 * @brief Implementation file for InttDefs.h
 */
#include "InttDefs.h"

#include <trackbase/TrkrDefs.h>

uint8_t
InttDefs::getLadderZId(TrkrDefs::hitsetkey key)
{
  TrkrDefs::hitsetkey tmp = (key >> InttDefs::kBitShiftLadderZIdOffset);
  // clear the bits not associated with the ladderZId
  uint8_t tmp1 = tmp;
  tmp1 = (tmp1 <<  (InttDefs::kBitShiftLadderZIdWidth - InttDefs::kBitShiftLadderZIdOffset));
  tmp1 = (tmp1 >>  (InttDefs::kBitShiftLadderZIdWidth - InttDefs::kBitShiftLadderZIdOffset));
  return tmp1;
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
  TrkrDefs::hitsetkey tmp = (key >> InttDefs::kBitShiftLadderPhiIdOffset);
  // clear the bits not associated with the ladderPhiId
  uint8_t tmp1 = tmp;
  tmp1 = (tmp1 <<  (InttDefs::kBitShiftLadderPhiIdWidth - InttDefs::kBitShiftLadderPhiIdOffset));
  tmp1 = (tmp1 >>  (InttDefs::kBitShiftLadderPhiIdWidth - InttDefs::kBitShiftLadderPhiIdOffset));
  return tmp1;
}

uint16_t
InttDefs::getTimeBucketId(TrkrDefs::hitsetkey key)
{
  TrkrDefs::hitsetkey tmp = (key >> InttDefs::kBitShiftTimeBucketIdOffset);
  // clear the bits not associated with the TimeBucketId
  uint16_t tmp1 = tmp;
  tmp1 = (tmp1 <<  (InttDefs::kBitShiftTimeBucketIdWidth - InttDefs::kBitShiftTimeBucketIdOffset));
  tmp1 = (tmp1 >>  (InttDefs::kBitShiftTimeBucketIdWidth - InttDefs::kBitShiftTimeBucketIdOffset));
  return tmp1;
}

uint8_t
InttDefs::getLadderPhiId(TrkrDefs::cluskey key)
{
  TrkrDefs::hitsetkey tmp = (key >> TrkrDefs::kBitShiftClusId);
  return getLadderPhiId(tmp);
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
InttDefs::genHitSetKey(const uint8_t lyr, const uint8_t ladder_z_index, uint8_t ladder_phi_index, uint16_t crossing)
{
  TrkrDefs::hitsetkey key = TrkrDefs::genHitSetKey(TrkrDefs::TrkrId::inttId, lyr);
  TrkrDefs::hitsetkey tmp = ladder_z_index;
  key |= (tmp << InttDefs::kBitShiftLadderZIdOffset);
  tmp = ladder_phi_index;
  key |= (tmp << InttDefs::kBitShiftLadderPhiIdOffset);
  tmp = crossing;
  key |= (tmp << InttDefs::kBitShiftTimeBucketIdOffset);
  return key;
}

TrkrDefs::cluskey
InttDefs::genClusKey(const uint8_t lyr, const uint8_t ladder_z_index, const uint8_t ladder_phi_index, const uint16_t crossing, const uint32_t clusid)
{
  TrkrDefs::cluskey tmp = genHitSetKey(lyr, ladder_z_index, ladder_phi_index, crossing);
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
