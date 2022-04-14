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
  TrkrDefs::hitsetkey tmp = (key >> InttDefs::kBitShiftLadderZIdOffset);
  // clear the bits not associated with the ladderZId
  uint8_t tmp1 = tmp;
  tmp1 = (tmp1 <<  (8 - InttDefs::kBitShiftLadderZIdWidth));
  tmp1 = (tmp1 >>  (8 - InttDefs::kBitShiftLadderZIdWidth));
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
  tmp1 = (tmp1 <<  (8 - InttDefs::kBitShiftLadderPhiIdWidth));
  tmp1 = (tmp1 >>  (8 - InttDefs::kBitShiftLadderPhiIdWidth));
  return tmp1;
}

uint8_t
InttDefs::getLadderPhiId(TrkrDefs::cluskey key)
{
  TrkrDefs::hitsetkey tmp = (key >> TrkrDefs::kBitShiftClusId);
  return getLadderPhiId(tmp);
}

int
InttDefs::getTimeBucketId(TrkrDefs::hitsetkey key)
{
  TrkrDefs::hitsetkey tmp = (key >> InttDefs::kBitShiftTimeBucketIdOffset);
  // clear the bits not associated with the TimeBucketId
  uint16_t tmp1 = tmp;
  tmp1 = (tmp1 <<  (16 - InttDefs::kBitShiftTimeBucketIdWidth));
  tmp1 = (tmp1 >>  (16 - InttDefs::kBitShiftTimeBucketIdWidth));

  int tmp2 = (int) tmp1 - crossingOffset;   // get back to signed crossing

  return tmp2;
}

int
InttDefs::getTimeBucketId(TrkrDefs::cluskey key)
{
  TrkrDefs::hitsetkey tmp = (key >> TrkrDefs::kBitShiftClusId);
  return getTimeBucketId(tmp);
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
InttDefs::genHitSetKey(const uint8_t lyr, const uint8_t ladder_z_index, uint8_t ladder_phi_index, const int crossing_in)
{
  TrkrDefs::hitsetkey key = TrkrDefs::genHitSetKey(TrkrDefs::TrkrId::inttId, lyr);

  // offset crossing to make it positive, fit inside 10 bits
  int crossing = crossing_in + crossingOffset;
  if(crossing < 0) crossing = 0;  
  if(crossing > 1023) crossing = 1023;
  unsigned int ucrossing = (unsigned int) crossing;

  TrkrDefs::hitsetkey tmp = ladder_z_index;
  key |= (tmp << InttDefs::kBitShiftLadderZIdOffset);
  tmp = ladder_phi_index;
  key |= (tmp << InttDefs::kBitShiftLadderPhiIdOffset);
  tmp = ucrossing;
  key |= (tmp << InttDefs::kBitShiftTimeBucketIdOffset);

  return key;
}

TrkrDefs::cluskey
InttDefs::genClusKey(const uint8_t lyr, const uint8_t ladder_z_index, const uint8_t ladder_phi_index, const int crossing, const uint32_t clusid)
{
  TrkrDefs::cluskey key = genHitSetKey(lyr, ladder_z_index, ladder_phi_index, crossing);
  return TrkrDefs::genClusKey( key, clusid );
}

TrkrDefs::hitsetkey
InttDefs::resetCrossingHitSetKey(const TrkrDefs::hitsetkey hitsetkey)
{
  // Note: this method uses the fact that the crossing is in the first 10 bits
   TrkrDefs::hitsetkey tmp = hitsetkey;
   // zero the crossing bits by shifting them out of the word, then shift back
   tmp = (tmp >>  InttDefs::kBitShiftTimeBucketIdWidth);
   tmp = (tmp << InttDefs::kBitShiftTimeBucketIdWidth);
   unsigned int zero_crossing = crossingOffset;
   tmp |= (zero_crossing << InttDefs::kBitShiftTimeBucketIdOffset);

  return tmp;
}
