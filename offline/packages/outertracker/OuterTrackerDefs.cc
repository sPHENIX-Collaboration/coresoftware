/**
 * @file outertracker/OuterTrackerDefs.cc
 * @author A Frawley
 * @date August 2019
 * @brief Implementation file for OuterTrackerDefs.h
 */
#include "OuterTrackerDefs.h"

#include <trackbase/TrkrDefs.h>

uint16_t
OuterTrackerDefs::getCol(TrkrDefs::hitkey key)
{
  TrkrDefs::hitkey tmp = (key >> OuterTrackerDefs::kBitShiftCol);
  return tmp;
}

uint16_t
OuterTrackerDefs::getRow(TrkrDefs::hitkey key)
{
  TrkrDefs::hitkey tmp = (key >> OuterTrackerDefs::kBitShiftRow);
  return tmp;
}

TrkrDefs::hitkey
OuterTrackerDefs::genHitKey(const uint16_t col, const uint16_t row)
{
  TrkrDefs::hitkey key = (col << OuterTrackerDefs::kBitShiftCol);
  TrkrDefs::hitkey tmp = (row << OuterTrackerDefs::kBitShiftRow);
  key |= tmp;
  return key;
}

TrkrDefs::hitsetkey
OuterTrackerDefs::genHitSetKey(const uint8_t lyr)
{
  TrkrDefs::hitsetkey key = TrkrDefs::genHitSetKey(TrkrDefs::TrkrId::outertrackerId, lyr);
  return key;
}


TrkrDefs::cluskey
OuterTrackerDefs::genClusKey(const uint8_t lyr, const uint32_t clusid)
{
  TrkrDefs::cluskey tmp = genHitSetKey(lyr);
  TrkrDefs::cluskey key = (tmp << TrkrDefs::kBitShiftClusId);
  key |= clusid;
  return key;
}

TrkrDefs::cluskey
OuterTrackerDefs::genClusKey(const TrkrDefs::hitsetkey hskey, const uint32_t clusid)
{
  TrkrDefs::cluskey tmp = hskey;
  TrkrDefs::cluskey key = (tmp << TrkrDefs::kBitShiftClusId);
  key |= clusid;
  return key;
}
