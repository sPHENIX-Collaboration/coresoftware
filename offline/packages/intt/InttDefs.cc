/**
 * @file intt/InttDefs.cc
 * @author D. McGlinchey
 * @date June 2018
 * @brief Implementation file for InttDefs.h
 */
#include "InttDefs.h"

uint8_t
InttDefs::getLadderId(TrkrDefs::hitsetkey key)
{
  TrkrDefs::hitsetkey tmp = (key >> InttDefs::kBitShiftLadderId);
  return tmp;
}

uint8_t
InttDefs::getLadderId(TrkrDefs::cluskey key)
{
  TrkrDefs::hitsetkey tmp = (key >> TrkrDefs::kBitShiftClusId);
  return getLadderId(tmp);
}

uint8_t
InttDefs::getSensorId(TrkrDefs::hitsetkey key)
{
  TrkrDefs::hitsetkey tmp = (key >> InttDefs::kBitShiftSensorId);
  return tmp;
}

uint8_t
InttDefs::getSensorId(TrkrDefs::cluskey key)
{
  TrkrDefs::hitsetkey tmp = (key >> TrkrDefs::kBitShiftClusId);
  return getSensorId(tmp);
}

TrkrDefs::hitkey 
InttDefs::genHitKey(const uint32_t strip)
{
  TrkrDefs::hitkey key = strip;
  return key;
}

TrkrDefs::hitsetkey
InttDefs::genHitSetKey(const uint8_t lyr, const uint8_t ladder, const uint8_t sensor)
{
  TrkrDefs::hitsetkey key = TrkrDefs::genHitSetKey(TrkrDefs::TrkrId::mvtxId, lyr);
  TrkrDefs::hitsetkey tmp = ladder;
  key |= (tmp << InttDefs::kBitShiftLadderId);
  tmp = sensor;
  key |= (tmp << InttDefs::kBitShiftSensorId);
  return key;
}

TrkrDefs::cluskey
InttDefs::genClusKey(const uint8_t lyr, const uint8_t ladder, const uint8_t sensor, const uint32_t clusid)
{
  TrkrDefs::cluskey tmp = genHitSetKey(lyr, ladder, sensor);
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
