#include <bitset>
#include "TrkrDefUtil.h"

void TrkrDefUtil::printBits(const TrkrDefs::hitsetkey key, std::ostream& os)
{
  os << "key: " << std::bitset<32>(key) << std::endl;
}

void TrkrDefUtil::printBits(const TrkrDefs::cluskey key, std::ostream& os)
{
  os << "key: " << std::bitset<64>(key) << std::endl;
}

uint8_t
TrkrDefUtil::getTrkrId(const TrkrDefs::hitsetkey key)
{
  TrkrDefs::hitsetkey tmp = (key >> kBitShiftTrkrId);
  return tmp;
}

uint8_t
TrkrDefUtil::getTrkrId(const TrkrDefs::cluskey key)
{
  TrkrDefs::hitsetkey tmp = (key >> kBitShiftClusId);
  return getTrkrId(tmp);
}

uint8_t
TrkrDefUtil::getLayer(const TrkrDefs::hitsetkey key)
{
  TrkrDefs::hitsetkey tmp = (key >> kBitShiftLayer);
  return tmp;
}

uint8_t
TrkrDefUtil::getLayer(const TrkrDefs::cluskey key)
{
  TrkrDefs::hitsetkey tmp = (key >> kBitShiftClusId);
  return getLayer(tmp);
}

uint32_t
TrkrDefUtil::getClusIndex(const TrkrDefs::cluskey key)
{
  return key;
}

uint32_t
TrkrDefUtil::getHitSetKeyFromClusKey(const TrkrDefs::cluskey key)
{
  return (key >> kBitShiftClusId);
}

TrkrDefs::hitsetkey
TrkrDefUtil::getHitSetKeyLo(const TrkrDefs::TrkrId trkrId)
{
  return genHitSetKey(trkrId, 0);
}

TrkrDefs::hitsetkey
TrkrDefUtil::getHitSetKeyHi(const TrkrDefs::TrkrId trkrId)
{
  return genHitSetKey(static_cast<TrkrDefs::TrkrId>(trkrId + 1), 0) - 1;
}

TrkrDefs::hitsetkey
TrkrDefUtil::getHitSetKeyLo(const TrkrDefs::TrkrId trkrId, const char lyr)
{
  return genHitSetKey(trkrId, lyr);
}

TrkrDefs::hitsetkey
TrkrDefUtil::getHitSetKeyHi(const TrkrDefs::TrkrId trkrId, const char lyr)
{
  return genHitSetKey(trkrId, lyr + 1) - 1;
}

TrkrDefs::cluskey
TrkrDefUtil::getClusKeyLo(const TrkrDefs::TrkrId trkrId)
{
  TrkrDefs::cluskey tmp = genHitSetKey(trkrId, 0);
  return (tmp << kBitShiftClusId);
}

TrkrDefs::cluskey
TrkrDefUtil::getClusKeyHi(const TrkrDefs::TrkrId trkrId)
{
  TrkrDefs::cluskey tmp = genHitSetKey(static_cast<TrkrDefs::TrkrId>(trkrId + 1), 0);
  return (tmp << kBitShiftClusId) - 1;
}

TrkrDefs::cluskey
TrkrDefUtil::getClusKeyLo(const TrkrDefs::TrkrId trkrId, const char lyr)
{
  TrkrDefs::cluskey tmp = genHitSetKey(trkrId, lyr);
  return (tmp << kBitShiftClusId);
}

TrkrDefs::cluskey
TrkrDefUtil::getClusKeyHi(const TrkrDefs::TrkrId trkrId, const char lyr)
{
  TrkrDefs::cluskey tmp = genHitSetKey(trkrId, lyr + 1);
  return (tmp << kBitShiftClusId) - 1;
}

TrkrDefs::hitsetkey
TrkrDefUtil::genHitSetKey(const TrkrDefs::TrkrId trkrId, const char lyr)
{
  TrkrDefs::hitsetkey tmp = trkrId;
  TrkrDefs::hitsetkey key = tmp << kBitShiftTrkrId;  // detector id
  tmp = lyr;
  key |= (tmp << kBitShiftLayer);  // layer
  return key;
}
