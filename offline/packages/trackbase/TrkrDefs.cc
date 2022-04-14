#include "TrkrDefs.h"

#include <bitset>


void TrkrDefs::printBits(const TrkrDefs::hitsetkey key, std::ostream& os)
{
  os << "key: " << std::bitset<32>(key) << std::endl;
}

void TrkrDefs::printBits(const TrkrDefs::cluskey key, std::ostream& os)
{
  os << "key: " << std::bitset<64>(key) << std::endl;
}

uint8_t
TrkrDefs::getTrkrId(const TrkrDefs::hitsetkey key)
{
  TrkrDefs::hitsetkey tmp = (key >> kBitShiftTrkrId);
  return tmp;
}

uint8_t
TrkrDefs::getTrkrId(const TrkrDefs::cluskey key)
{
  TrkrDefs::hitsetkey tmp = (key >> kBitShiftClusId);
  return getTrkrId(tmp);
}

uint8_t
TrkrDefs::getLayer(const TrkrDefs::hitsetkey key)
{
  TrkrDefs::hitsetkey tmp = (key >> kBitShiftLayer);
  return tmp;
}

uint8_t
TrkrDefs::getLayer(const TrkrDefs::cluskey key)
{
  TrkrDefs::hitsetkey tmp = (key >> kBitShiftClusId);
  return getLayer(tmp);
}

uint32_t
TrkrDefs::getClusIndex(const TrkrDefs::cluskey key)
{
  return key;
}

uint32_t
TrkrDefs::getHitSetKeyFromClusKey(const TrkrDefs::cluskey key)
{
  return (key >> kBitShiftClusId);
}

TrkrDefs::hitsetkey
TrkrDefs::getHitSetKeyLo(const TrkrDefs::TrkrId trkrId)
{
  return genHitSetKey(trkrId, 0);
}

TrkrDefs::hitsetkey
TrkrDefs::getHitSetKeyHi(const TrkrDefs::TrkrId trkrId)
{
  return genHitSetKey(static_cast<TrkrDefs::TrkrId>(trkrId + 1), 0) - 1;
}

TrkrDefs::hitsetkey
TrkrDefs::getHitSetKeyLo(const TrkrDefs::TrkrId trkrId, const uint8_t lyr)
{
  return genHitSetKey(trkrId, lyr);
}

TrkrDefs::hitsetkey
TrkrDefs::getHitSetKeyHi(const TrkrDefs::TrkrId trkrId, const uint8_t lyr)
{
  return genHitSetKey(trkrId, lyr + 1) - 1;
}

TrkrDefs::cluskey
TrkrDefs::getClusKeyLo(const TrkrDefs::TrkrId trkrId)
{
  const TrkrDefs::cluskey tmp = genHitSetKey(trkrId, 0);
  return (tmp << kBitShiftClusId);
}

TrkrDefs::cluskey
TrkrDefs::getClusKeyHi(const TrkrDefs::TrkrId trkrId)
{
  const TrkrDefs::cluskey tmp = genHitSetKey(static_cast<TrkrDefs::TrkrId>(trkrId + 1), 0);
  return (tmp << kBitShiftClusId) - 1;
}

TrkrDefs::cluskey
TrkrDefs::getClusKeyLo(const TrkrDefs::TrkrId trkrId, const uint8_t lyr)
{
  const TrkrDefs::cluskey tmp = genHitSetKey(trkrId, lyr);
  return (tmp << kBitShiftClusId);
}

TrkrDefs::cluskey
TrkrDefs::getClusKeyHi(const TrkrDefs::TrkrId trkrId, const uint8_t lyr)
{
  const TrkrDefs::cluskey tmp = genHitSetKey(trkrId, lyr + 1);
  return (tmp << kBitShiftClusId) - 1;
}

TrkrDefs::hitsetkey
TrkrDefs::genHitSetKey(const TrkrDefs::TrkrId trkrId, const uint8_t lyr)
{
  TrkrDefs::hitsetkey tmp = trkrId;
  TrkrDefs::hitsetkey key = tmp << kBitShiftTrkrId;  // detector id
  tmp = lyr;
  key |= (tmp << kBitShiftLayer);  // layer
  return key;
}

TrkrDefs::cluskey
TrkrDefs::genClusKey(const TrkrDefs::hitsetkey hskey, const uint32_t clusid)
{
  const TrkrDefs::cluskey tmp = hskey;
  TrkrDefs::cluskey key = (tmp << TrkrDefs::kBitShiftClusId);
  key |= clusid;
  return key;
}

uint8_t TrkrDefs::getPhiElement(TrkrDefs::hitsetkey key)
{
  const TrkrDefs::hitsetkey tmp = (key >> TrkrDefs::kBitShiftPhiElement);
  return tmp;
}



uint8_t TrkrDefs::getZElement(TrkrDefs::hitsetkey key)
{
  const TrkrDefs::hitsetkey tmp = (key >> TrkrDefs::kBitShiftZElement);
  return tmp;
}

uint8_t TrkrDefs::getPhiElement(TrkrDefs::cluskey key)
{
  const TrkrDefs::hitsetkey tmp = (key >> TrkrDefs::kBitShiftClusId);
  return getPhiElement(tmp);
}

uint8_t TrkrDefs::getZElement(TrkrDefs::cluskey key)//side
{
  const TrkrDefs::hitsetkey tmp = (key >> TrkrDefs::kBitShiftClusId);
  return getZElement(tmp);
}
