/**
 * @file trackbase/TpcDefs.cc
 * @author D. McGlinchey
 * @date June 2018
 * @brief Implementation file for TpcDefs.h
 */
#include "TpcDefs.h"

#include "TrkrDefs.h"  // for hitsetkey, cluskey, hitkey, kBitShif...

uint8_t
TpcDefs::getSectorId(TrkrDefs::hitsetkey key)
{
  TrkrDefs::hitsetkey tmp = (key >> TpcDefs::kBitShiftSectorId);
  return tmp;
}

uint8_t
TpcDefs::getSectorId(TrkrDefs::cluskey key)
{
  TrkrDefs::hitsetkey tmp = (key >> TrkrDefs::kBitShiftClusId);
  return getSectorId(tmp);
}

uint8_t
TpcDefs::getSide(TrkrDefs::hitsetkey key)
{
  TrkrDefs::hitsetkey tmp = (key >> TpcDefs::kBitShiftSide);
  return tmp;
}

uint8_t
TpcDefs::getSide(TrkrDefs::cluskey key)
{
  TrkrDefs::hitsetkey tmp = (key >> TrkrDefs::kBitShiftClusId);
  return getSide(tmp);
}

uint16_t
TpcDefs::getPad(TrkrDefs::hitkey key)
{
  TrkrDefs::hitkey tmp = (key >> TpcDefs::kBitShiftPad);
  return tmp;
}

uint16_t
TpcDefs::getTBin(TrkrDefs::hitkey key)
{
  TrkrDefs::hitkey tmp = (key >> TpcDefs::kBitShiftTBin);
  return tmp;
}

TrkrDefs::hitkey
TpcDefs::genHitKey(const uint16_t pad, const uint16_t tbin)
{
  TrkrDefs::hitkey key = (pad << TpcDefs::kBitShiftPad);
  TrkrDefs::hitkey tmp = (tbin << TpcDefs::kBitShiftTBin);
  key |= tmp;
  return key;
}

TrkrDefs::hitsetkey
TpcDefs::genHitSetKey(const uint8_t lyr, const uint8_t sector, const uint8_t side)
{
  TrkrDefs::hitsetkey key = TrkrDefs::genHitSetKey(TrkrDefs::TrkrId::tpcId, lyr);
  TrkrDefs::hitsetkey tmp = sector;
  key |= (tmp << TpcDefs::kBitShiftSectorId);
  tmp = side;
  key |= (tmp << TpcDefs::kBitShiftSide);
  return key;
}

TrkrDefs::cluskey
TpcDefs::genClusKey(const uint8_t lyr, const uint8_t sector, const uint8_t side, const uint32_t clusid)
{
  const auto key = genHitSetKey(lyr, sector, side);
  return TrkrDefs::genClusKey( key, clusid );
}
