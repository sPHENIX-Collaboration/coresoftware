/**
 * @file trackbase/TpcDefs.cc
 * @author D. McGlinchey
 * @date June 2018
 * @brief Implementation file for TpcDefs.h
 */
#include "TpcDefs.h"

#include "TrkrDefs.h"  // for hitsetkey, cluskey, hitkey, kBitShif...

namespace
{

  // hitsetkey layout:
  //  Tpc specific lower 16 bits
  //   24 - 32  tracker id
  //   16 - 24  layer
  //   8  - 16  sector id
  //   0  -  8  side
  static constexpr unsigned int kBitShiftSectorId = 8;
  static constexpr unsigned int kBitShiftSide = 0;

  // bit shift for hitkey
  //  16 - 32 pad id
  //  0  - 16 time bin
  static constexpr unsigned int kBitShiftPad = 16;
  static constexpr unsigned int kBitShiftTBin = 0;
}

uint8_t
TpcDefs::getSectorId(TrkrDefs::hitsetkey key)
{
  TrkrDefs::hitsetkey tmp = (key >> kBitShiftSectorId);
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
  TrkrDefs::hitsetkey tmp = (key >> kBitShiftSide);
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
  TrkrDefs::hitkey tmp = (key >> kBitShiftPad);
  return tmp;
}

uint16_t
TpcDefs::getTBin(TrkrDefs::hitkey key)
{
  TrkrDefs::hitkey tmp = (key >> kBitShiftTBin);
  return tmp;
}

TrkrDefs::hitkey
TpcDefs::genHitKey(const uint16_t pad, const uint16_t tbin)
{
  TrkrDefs::hitkey key = (pad << kBitShiftPad);
  TrkrDefs::hitkey tmp = (tbin << kBitShiftTBin);
  key |= tmp;
  return key;
}

TrkrDefs::hitsetkey
TpcDefs::genHitSetKey(const uint8_t lyr, const uint8_t sector, const uint8_t side)
{
  TrkrDefs::hitsetkey key = TrkrDefs::genHitSetKey(TrkrDefs::TrkrId::tpcId, lyr);
  TrkrDefs::hitsetkey tmp = sector;
  key |= (tmp << kBitShiftSectorId);
  tmp = side;
  key |= (tmp << kBitShiftSide);
  return key;
}

TrkrDefs::cluskey
TpcDefs::genClusKey(const uint8_t lyr, const uint8_t sector, const uint8_t side, const uint32_t clusid)
{
  const auto key = genHitSetKey(lyr, sector, side);
  return TrkrDefs::genClusKey(key, clusid);
}
