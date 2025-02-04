#include "MvtxPixelDefs.h"

#include <trackbase/MvtxDefs.h>
#include <trackbase/TrkrDefs.h>

#include <cstdint>

// MvtPixelDefs namespace
//==============================================================================
MvtxPixelDefs::pixelkey MvtxPixelDefs::gen_pixelkey(const uint32_t hitset_key, const uint32_t hitkey)
{
  return ((uint64_t) hitset_key << MvtxPixelDefs::kBitShiftLadder) | ((uint64_t) hitkey & 0xFFFFFFFF);
}

uint32_t MvtxPixelDefs::get_hitsetkey(const MvtxPixelDefs::pixelkey pkey)
{
  return (uint32_t) (pkey >> MvtxPixelDefs::kBitShiftLadder);
}

uint32_t MvtxPixelDefs::get_hitkey(const MvtxPixelDefs::pixelkey pkey)
{
  return (uint32_t) pkey & 0xFFFFFFFF;
}

MvtxPixelDefs::pixelkey MvtxPixelDefs::gen_pixelkey_from_coors(const uint8_t layer, const uint8_t stave, const uint8_t chip, const uint16_t row, const uint16_t col)
{
  return MvtxPixelDefs::gen_pixelkey(MvtxDefs::genHitSetKey(layer, stave, chip, 0), MvtxDefs::genHitKey(col, row));
}

MvtxPixelDefs::pixelkey MvtxPixelDefs::gen_pixelkey(MvtxRawHit *hit)
{
  if (!hit)
  {
    return MvtxPixelDefs::VOID_PIXEL;
  }
  const TrkrDefs::hitkey this_hitkey = MvtxDefs::genHitKey(hit->get_col(), hit->get_row());
  const TrkrDefs::hitsetkey this_hitsetkey = MvtxDefs::genHitSetKey(hit->get_layer_id(), hit->get_stave_id(), hit->get_chip_id(), 0);
  MvtxPixelDefs::pixelkey this_pixelkey = MvtxPixelDefs::gen_pixelkey(this_hitsetkey, this_hitkey);
  return this_pixelkey;
}

unsigned int MvtxPixelDefs::get_layer(const MvtxPixelDefs::pixelkey pkey)
{
  uint32_t hitsetkey = get_hitsetkey(pkey);
  return TrkrDefs::getLayer(hitsetkey);
}

unsigned int MvtxPixelDefs::get_stave(const MvtxPixelDefs::pixelkey pkey)
{
  uint32_t hitsetkey = get_hitsetkey(pkey);
  return MvtxDefs::getStaveId(hitsetkey);
}

unsigned int MvtxPixelDefs::get_chip(const MvtxPixelDefs::pixelkey pkey)
{
  uint32_t hitsetkey = get_hitsetkey(pkey);
  return MvtxDefs::getChipId(hitsetkey);
}

unsigned int MvtxPixelDefs::get_row(const MvtxPixelDefs::pixelkey pkey)
{
  uint32_t hitkey = get_hitkey(pkey);
  return MvtxDefs::getRow(hitkey);
}

unsigned int MvtxPixelDefs::get_col(const MvtxPixelDefs::pixelkey pkey)
{
  uint32_t hitkey = get_hitkey(pkey);
  return MvtxDefs::getCol(hitkey);
}