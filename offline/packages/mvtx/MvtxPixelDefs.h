/*!
 * \file MvtxPixelDefs.h
 * \brief Defines pixel ids and keys
 * \author Tanner Mengel <tmengel@bnl.gov>
 * \version $Version: 2.0.1 $
 * \date $Date: 05/23/2025.
 */


#ifndef MVTX_MVTXPIXELDEFS_H
#define MVTX_MVTXPIXELDEFS_H

#include <climits>
#include <cstdint>
#include <memory>

#include <ffarawobjects/MvtxRawHit.h>

namespace MvtxPixelDefs
{

  typedef uint64_t pixelkey;  // (hitsetkey << 32) | hitkey
  static MvtxPixelDefs::pixelkey VOID_PIXEL __attribute__((unused)) = UINT64_MAX;
  static const unsigned int kBitShiftLadder __attribute__((unused)) = 32;

  MvtxPixelDefs::pixelkey gen_pixelkey(const uint32_t hitset_key, const uint32_t hitkey);
  uint32_t get_hitsetkey(const MvtxPixelDefs::pixelkey key);
  uint32_t get_hitkey(const MvtxPixelDefs::pixelkey key);
  MvtxPixelDefs::pixelkey gen_pixelkey_from_coors(const uint8_t layer, const uint8_t stave, const uint8_t chip, const uint16_t row, const uint16_t col);
  MvtxPixelDefs::pixelkey gen_pixelkey(MvtxRawHit *hit);

  unsigned int get_layer(const MvtxPixelDefs::pixelkey key);
  unsigned int get_stave(const MvtxPixelDefs::pixelkey key);
  unsigned int get_chip(const MvtxPixelDefs::pixelkey key);
  unsigned int get_row(const MvtxPixelDefs::pixelkey key);
  unsigned int get_col(const MvtxPixelDefs::pixelkey key);
}  // namespace MvtxPixelDefs

#endif  // MVTX_MVTXPIXELDEFS_H
