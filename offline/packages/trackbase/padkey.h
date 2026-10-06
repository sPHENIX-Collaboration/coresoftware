#ifndef TRACKBASE_PADKEY_H
#define TRACKBASE_PADKEY_H

#include <cassert>
#include <cstdint>
#include <limits>

/**
 * @brief packed key for aggregated TPC laser pads
 *
 * bit layout (uint32_t, bits 28-31 unused):
 *
 *   31 ... 27 | 26   | 25 .. 22 | 21 .. 16 | 15 .. 0
 *   unused    | side | sector   | rlayer   | phibin
 *                                 ^^^^^^
 *                    top 2 bits of rlayer == module (R1/R2/R3)
 *
 * rlayer is the layer relative to the first TPC layer (0-47), so the
 * module needs no field of its own: rlayer>>4 is the module, and the
 * key remains ordered side -> sector -> module -> layer -> phi.
 *
 * Consequence: each of the 72 side/sector/module blocks is a contiguous
 * range of keys, so a block can be handed to a thread as a pair of map
 * iterators via [block_begin(blk), block_end(blk)).
 */
namespace padkey
{
  using key = uint32_t;

  //! first TPC layer in the global layer numbering
  static constexpr uint16_t kFirstTpcLayer   = 7;
  static constexpr uint16_t kNTpcLayers      = 48;
  static constexpr uint16_t kLayersPerModule = 16;
  static constexpr uint8_t  kNModules        = 3;
  static constexpr uint8_t  kNSectors        = 12;
  static constexpr uint8_t  kNSides          = 2;
  static constexpr uint8_t  kNBlocks         = kNSides * kNSectors * kNModules;  // 72

  // the module-from-rlayer shortcut below is only valid for a power-of-two
  // number of layers per module
  static_assert(kLayersPerModule == 16, "module extraction assumes 16 layers/module");
  static_assert(kNTpcLayers == kNModules * kLayersPerModule, "layer/module mismatch");

  static constexpr unsigned kPhiBits    = 16;
  static constexpr unsigned kLayerBits  = 6;   // relative layer, 0-47
  static constexpr unsigned kSectorBits = 4;   // 0-11
  static constexpr unsigned kSideBits   = 1;

  static constexpr unsigned kPhiShift    = 0;
  static constexpr unsigned kLayerShift  = kPhiShift + kPhiBits;        // 16
  static constexpr unsigned kBlockShift  = 20;                          // module = top 2 bits of rlayer
  static constexpr unsigned kSectorShift = kLayerShift + kLayerBits;    // 22
  static constexpr unsigned kSideShift   = kSectorShift + kSectorBits;  // 26

  static constexpr unsigned kModuleShiftInLayer = kBlockShift - kLayerShift;  // 4

  static constexpr key kPhiMask       = (1U << kPhiBits) - 1;
  static constexpr key kLayerMask     = (1U << kLayerBits) - 1;
  static constexpr key kSectorMask    = (1U << kSectorBits) - 1;
  static constexpr key kSideMask      = (1U << kSideBits) - 1;
  static constexpr key kLayerFieldMask = kLayerMask << kLayerShift;
  static constexpr key kBlockMask     = (1U << (kSideBits + kSectorBits + 2)) - 1;  // 7 bits

  //! step between adjacent layers within a module
  static constexpr key kLayerStep = 1U << kLayerShift;

  //! sentinel — the unused high bits guarantee this is not a legal key
  static constexpr key kInvalidKey = std::numeric_limits<key>::max();

  // ---------------------------------------------------------------- encode

  //! block index 0-71, ordered side -> sector -> module
  inline uint8_t genblock(uint8_t side, uint8_t sector, uint8_t module)
  {
    return static_cast<uint8_t>(((side & kSideMask) << 6U)
                              | ((sector & kSectorMask) << 2U)
                              | (module & 0x3U));
  }

  //! key from absolute layer (7-54), side (0-1), sector (0-11), phi bin
  inline key genkey(uint16_t layer, uint8_t side, uint8_t sector, uint16_t phibin)
  {
    assert(layer >= kFirstTpcLayer && layer < kFirstTpcLayer + kNTpcLayers);
    assert(sector < kNSectors);
    const key rlayer = static_cast<key>(layer - kFirstTpcLayer);
    return (static_cast<key>(side & kSideMask)     << kSideShift)
         | (static_cast<key>(sector & kSectorMask) << kSectorShift)
         | (rlayer                                 << kLayerShift)
         | (static_cast<key>(phibin) & kPhiMask);
  }

  // ---------------------------------------------------------------- decode

  inline uint8_t  get_side(key k)   { return (k >> kSideShift) & kSideMask; }
  inline uint8_t  get_sector(key k) { return (k >> kSectorShift) & kSectorMask; }
  inline uint16_t get_rlayer(key k) { return (k >> kLayerShift) & kLayerMask; }
  inline uint8_t  get_module(key k) { return get_rlayer(k) >> kModuleShiftInLayer; }
  inline uint16_t get_layer(key k)  { return kFirstTpcLayer + get_rlayer(k); }
  inline uint16_t get_phibin(key k) { return k & kPhiMask; }

  //! block index 0-71 this key belongs to
  inline uint8_t get_block(key k) { return (k >> kBlockShift) & kBlockMask; }

  // ---------------------------------------------------- block partitioning

  inline key block_begin(uint8_t blk) { return static_cast<key>(blk) << kBlockShift; }
  inline key block_end(uint8_t blk)   { return static_cast<key>(blk + 1) << kBlockShift; }

  inline uint8_t block_side(uint8_t blk)   { return (blk >> 6U) & kSideMask; }
  inline uint8_t block_sector(uint8_t blk) { return (blk >> 2U) & kSectorMask; }
  inline uint8_t block_module(uint8_t blk) { return blk & 0x3U; }

  //! blocks are 0-71 but the 7-bit field allows 0-127; module 3 never occurs
  inline bool is_valid_block(uint8_t blk) { return block_module(blk) < kNModules; }

  // ------------------------------------------------------------- neighbors
  // Both stop at physical boundaries: sector edges in phi, module gaps in
  // layer. Neither ever leaves the block, so they are safe to call from a
  // per-block worker thread.

  //! phi neighbor within the same sector; false at a sector edge
  inline bool phi_neighbor(key k, int d, uint16_t nphi_per_sector, key& out)
  {
    const int p = static_cast<int>(get_phibin(k)) + d;
    if (p < 0 || p >= static_cast<int>(nphi_per_sector)) return false;
    out = (k & ~kPhiMask) | static_cast<key>(p);
    return true;
  }

  //! layer neighbor within the same module; false at a module boundary
  inline bool layer_neighbor(key k, int d, key& out)
  {
    const int L = static_cast<int>(get_rlayer(k)) + d;
    if (L < 0 || L >= static_cast<int>(kNTpcLayers)) return false;
    if ((static_cast<unsigned>(L) >> kModuleShiftInLayer) != get_module(k)) return false;
    out = (k & ~kLayerFieldMask) | (static_cast<key>(L) << kLayerShift);
    return true;
  }
}

#endif  // TRACKBASE_PADKEY_H