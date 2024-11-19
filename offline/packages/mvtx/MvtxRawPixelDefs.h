#ifndef MVTX_MVTXRAWPIXELDEFS_H
#define MVTX_MVTXRAWPIXELDEFS_H

#include <cstdint>
#include <climits>
#include <vector>

class MvtxRawHit;

namespace MvtxRawPixelDefs 
{

    typedef uint64_t mvtx_pixelkey; // (hitsetkey << 32) | hitkey
    static MvtxRawPixelDefs::mvtx_pixelkey MVTX_VOID_PIXEL __attribute__((unused)) = UINT64_MAX;
    static const unsigned int kBitShiftLadder __attribute__((unused)) = 32;

    static const unsigned int kBitShiftLayer __attribute__((unused)) = 16;   // bitshift_trackerid - 8

    static const unsigned int kBitDummyShift __attribute__((unused)) = 24;  // 32 - 8
    static const unsigned int kMvtxIdDummy __attribute__((unused)) = 0;

    static const unsigned int kStrobeOffset __attribute__((unused)) = 16;
    static const unsigned int kBitShiftStrobeIdOffset __attribute__((unused)) = 0;
    static const uint32_t kStrobeDummy __attribute__((unused)) = kStrobeOffset << kBitShiftStrobeIdOffset;


    static const unsigned int kBitShiftStaveIdOffset __attribute__((unused)) = 9;
    static const unsigned int kBitShiftStaveIdWidth __attribute__((unused)) = 7;

    static const unsigned int kBitShiftChipIdOffset __attribute__((unused)) = 5;
    static const unsigned int kBitShiftChipIdWidth __attribute__((unused)) = 4;
        

    // bit shift for hitkey
    static const unsigned int kBitShiftCol __attribute__((unused)) = 16;
    static const unsigned int kBitShiftRow __attribute__((unused)) = 0;


    // encoding
    MvtxRawPixelDefs::mvtx_pixelkey gen_pixelkey(const uint32_t hitset_key, const uint32_t hitkey);
    MvtxRawPixelDefs::mvtx_pixelkey gen_pixelkey(MvtxRawHit *hit);
    MvtxRawPixelDefs::mvtx_pixelkey gen_pixelkey(const uint8_t layer, const uint8_t stave, const uint8_t chip, const uint16_t row, const uint16_t col);

    // decodingS
    uint32_t get_hitsetkey(const MvtxRawPixelDefs::mvtx_pixelkey key);
    uint32_t get_hitkey(const MvtxRawPixelDefs::mvtx_pixelkey key);

    uint8_t get_layer(const MvtxRawPixelDefs::mvtx_pixelkey key);
    uint8_t get_stave(const MvtxRawPixelDefs::mvtx_pixelkey key);
    uint8_t get_chip(const MvtxRawPixelDefs::mvtx_pixelkey key);
    uint16_t get_row(const MvtxRawPixelDefs::mvtx_pixelkey key);
    uint16_t get_col(const MvtxRawPixelDefs::mvtx_pixelkey key);

    typedef std::vector<MvtxRawPixelDefs::mvtx_pixelkey> pixelkeyvec_t;
    static const MvtxRawPixelDefs::pixelkeyvec_t VOID_PIXELVEC __attribute__((unused)) = {};

    typedef std::pair<MvtxRawPixelDefs::mvtx_pixelkey, uint64_t> pixelhitpair_t;
    typedef std::vector<MvtxRawPixelDefs::pixelhitpair_t> pixelhitpair_vec_t;

    static const unsigned int kNpixelsTotal __attribute__((unused)) = 226492416;
}

#endif // MVTX_MVTXRAWPIXELDEFS_H
