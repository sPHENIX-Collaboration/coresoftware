#include "MvtxRawPixelDefs.h"

#include <ffarawobjects/MvtxRawHit.h>


MvtxRawPixelDefs::mvtx_pixelkey MvtxRawPixelDefs::gen_pixelkey(const uint32_t hitset_key, const uint32_t hitkey)
{
    return ((uint64_t) hitset_key << MvtxRawPixelDefs::kBitShiftLadder) | ((uint64_t) hitkey & 0xFFFFFFFF);
}

MvtxRawPixelDefs::mvtx_pixelkey MvtxRawPixelDefs::gen_pixelkey(const uint8_t layer, 
                    const uint8_t stave, 
                    const uint8_t chip, 
                    const uint16_t row, 
                    const uint16_t col)
{
    const uint32_t this_hitkey = ((uint32_t) col << MvtxRawPixelDefs::kBitShiftCol) | ((uint32_t) row << MvtxRawPixelDefs::kBitShiftRow);
    const uint32_t this_hitsetkey =  ( MvtxRawPixelDefs::kMvtxIdDummy << MvtxRawPixelDefs::kBitDummyShift )  
                            | ((uint32_t) layer << MvtxRawPixelDefs::kBitShiftLayer) 
                            | ((uint32_t) stave << MvtxRawPixelDefs::kBitShiftStaveIdOffset) 
                            | ((uint32_t) chip << MvtxRawPixelDefs::kBitShiftChipIdOffset) 
                            | ((uint32_t) MvtxRawPixelDefs::kStrobeDummy << MvtxRawPixelDefs::kBitShiftStrobeIdOffset);   
    return MvtxRawPixelDefs::gen_pixelkey(this_hitsetkey, this_hitkey);
}

MvtxRawPixelDefs::mvtx_pixelkey MvtxRawPixelDefs::gen_pixelkey( MvtxRawHit * hit )
{
    if( !hit ) { return MvtxRawPixelDefs::MVTX_VOID_PIXEL; }
    return MvtxRawPixelDefs::gen_pixelkey(hit->get_layer_id(), 
                            hit->get_stave_id(), 
                            hit->get_chip_id(), 
                            hit->get_row(), 
                            hit->get_col());
}

uint32_t MvtxRawPixelDefs::get_hitsetkey(const MvtxRawPixelDefs::mvtx_pixelkey pkey)
{
    return (uint32_t) (pkey >> MvtxRawPixelDefs::kBitShiftLadder);
}

uint32_t MvtxRawPixelDefs::get_hitkey(const MvtxRawPixelDefs::mvtx_pixelkey pkey)
{
    return (uint32_t) pkey & 0xFFFFFFFF;
}

uint8_t MvtxRawPixelDefs::get_layer(const MvtxRawPixelDefs::mvtx_pixelkey pkey)
{
    uint32_t hitsetkey = MvtxRawPixelDefs::get_hitsetkey(pkey);
    uint32_t layer = (hitsetkey >> MvtxRawPixelDefs::kBitShiftLayer) & 0xFF;
    return (uint8_t) layer;
}

uint8_t MvtxRawPixelDefs::get_stave(const MvtxRawPixelDefs::mvtx_pixelkey pkey) 
{
    uint32_t hitsetkey = MvtxRawPixelDefs::get_hitsetkey(pkey);
    return (hitsetkey >> MvtxRawPixelDefs::kBitShiftStaveIdOffset) & ((1 << MvtxRawPixelDefs::kBitShiftStaveIdWidth) - 1);
}

uint8_t MvtxRawPixelDefs::get_chip(const MvtxRawPixelDefs::mvtx_pixelkey pkey) 
{
    uint32_t hitsetkey = MvtxRawPixelDefs::get_hitsetkey(pkey);
    return (hitsetkey >> MvtxRawPixelDefs::kBitShiftChipIdOffset) & ((1 << MvtxRawPixelDefs::kBitShiftChipIdWidth) - 1);
}

uint16_t MvtxRawPixelDefs::get_row(const MvtxRawPixelDefs::mvtx_pixelkey pkey)
{
    uint32_t hitkey = get_hitkey(pkey);
    return (uint16_t) (hitkey >> MvtxRawPixelDefs::kBitShiftRow);
}

uint16_t MvtxRawPixelDefs::get_col(const MvtxRawPixelDefs::mvtx_pixelkey pkey)
{
    uint32_t hitkey = get_hitkey(pkey);
    return (uint16_t) (hitkey >> MvtxRawPixelDefs::kBitShiftCol);
}
