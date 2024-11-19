#ifndef MVTX_MVTXPIXELMASKV1_H
#define MVTX_MVTXPIXELMASKV1_H

#include "MvtxRawPixelDefs.h"
#include "MvtxPixelMask.h"

#include <iostream>

class MvtxRawHit;

class MvtxPixelMaskv1 : public MvtxPixelMask
{
    public:

        MvtxPixelMaskv1() {}
        ~MvtxPixelMaskv1() override {}

        void Reset() override {clear(); }
        void identify(std::ostream& os = std::cout) const override;
        int isValid() const override { return 1; }


        void load_from_CDB(std::string calibfile = "MVTX_HotPixelMap") override;

        void add_pixel(MvtxRawPixelDefs::mvtx_pixelkey key) override;
        void remove_pixel(MvtxRawPixelDefs::mvtx_pixelkey key) override;
        
        void clear() override { m_pixel_mask.clear(); }

        bool is_masked(MvtxRawHit* hit) const override;

        MvtxRawPixelDefs::pixelkeyvec_t get_pixel_mask() const override { return m_pixel_mask; }

    private:

        MvtxRawPixelDefs::pixelkeyvec_t m_pixel_mask = MvtxRawPixelDefs::VOID_PIXELVEC;

        ClassDefOverride(MvtxPixelMaskv1, 1);
};

#endif // MVTX_MVTXPIXELMASKV1_H
