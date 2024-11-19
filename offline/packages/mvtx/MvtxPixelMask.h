#ifndef MVTX_MVTXPIXELMASK_H
#define MVTX_MVTXPIXELMASK_H

#include "MvtxRawPixelDefs.h"

#include <phool/PHObject.h>

#include <iostream>

class MvtxRawHit;

class MvtxPixelMask : public PHObject
{
    public:

        MvtxPixelMask() {}
        ~MvtxPixelMask() override {}

        virtual void identify(std::ostream& os) const override { os << "MvtxPixelMask" << std::endl; }
        int isValid() const override { return 0; }
        PHObject* CloneMe() const override { return nullptr; }

        virtual void load_from_CDB(std::string /*calibfile = "MVTX_HotPixelMap"*/) { return; }

        virtual void add_pixel(MvtxRawPixelDefs::mvtx_pixelkey /*key*/) { return; }
        virtual void remove_pixel(MvtxRawPixelDefs::mvtx_pixelkey /*key*/) { return; }
        
        virtual void clear() { return; }

        virtual bool is_masked(MvtxRawHit* /*hit*/) const { return false; }

        virtual MvtxRawPixelDefs::pixelkeyvec_t get_pixel_mask() const { return MvtxRawPixelDefs::VOID_PIXELVEC; }

    private:

        ClassDefOverride(MvtxPixelMask, 1);
};

#endif // MVTX_MVTXPIXELMASK_H
