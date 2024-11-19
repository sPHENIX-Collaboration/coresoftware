#ifndef MVTX_MVTXPIXELMASKER_H
#define MVTX_MVTXPIXELMASKER_H


#include "MvtxRawPixelDefs.h"

#include <fun4all/SubsysReco.h>

#include <cstdint>
#include <iostream>
#include <string>

class PHCompositeNode;
class MvtxPixelMaskv1;

class MvtxRawPixelMasker : public SubsysReco
{
    public:

        MvtxRawPixelMasker(const std::string& name = "MvtxRawPixelMasker") : SubsysReco(name) {}
        ~MvtxRawPixelMasker() override {}

        // standard Fun4All functions
        int InitRun(PHCompositeNode* topNode) override;
        int process_event(PHCompositeNode* topNode) override;
        int End(PHCompositeNode* /*topNode*/ ) override;

        void set_threshold(double threshold) { m_threshold = threshold; }
        void set_update_every_nstrobe(unsigned int n) { m_strobe_update = n; }
        void load_starting_from_CDB(bool b, const std::string& calibfile = "MVTX_HotPixelMap" ) { 
            m_load_from_cdb = b;
            m_calibfile = calibfile; 
        }

        void enable_dynamic_masking(bool b) { m_dynamic_mode = b; }
                
    private:

        bool m_dynamic_mode{false};

        int m_max_masked_pixels{0};
        uint64_t m_last_strobe{0};
        uint64_t m_current_strobe_counter{0};
        MvtxRawPixelDefs::pixelhitpair_vec_t m_current_pixel_hits{};

        double m_threshold{1e-8};
        unsigned int m_strobe_update{100000};

        bool m_load_from_cdb{false};
        std::string m_calibfile{"MVTX_HotPixelMap"};

        std::string m_hot_pixel_node_name{"MvtxHotPixelMask"};

        unsigned int m_npixels_masked{0};
        int CreateNode(PHCompositeNode* topNode);
        void UpdateMaskedPixels();
        void addHit(MvtxRawHit* hit, const unsigned int nhits=1);

};

#endif // MVTX_MVTXPIXELMASKER_H
