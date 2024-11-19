#include "MvtxPixelMaskv1.h"

#include <ffamodules/CDBInterface.h> 
#include <cdbobjects/CDBTTree.h>

#include <ffarawobjects/MvtxRawHit.h>
#include <ffarawobjects/MvtxRawHitv1.h>

#include <algorithm>
#include <set>

void MvtxPixelMaskv1::identify(std::ostream& os) const
{
    os << "MvtxPixelMaskv1" << std::endl;
    os << "Number of masked pixels: " << m_pixel_mask.size() << std::endl;
    return;
}

// MvtxPixelMask class
//==============================================================================
void MvtxPixelMaskv1::load_from_CDB(std::string calibfile)
{
    if( !m_pixel_mask.empty() ) clear();


    // Load the hot pixel file from the CDB
    std::string database = CDBInterface::instance()->getUrl(calibfile);  // This is specifically for MVTX Hot Pixels
    CDBTTree *cdbttree = new CDBTTree(database);
    int n_total_masked =  cdbttree->GetSingleIntValue("TotalHotPixels");

    // Load the hot pixel map
    std::set<MvtxRawPixelDefs::mvtx_pixelkey> masked_pixels {};
    for (int i = 0; i < n_total_masked; i++)
    {
        int layer = cdbttree->GetIntValue(i, "layer");
        int stave = cdbttree->GetIntValue(i, "stave");
        int chip = cdbttree->GetIntValue(i, "chip");
        int row = cdbttree->GetIntValue(i, "row");
        int col = cdbttree->GetIntValue(i, "col");

        MvtxRawPixelDefs::mvtx_pixelkey this_pixel_key = MvtxRawPixelDefs::gen_pixelkey(layer, stave, chip, row, col);
        masked_pixels.insert(this_pixel_key);
    }

    // Copy the masked pixels to the hot pixel map
    // m_hot_pixel_map.assign(masked_pixels.begin(), masked_pixels.end());
    for (auto it = masked_pixels.begin(); it != masked_pixels.end(); it++){
        m_pixel_mask.push_back(*it);
    }

    return;

}

void MvtxPixelMaskv1::add_pixel(MvtxRawPixelDefs::mvtx_pixelkey key)
{
    if (std::find(m_pixel_mask.begin(), m_pixel_mask.end(), key) == m_pixel_mask.end())
    {
        m_pixel_mask.push_back(key);
    }

    return;
}

void MvtxPixelMaskv1::remove_pixel(MvtxRawPixelDefs::mvtx_pixelkey key)
{
    auto it = std::find(m_pixel_mask.begin(), m_pixel_mask.end(), key);
    if (it != m_pixel_mask.end())
    {
        m_pixel_mask.erase(it);
    }

    return;
}

bool MvtxPixelMaskv1::is_masked(MvtxRawHit* hit) const
{
    MvtxRawPixelDefs::mvtx_pixelkey this_pixel_key = MvtxRawPixelDefs::gen_pixelkey(hit);

    return std::find(m_pixel_mask.begin(), m_pixel_mask.end(), this_pixel_key) != m_pixel_mask.end();
}