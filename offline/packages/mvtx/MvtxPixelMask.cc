#include "MvtxPixelMask.h"
#include "MvtxPixelDefs.h"

#include <cdbobjects/CDBTTree.h>
#include <ffamodules/CDBInterface.h>

#include <trackbase/MvtxDefs.h>
#include <trackbase/TrkrDefs.h>

#include <ffarawobjects/MvtxRawHit.h>
#include <ffarawobjects/MvtxRawHitv1.h>

#include <algorithm>
#include <iostream>
#include <set>
#include <string>

// MvtxPixelMask class
//==============================================================================
void MvtxPixelMask::load_from_CDB()
{
  if (!m_hot_pixel_map.empty())
  {
    clear();
  }

  // Load the hot pixel file from the CDB
  std::string database = CDBInterface::instance()->getUrl("MVTX_HotPixelMap");  // This is specifically for MVTX Hot Pixels
  CDBTTree* cdbttree = new CDBTTree(database);
  int n_total_masked = cdbttree->GetSingleIntValue("TotalHotPixels");

  // Load the hot pixel map
  std::set<MvtxPixelDefs::pixelkey> masked_pixels{};
  for (int i = 0; i < n_total_masked; i++)
  {
    int layer = cdbttree->GetIntValue(i, "layer");
    int stave = cdbttree->GetIntValue(i, "stave");
    int chip = cdbttree->GetIntValue(i, "chip");
    int row = cdbttree->GetIntValue(i, "row");
    int col = cdbttree->GetIntValue(i, "col");

    MvtxPixelDefs::pixelkey this_pixel_key = MvtxPixelDefs::gen_pixelkey_from_coors(layer, stave, chip, row, col);
    masked_pixels.insert(this_pixel_key);
  }

  // Copy the masked pixels to the hot pixel map
  // m_hot_pixel_map.assign(masked_pixels.begin(), masked_pixels.end());
  for (unsigned long masked_pixel : masked_pixels)
  {
    m_hot_pixel_map.push_back(masked_pixel);
  }

  return;
}

void MvtxPixelMask::add_pixel(MvtxPixelDefs::pixelkey key)
{
  if (std::find(m_hot_pixel_map.begin(), m_hot_pixel_map.end(), key) == m_hot_pixel_map.end())
  {
    m_hot_pixel_map.push_back(key);
  }

  return;
}

void MvtxPixelMask::remove_pixel(MvtxPixelDefs::pixelkey key)
{
  auto it = std::find(m_hot_pixel_map.begin(), m_hot_pixel_map.end(), key);
  if (it != m_hot_pixel_map.end())
  {
    m_hot_pixel_map.erase(it);
  }

  return;
}

void MvtxPixelMask::clear()
{
  m_hot_pixel_map.clear();
  return;
}

bool MvtxPixelMask::is_masked(MvtxRawHit* hit) const
{
  uint8_t layer = hit->get_layer_id();
  uint8_t stave = hit->get_stave_id();
  uint8_t chip = hit->get_chip_id();
  uint16_t row = hit->get_row();
  uint16_t col = hit->get_col();

  // generate hit key and hitset key
  const TrkrDefs::hitkey this_pixel_hitkey = MvtxDefs::genHitKey(col, row);
  const TrkrDefs::hitsetkey this_pixel_hitsetkey = MvtxDefs::genHitSetKey(layer, stave, chip, 0);
  MvtxPixelDefs::pixelkey this_pixel_key = MvtxPixelDefs::gen_pixelkey(this_pixel_hitsetkey, this_pixel_hitkey);

  return std::find(m_hot_pixel_map.begin(), m_hot_pixel_map.end(), this_pixel_key) != m_hot_pixel_map.end();
}