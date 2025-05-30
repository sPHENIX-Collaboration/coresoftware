#include "MvtxHitMap.h"
#include "MvtxPixelDefs.h"

#include <trackbase/MvtxDefs.h>
#include <trackbase/TrkrDefs.h>

#include <ffarawobjects/MvtxRawHit.h>
#include <ffarawobjects/MvtxRawHitv1.h>

#include <algorithm>
#include <iostream>
#include <set>
#include <string>

// MvtxHitMap class
//==============================================================================
void MvtxHitMap::add_hit(MvtxPixelDefs::pixelkey key, uint32_t nhits)
{
  // Check if the pixel is already in the map
  is_sorted = false;
  auto it = std::find_if(m_pixel_hit_vector.begin(), m_pixel_hit_vector.end(), [key](const pixel_hits_pair_t& element)
                         { return element.first == key; });

  if (it != m_pixel_hit_vector.end())
  {
    // If the pixel is already in the map, increment the hit count
    it->second += nhits;
  }
  else
  {
    // If the pixel is not in the map, add it
    m_pixel_hit_vector.emplace_back(key, nhits);
  }

  return;
}

uint32_t MvtxHitMap::get_nhits(MvtxPixelDefs::pixelkey key) const
{
  // Check if the pixel is in the map
  auto it = std::find_if(m_pixel_hit_vector.begin(), m_pixel_hit_vector.end(), [key](const pixel_hits_pair_t& element)
                         { return element.first == key; });

  if (it != m_pixel_hit_vector.end())
  {
    // If the pixel is in the map, return the hit count
    return it->second;
  }

  // If the pixel is not in the map, return 0
  return 0;
}

MvtxPixelDefs::pixelkey MvtxHitMap::get_most_significant_pixel()
{
  // Find the pixel with the most hits
  if (m_pixel_hit_vector.empty())
  {
    return MvtxPixelDefs::VOID_PIXEL;
  }

  sort_by_hits();
  auto it = m_pixel_hit_vector.begin();
  return it->first;
}

void MvtxHitMap::sort_by_hits()
{
  // Sort the pixel hit vector by hit count
  if (is_sorted)
  {
    return;
  }

  std::sort(m_pixel_hit_vector.begin(), m_pixel_hit_vector.end(), [](const pixel_hits_pair_t& a, const pixel_hits_pair_t& b)
            { return a.second > b.second; });
  is_sorted = true;
  return;
}

uint32_t MvtxHitMap::sum_hits(unsigned int nmasked)
{
  // Sum the hit counts in the pixel hit vector
  uint32_t sum = 0;
  sort_by_hits();
  // sum all hits after the first nmasked pixels
  for (auto it = m_pixel_hit_vector.begin() + nmasked; it != m_pixel_hit_vector.end(); ++it)
  {
    sum += it->second;
  }

  return sum;
}
