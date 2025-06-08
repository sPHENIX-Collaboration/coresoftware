/*!
 * \file MvtxHitMap.h
 * \brief Holds pixel hit info until processing
 * \author Tanner Mengel <tmengel@bnl.gov>
 * \version $Version: 2.0.1 $
 * \date $Date: 05/23/2025.
 */

#ifndef MVTX_MVTXHITMAP_H
#define MVTX_MVTXHITMAP_H

#include "MvtxPixelDefs.h"

#include <cstdint>
#include <utility>
#include <vector>

class MvtxRawHit;
class MvtxHitMap
{
 public:
  MvtxHitMap() = default;
  ~MvtxHitMap() { clear(); }

  typedef std::pair<MvtxPixelDefs::pixelkey, uint32_t> pixel_hits_pair_t;
  typedef std::vector<MvtxHitMap::pixel_hits_pair_t> pixel_hit_vector_t;

  void add_hit(MvtxPixelDefs::pixelkey key, uint32_t nhits = 1);
  void clear()
  {
    m_pixel_hit_vector.clear();
    is_sorted = false;
  }

  int npixels() const { return m_pixel_hit_vector.size(); }

  uint32_t get_nhits(MvtxPixelDefs::pixelkey key) const;
  const MvtxHitMap::pixel_hit_vector_t &get_pixel_hit_vector() const { return m_pixel_hit_vector; }
  MvtxPixelDefs::pixelkey get_most_significant_pixel();

  uint32_t sum_hits(unsigned int nmasked = 0);

 private:
  MvtxHitMap::pixel_hit_vector_t m_pixel_hit_vector{};

  void sort_by_hits();
  bool is_sorted{false};
};

#endif  // MVTX_MVTXHITMAP_H
