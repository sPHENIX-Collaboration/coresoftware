#ifndef MVTX_MVTXPIXELMASK_H
#define MVTX_MVTXPIXELMASK_H

#include "MvtxPixelDefs.h"

#include <climits>
#include <memory>
#include <vector>

class MvtxRawHit;

class MvtxPixelMask
{
 public:
  MvtxPixelMask() {}
  ~MvtxPixelMask() { clear(); }

  typedef std::vector<MvtxPixelDefs::pixelkey> hot_pixel_map_t;

  void load_from_CDB();

  void add_pixel(MvtxPixelDefs::pixelkey key);
  void remove_pixel(MvtxPixelDefs::pixelkey key);

  void clear();

  bool is_masked(MvtxRawHit* hit) const;

  hot_pixel_map_t get_hot_pixel_map() const { return m_hot_pixel_map; }

 private:
  hot_pixel_map_t m_hot_pixel_map{};
};

#endif  // MVTX_MVTXPIXELMASK_H