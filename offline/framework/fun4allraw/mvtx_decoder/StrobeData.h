#ifndef MVTXDECODER_STROBEDATA_H
#define MVTXDECODER_STROBEDATA_H

#include "GBTWord.h"
#include "InteractionRecord.h"

#include <vector>

namespace mvtx_offline
{
  using mvtx_hit = struct mvtx_hit
  {
    uint8_t chip_id{0xf};
    uint16_t bunchcounter{0xFFFF};
    uint16_t row_pos{0xFFFF};
    uint16_t col_pos{0xFFFF};
  };

  struct StrobeData
  {
    StrobeData(uint64_t orb, uint16_t b)
      : ir(orb, b) {};
    ~StrobeData() = default;

    void clear();

    InteractionRecord ir = {};
    bool hasCDW = false;
    GBTCalibDataWord calWord = {};
    uint32_t detectorField = 0;

    std::vector<mvtx_hit *> hit_vector = {};
  };

}  // namespace mvtx_offline

#endif  // _MVTXDECODER_STROBEDATA
