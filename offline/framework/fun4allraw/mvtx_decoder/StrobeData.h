#ifndef MVTXDECODER_STROBEDATA_H
#define MVTXDECODER_STROBEDATA_H

#include "InteractionRecord.h"
#include "GBTWord.h"

#include <vector>

namespace mvtx
{
  typedef struct mvtx_hit
  {
    uint8_t chip_id = 0xf;
    uint16_t bunchcounter = 0xFFFF;
    uint16_t row_pos = 0xFFFF;
    uint16_t col_pos = 0xFFFF;
  } mvtx_hit;

  struct StrobeData
  {
    StrobeData(uint64_t orb, uint16_t b) : ir(orb, b) {};
    ~StrobeData();

    void clear();

    InteractionRecord ir = {};
    bool hasCDW = false;
    GBTCalibDataWord calWord = {};
    uint32_t detectorField = 0;

    std::vector<mvtx_hit *> hit_vector = {};
  };

} // namespace mvtx

#endif // _MVTXDECODER_STROBEDATA
