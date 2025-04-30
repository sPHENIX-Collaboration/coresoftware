#ifndef MVTXRAWDEFS_H
#define MVTXRAWDEFS_H

#include <array>
#include <cstdint>
#include <map>

namespace MvtxRawDefs
{

  static constexpr uint8_t nMvtxLayers = 3;
  static constexpr uint8_t nGbtPerStave = 3;
  static constexpr uint8_t nChipsPerGbt = 3;
  static constexpr uint8_t nStavesPerFelix = 8;

  static const std::array<uint8_t, nMvtxLayers> firstStaveIndex = {{0, 12, 28}};

  static const std::array<std::array<uint8_t, nChipsPerGbt>, nGbtPerStave> gbtChipId_to_staveChipId = {{{{0, 1, 2}}, {{3, 4, 5}}, {{6, 7, 8}}}};

  // 07/03/2024
  // {staveIndex[0-47], { FLX[0-5], EndPoint[0-1] } }
  const std::map<uint8_t, std::pair<uint8_t, uint8_t>> stave_felix_map =
      {
          {0, {2, 0}},
          {1, {2, 0}},
          {2, {2, 1}},
          {3, {0, 0}},
          {4, {0, 0}},
          {5, {0, 1}},
          {6, {4, 0}},
          {7, {4, 0}},
          {8, {4, 1}},
          {9, {3, 0}},
          {10, {3, 0}},
          {11, {3, 1}},
          {12, {1, 0}},
          {13, {1, 0}},
          {14, {1, 1}},
          {15, {1, 1}},
          {16, {2, 0}},
          {17, {2, 1}},
          {18, {2, 1}},
          {19, {2, 1}},
          {20, {5, 0}},
          {21, {5, 0}},
          {22, {5, 1}},
          {23, {5, 1}},
          {24, {4, 0}},
          {25, {4, 1}},
          {26, {4, 1}},
          {27, {4, 1}},
          {28, {4, 0}},
          {29, {0, 0}},
          {30, {0, 0}},
          {31, {0, 1}},
          {32, {0, 1}},
          {33, {3, 0}},
          {34, {1, 0}},
          {35, {1, 0}},
          {36, {1, 1}},
          {37, {1, 1}},
          {38, {2, 0}},
          {39, {3, 0}},
          {40, {3, 1}},
          {41, {3, 1}},
          {42, {3, 1}},
          {43, {0, 1}},
          {44, {5, 0}},
          {45, {5, 0}},
          {46, {5, 1}},
          {47, {5, 1}},
  };

  using linkId_t = struct linkId
  {
    uint32_t layer{0xFF};
    uint32_t stave{0xFF};
    uint32_t gbtid{0xFF};
  };

  uint8_t getStaveIndex(const uint8_t& lyrId, const uint8_t& stvId);
  std::pair<uint8_t, uint8_t> const& get_flx_endpoint(const uint8_t& lyrId, const uint8_t& stvId);

  linkId_t decode_feeid(uint16_t feeid);

  float getStrobeLength(const int& runNumber);

}  // namespace MvtxRawDefs

#endif  // MVTXRAWDEFS_H
