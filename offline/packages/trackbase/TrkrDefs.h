/**
 * @file trackbase/TrkrDefs.h
 * @author D. McGlinchey
 * @date June 2018
 * @brief Namespace with Trkr key types and utility functions
 */
#ifndef TRACKBASE_TRKRDEFS_H
#define TRACKBASE_TRKRDEFS_H

#include <cstdint>
#include <iostream>
#include <map>
#include <string>

/**
 * @brief Define a namespace for Trkr typedefs
 */
namespace TrkrDefs
{
  [[maybe_unused]] static constexpr double EdepScaleFactor = 0.25;
  [[maybe_unused]] static constexpr double MvtxEnergyScaleup = 5.0e8;
  [[maybe_unused]] static constexpr double InttEnergyScaleup = 5.0e7;

  /// Key types
  using hitkey = uint32_t;      // 32 bit TrkrHit key type
  using hitsetkey = uint32_t;   // 32 bit TrkrHitSet key type
  using cluskey = uint64_t;     // 64 but TrkrCluster id type
  using clushitkey = uint32_t;  // 32 bit hit id type in TrkrCluster
  using subsurfkey = uint16_t;  // 16 bit sub surface key type

  /// Max values for keys (used as defaults or invalid values)
  [[maybe_unused]] static constexpr hitkey HITKEYMAX = UINT32_MAX;
  [[maybe_unused]] static constexpr hitsetkey HITSETKEYMAX = UINT32_MAX;
  [[maybe_unused]] static constexpr cluskey CLUSKEYMAX = UINT64_MAX;
  [[maybe_unused]] static constexpr clushitkey CLUSHITKEYMAX = UINT32_MAX;
  [[maybe_unused]] static constexpr subsurfkey SUBSURFKEYMAX = UINT16_MAX;

  /// Enumeration for tracker id to easily maintain consistency
  enum TrkrId
  {
    mvtxId = 0,
    inttId = 1,
    tpcId = 2,
    micromegasId = 3
  };

  //! Standard names for trackers
  static const std::map<TrkrId, std::string> TrkrNames =
  {
      {mvtxId, "MVTX"},
      {inttId, "INTT"},
      {tpcId, "TPC"},
      {micromegasId, "MICROMEGAS"}
  };

  /// Print the bits for each key type
  void printBits(const TrkrDefs::hitsetkey key, std::ostream& os = std::cout);
  void printBits(const TrkrDefs::cluskey key, std::ostream& os = std::cout);
  // void print_bits(const TrkrDefs::hitkey key, std::ostream& os = std::cout);

  /// Get the tracker ID from either key type
  uint8_t getTrkrId(const TrkrDefs::hitsetkey key);
  uint8_t getTrkrId(const TrkrDefs::cluskey key);

  /// Get the layer number from either key type
  uint8_t getLayer(const TrkrDefs::hitsetkey key);
  uint8_t getLayer(const TrkrDefs::cluskey key);

  /// Get the lower 32 bits for cluster keys only
  uint32_t getClusIndex(const TrkrDefs::cluskey key);

  /// generate the common upper 16 bits for hitsetkey
  TrkrDefs::hitsetkey genHitSetKey(const TrkrDefs::TrkrId trkrId, const uint8_t lyr);

  /// generate cluster key from hitset key and cluster index
  TrkrDefs::cluskey genClusKey(const TrkrDefs::hitsetkey hskey, const uint32_t clusid);

  /// Get the upper 32 bits from cluster keys
  uint32_t getHitSetKeyFromClusKey(const TrkrDefs::cluskey key);

  /// Get a valid low / hi range for hitsetkey given tracker id & layer
  TrkrDefs::hitsetkey getHitSetKeyLo(const TrkrDefs::TrkrId trkrId);
  TrkrDefs::hitsetkey getHitSetKeyHi(const TrkrDefs::TrkrId trkrId);
  TrkrDefs::hitsetkey getHitSetKeyLo(const TrkrDefs::TrkrId trkrId, const uint8_t lyr);
  TrkrDefs::hitsetkey getHitSetKeyHi(const TrkrDefs::TrkrId trkrId, const uint8_t lyr);

  /// Get a valid low / hi range for cluskey given tracker id & layer
  TrkrDefs::cluskey getClusKeyLo(const TrkrDefs::TrkrId trkrId);
  TrkrDefs::cluskey getClusKeyHi(const TrkrDefs::TrkrId trkrId);
  TrkrDefs::cluskey getClusKeyLo(const TrkrDefs::TrkrId trkrId, const uint8_t lyr);
  TrkrDefs::cluskey getClusKeyHi(const TrkrDefs::TrkrId trkrId, const uint8_t lyr);

  [[maybe_unused]] static constexpr unsigned int kBitShiftPhiElement = 8;  // sector
  [[maybe_unused]] static constexpr unsigned int kBitShiftZElement = 0;    // side

  uint8_t getPhiElement(TrkrDefs::hitsetkey key);  // sector
  uint8_t getZElement(TrkrDefs::hitsetkey key);    // side
  uint8_t getPhiElement(TrkrDefs::cluskey key);    // sector
  uint8_t getZElement(TrkrDefs::cluskey key);      // side

}  // namespace TrkrDefs

#endif  // TRACKBASE_TRKRDEFS_H
