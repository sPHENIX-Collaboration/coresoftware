/**
 * @file intt/InttDefs.h
 * @author D. McGlinchey
 * @date June 2018
 * @brief Utility functions for INTT
 */
#ifndef INTT_INTTDEFS_H
#define INTT_INTTDEFS_H

#include <trackbase/TrkrDefs.h>

#include <cstdint>              // for uint8_t, uint16_t, uint32_t

/**
 * @brief Utility functions for INTT
 *
 * Contains the functions for manipulating the various keys
 * used by the intt for hits, hit sets, and clusters
 */
namespace InttDefs
{
// hitsetkey layout:
//  Intt specific lower 16 bits
//   24 - 31  tracker id  // 8 bits
//   16 - 23  layer       // 8 bits
//   14 - 15  ladder z id   // 2 bits
//   10 - 13       ladder phi id  // 4 bits
//   0 - 9     time bucket  // 10 bits

// old
//   8  - 16  ladder z id
//   0  -  8  ladder phi id

static const unsigned int kBitShiftTimeBucketIdOffset __attribute__((unused)) = 0;
static const unsigned int kBitShiftTimeBucketIdWidth __attribute__((unused)) = 10;
static const unsigned int kBitShiftLadderPhiIdOffset __attribute__((unused)) = 10;
static const unsigned int kBitShiftLadderPhiIdWidth __attribute__((unused)) = 4;
static const unsigned int kBitShiftLadderZIdOffset __attribute__((unused)) = 14;
static const unsigned int kBitShiftLadderZIdWidth __attribute__((unused)) = 2;

// bit shift for hitkey
static const unsigned int kBitShiftCol __attribute__((unused)) = 16;
static const unsigned int kBitShiftRow __attribute__((unused)) = 0;

/**
   * @brief Get the ladder id from hitsetkey
   * @param[in] hitsetkey
   * @param[out] ladder id
   */
uint8_t getLadderZId(TrkrDefs::hitsetkey key);

/**
   * @brief Get the ladder id from cluskey
   * @param[in] cluskey
   * @param[out] ladder id
   */
uint8_t getLadderZId(TrkrDefs::cluskey key);

/**
   * @brief Get the sensor id from hitsetkey
   * @param[in] hitsetkey
   * @param[out] sensor id
   */
uint8_t getLadderPhiId(TrkrDefs::hitsetkey key);

/**
   * @brief Get the sensor id from cluskey
   * @param[in] cluskey
   * @param[out] sensor id
   */

uint8_t getLadderPhiId(TrkrDefs::cluskey key);

/**
   * @brief Generate a hitkey from a strip id
   * @param[in] strip Strip id
   * @param[out] hitkey
   */

uint16_t getTimeBucketId(TrkrDefs::hitsetkey key);

/**
   * @brief Get the time bucket id from the hitsetkey
   * @param[in] hitsetkey
   * @param[out] time bucket id
   */

uint16_t getTimeBucketId(TrkrDefs::cluskey key);

/**
   * @brief Get the time bucket id from the cluskey
   * @param[in] cluskey
   * @param[out] time bucket id
   */

/**
   * @brief Get the column index from hitkey
   * @param[in] hitkey
   * @param[out] column index
   */
uint16_t getCol(TrkrDefs::hitkey key);

/**
   * @brief Get the row index from hitkey
   * @param[in] hitkey
   * @param[out] row index
   */
uint16_t getRow(TrkrDefs::hitkey key);

 TrkrDefs::hitkey genHitKey(const uint16_t col, const uint16_t row);

/**
   * @brief Generate a hitsetkey for the intt
   * @param[in] lyr Layer index
   * @param[in] ladder Ladder index
   * @param[in] sensor Sensor index
   * @param[out] hitsetkey
   *
   * Generate a hitsetkey for the intt. The tracker id is known
   * implicitly and used in the function.
   */
 TrkrDefs::hitsetkey genHitSetKey(const uint8_t lyr, const uint8_t ladder_z_index, const uint8_t ladder_phi_index, const uint16_t time_bucket);

/**
   * @brief Generate a cluster key from indices 
   * @param[in] lyr Layer index
   * @param[in] ladder_z_index z index of sensor in ladder
   * @param[in] ladder_phi_index phi index of ladder in layer
   * @param[in] time bunch index of time bucket
   * @param[out] hitsetkey
   */

 TrkrDefs::cluskey genClusKey(const uint8_t lyr, const uint8_t LadderZId, const uint8_t LadderPhiId, const uint16_t time_bucket, const uint32_t clusid);

/**
   * @brief Generate a cluster key using a hitsetkey and cluster id
   * @param[in] hskey hitsetkey
   * @param[in] clusid Cluster id
   * @param[out] cluskey
   */
TrkrDefs::cluskey genClusKey(const TrkrDefs::hitsetkey hskey, const uint32_t clusid);

/**
   * @brief Zero the crossing bits in a copy of the  hitsetkey
   * @param[in] hitsetkey
   * @param[out] hitsetkey
   */
TrkrDefs::hitsetkey resetCrossingHitSetKey(const TrkrDefs::hitsetkey hitsetkey);

}  // namespace InttDefs

#endif  //INTT_INTTDEFS_H
