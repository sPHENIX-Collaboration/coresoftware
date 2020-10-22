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
//   24 - 32  tracker id
//   16 - 24  layer
//   8  - 16  ladder z id
//   0  -  8  ladder phi id
static const unsigned int kBitShiftLadderZId __attribute__((unused)) = 8;
static const unsigned int kBitShiftLadderPhiId __attribute__((unused)) = 0;
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
TrkrDefs::hitsetkey genHitSetKey(const uint8_t lyr, const uint8_t ladder_z_index, const uint8_t ladder_phi_index);

/**
   * @brief Generate a cluster key from indeces 
   * @param[in] lyr Layer index
   * @param[in] ladder_z_index z index of sensor in ladder
   * @param[in] ladder_phi_ndex phi index of ladder in layer
   * @param[in] clusid Cluster id
   * @param[out] cluskey
   */
TrkrDefs::cluskey genClusKey(const uint8_t lyr, const uint8_t LadderZId, const uint8_t LadderPhiId, const uint32_t clusid);

/**
   * @brief Generate a cluster key using a hitsetkey and cluster id
   * @param[in] hskey hitsetkey
   * @param[in] clusid Cluster id
   * @param[out] cluskey
   */
TrkrDefs::cluskey genClusKey(const TrkrDefs::hitsetkey hskey, const uint32_t clusid);

}  // namespace InttDefs

#endif  //INTT_INTTDEFS_H
