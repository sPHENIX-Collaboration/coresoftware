/**
 * @file outertracker/OuterTrackerDefs.h
 * @author D. McGlinchey
 * @date June 2018
 * @brief Utility functions for OuterTracker
 */
#ifndef OTRACK_OUTERTRACKERDEFS_H
#define OTRACK_OUTERTRACKERDEFS_H

#include <trackbase/TrkrDefs.h>

#include <cstdint>              // for uint8_t, uint16_t, uint32_t

/**
 * @brief Utility functions for OUTERTRACKER
 *
 * Contains the functions for manipulating the various keys
 * used by the OuterTracker for hits, hit sets, and clusters
 */
namespace OuterTrackerDefs
{
#ifndef __CINT__
// hitsetkey layout:
//  OuterTracker specific lower 16 bits
//   24 - 32  tracker id
//   16 - 24  layer

// bit shift for hitkey
static const unsigned int kBitShiftCol __attribute__((unused)) = 16;
static const unsigned int kBitShiftRow __attribute__((unused)) = 0;

#endif  // __CINT__


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
   * @brief Generate a hitsetkey for the OuterTracker
   * @param[in] lyr Layer index
   * @param[in] ladder Ladder index
   * @param[in] sensor Sensor index
   * @param[out] hitsetkey
   *
   * Generate a hitsetkey for the OuterTracker. The tracker id is known
   * implicitly and used in the function.
   */
TrkrDefs::hitsetkey genHitSetKey(const uint8_t lyr);

/**
   * @brief Generate a cluster key from indices 
   * @param[in] lyr Layer index
   * @param[in] clusid Cluster id
   * @param[out] cluskey
   */
TrkrDefs::cluskey genClusKey(const uint8_t lyr, const uint32_t clusid);

/**
   * @brief Generate a cluster key using a hitsetkey and cluster id
   * @param[in] hskey hitsetkey
   * @param[in] clusid Cluster id
   * @param[out] cluskey
   */
TrkrDefs::cluskey genClusKey(const TrkrDefs::hitsetkey hskey, const uint32_t clusid);

}  // namespace OuterTrackerDefs

#endif  //OTRACK_OUTERTRACKERDEFS_H
