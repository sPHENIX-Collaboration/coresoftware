/**
 * @file mvtx/MvtxDefs.h
 * @author D. McGlinchey
 * @date June 2018
 * @brief Utility functions for MVTX
 */
#ifndef MVTX_MVTXDEFUTIL_H
#define MVTX_MVTXDEFUTIL_H

#include "TrkrDefs.h"

#include <cstdint>  // for uint8_t, uint16_t, uint32_t

/**
 * @brief Utility functions for MVTX
 *
 * Contains the functions for manipulating the various keys
 * used by the mvtx for hits, hit sets, and clusters
 */
namespace MvtxDefs
{
  // max values for col and row index in chip
  static constexpr uint16_t MAXCOL __attribute__((unused)) = 1024;
  static constexpr uint16_t MAXROW __attribute__((unused)) = 512;

  /**
   * @brief Get the stave id from hitsetkey
   * @param[in] hitsetkey
   * @param[out] stave id
   */
  uint8_t getStaveId(TrkrDefs::hitsetkey key);

  /**
   * @brief Get the stave id from cluskey
   * @param[in] cluskey
   * @param[out] stave id
   */
  uint8_t getStaveId(TrkrDefs::cluskey key);

  /**
   * @brief Get the chip id from hitsetkey
   * @param[in] hitsetkey
   * @param[out] chip id
   */
  uint8_t getChipId(TrkrDefs::hitsetkey key);

  /**
   * @brief Get the chip id from cluskey
   * @param[in] cluskey
   * @param[out] chip id
   */
  uint8_t getChipId(TrkrDefs::cluskey key);

  /**
   * @brief Get the chip id from hitsetkey
   * @param[in] hitsetkey
   * @param[out] chip id
   */
  int getStrobeId(TrkrDefs::hitsetkey key);

  /**
   * @brief Get the strobe id from hitsetkey
   * @param[in] hitsetkey
   * @param[out] strobe id
   */
  int getStrobeId(TrkrDefs::cluskey key);

  /**
   * @brief Get the strobe id from the cluskey
   * @param[in] cluskey
   * @param[out] strobe id
   */
  uint16_t getCol(TrkrDefs::hitkey key);

  /**
   * @brief Get the row index from hitkey
   * @param[in] hitkey
   * @param[out] row index
   */
  uint16_t getRow(TrkrDefs::hitkey key);

  /**
   * @brief Generate a hitkey from a pixels column and row index
   * @param[in] col Column index
   * @param[in] row Row index
   * @param[out] hitkey
   */
  TrkrDefs::hitkey genHitKey(const uint16_t col, const uint16_t row);

  /**
   * @brief Generate a hitsetkey for the mvtx
   * @param[in] lyr Layer index
   * @param[in] stave Stave index
   * @param[in] chip Chip index
   * @param[out] hitsetkey
   *
   * Generate a hitsetkey for the mvtx. The tracker id is known
   * implicitly and used in the function.
   */
  TrkrDefs::hitsetkey genHitSetKey(const uint8_t lyr, const uint8_t stave, const uint8_t chip, const int strobe);

  /**
   * @brief Generate a cluster key from indeces
   * @param[in] lyr Layer index
   * @param[in] stave Stave index
   * @param[in] chip Chip index
   * @param[in] crossing bunch crossing
   * @param[in] clusid Cluster id
   * @param[out] cluskey
   */
  TrkrDefs::cluskey genClusKey(const uint8_t lyr, const uint8_t stave, const uint8_t chip, const int strobe, const uint32_t clusid);

  /**
   * @brief Zero the strobe bits in the hitsetkey
   * @param[in] hskey hitsetkey
   * @param[out] hitsetkey with strobe bits set to zero
   */
  TrkrDefs::hitsetkey resetStrobeHitSetKey(const TrkrDefs::hitsetkey hitsetkey);

}  // namespace MvtxDefs

#endif  // MVTX_MVTXDEFUTIL_H
