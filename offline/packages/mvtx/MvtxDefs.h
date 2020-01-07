/**
 * @file mvtx/MvtxDefs.h
 * @author D. McGlinchey
 * @date June 2018
 * @brief Utility functions for MVTX
 */
#ifndef MVTX_MVTXDEFUTIL_H
#define MVTX_MVTXDEFUTIL_H

#include <trackbase/TrkrDefs.h>

#include <cstdint>              // for uint8_t, uint16_t, uint32_t

/**
 * @brief Utility functions for MVTX
 *
 * Contains the functions for manipulating the various keys
 * used by the mvtx for hits, hit sets, and clusters
 */
namespace MvtxDefs
{
#if !defined(__CINT__) || defined(__CLING__)
// hitsetkey layout:
//  Mvtx specific lower 16 bits
//   24 - 32  tracker id
//   16 - 24  layer
//   8  - 16  stave id
//   0  -  8  chip id
static const unsigned int kBitShiftStaveId __attribute__((unused)) = 8;
static const unsigned int kBitShiftChipId __attribute__((unused)) = 0;

// bit shift for hitkey
static const unsigned int kBitShiftCol __attribute__((unused)) = 16;
static const unsigned int kBitShiftRow __attribute__((unused)) = 0;

// max values for col and row index in chip
static const uint16_t MAXCOL __attribute__((unused)) = 1024;
static const uint16_t MAXROW __attribute__((unused)) = 512;

#endif  // __CINT__

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
TrkrDefs::hitsetkey genHitSetKey(const uint8_t lyr, const uint8_t stave, const uint8_t chip);

/**
   * @brief Generate a cluster key from indeces 
   * @param[in] lyr Layer index
   * @param[in] stave Stave index
   * @param[in] chip Chip index
   * @param[in] clusid Cluster id
   * @param[out] cluskey
   */
TrkrDefs::cluskey genClusKey(const uint8_t lyr, const uint8_t stave, const uint8_t chip, const uint32_t clusid);

/**
   * @brief Generate a cluster key using a hitsetkey and cluster id
   * @param[in] hskey hitsetkey
   * @param[in] clusid Cluster id
   * @param[out] cluskey
   */
TrkrDefs::cluskey genClusKey(const TrkrDefs::hitsetkey hskey, const uint32_t clusid);

}  // namespace MvtxDefs

#endif  //MVTX_MVTXDEFUTIL_H
