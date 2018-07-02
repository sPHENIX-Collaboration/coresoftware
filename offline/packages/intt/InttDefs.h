/**
 * @file intt/InttDefs.h
 * @author D. McGlinchey
 * @date June 2018
 * @brief Utility functions for INTT
 */
#ifndef INTT_INTTDEFS_H
#define INTT_INTTDEFS_H

#include <trackbase/TrkrDefs.h>

/**
 * @brief Utility functions for INTT
 *
 * Contains the functions for manipulating the various keys
 * used by the intt for hits, hit sets, and clusters
 */
namespace InttDefs
{

#ifndef __CINT__
  // hitsetkey layout:
  //  Intt specific lower 16 bits
  //   24 - 32  tracker id
  //   16 - 24  layer
  //   8  - 16  ladder id
  //   0  -  8  sensor id
  static const unsigned int kBitShiftLadderId __attribute__((unused)) = 8;
  static const unsigned int kBitShiftSensorId __attribute__((unused)) = 0;

#endif // __CINT__

  /**
   * @brief Get the ladder id from hitsetkey
   * @param[in] hitsetkey
   * @param[out] ladder id
   */
  uint8_t getLadderId(TrkrDefs::hitsetkey key);

  /**
   * @brief Get the ladder id from cluskey
   * @param[in] cluskey
   * @param[out] ladder id
   */
  uint8_t getLadderId(TrkrDefs::cluskey key);

  /**
   * @brief Get the sensor id from hitsetkey
   * @param[in] hitsetkey
   * @param[out] sensor id
   */
  uint8_t getSensorId(TrkrDefs::hitsetkey key);

  /**
   * @brief Get the sensor id from cluskey
   * @param[in] cluskey
   * @param[out] sensor id
   */
  uint8_t getSensorId(TrkrDefs::cluskey key);

  /**
   * @brief Generate a hitkey from a strip id
   * @param[in] strip Strip id
   * @param[out] hitkey
   */
  TrkrDefs::hitkey genHitKey(const uint32_t strip);

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
  TrkrDefs::hitsetkey genHitSetKey(const uint8_t lyr, const uint8_t ladder, const uint8_t sensor);

  /**
   * @brief Generate a cluster key from indeces 
   * @param[in] lyr Layer index
   * @param[in] ladder Ladder index
   * @param[in] sensor Sensor index
   * @param[in] clusid Cluster id
   * @param[out] cluskey
   */
  TrkrDefs::cluskey genClusKey(const uint8_t lyr, const uint8_t ladder, const uint8_t sensor, const uint32_t clusid);

  /**
   * @brief Generate a cluster key using a hitsetkey and cluster id
   * @param[in] hskey hitsetkey
   * @param[in] clusid Cluster id
   * @param[out] cluskey
   */
  TrkrDefs::cluskey genClusKey(const TrkrDefs::hitsetkey hskey, const uint32_t clusid);


}

#endif  //INTT_INTTDEFS_H
