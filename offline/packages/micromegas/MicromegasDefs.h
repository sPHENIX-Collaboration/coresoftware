// Tell emacs that this is a C++ source
// This file is really -*- C++ -*-.

#ifndef MICROMEGAS_MICROMEGASDEFS_H
#define MICROMEGAS_MICROMEGASDEFS_H

/*!
 * \file MicromegasDefs.h
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include <trackbase/TrkrDefs.h>

namespace MicromegasDefs
{

  //* tells the direction along which a given cylinder is segmented
  enum class SegmentationType
  {
    SEGMENTATION_Z,
    SEGMENTATION_PHI
  };

  /*!
   * hitsetkey layout:
   * Micromegas specific lower 16 bits
   * 24 - 32  tracker id
   * 16 - 24  layer
   * 0 - 16 tile id
   */
  static constexpr unsigned int kBitShiftTileId __attribute__((unused)) = 0;

  //! bit shift for hit key
  static constexpr unsigned int kBitShiftStrip __attribute__((unused)) = 0;

  /*!
   * @brief Generate a hitsetkey for the micromegas
   * @param[in] layer Layer index
   * @param[in] tile tile index
   * @param[out] hitsetkey
   *
   * Generate a hitsetkey for the mvtx. The tracker id is known
   * implicitly and used in the function.
   */
  TrkrDefs::hitsetkey genHitSetKey(uint8_t layer, uint8_t tile );

  /*!
   * @brief Get the stave id from hitsetkey
   * @param[in] hitsetkey
   * @param[out] stave id
   s*/
  uint8_t getTileId(TrkrDefs::hitsetkey);

  /*!
   * @brief Generate a hitkey from strip index inside tile
   * @param[in] strip strip index
   */
  TrkrDefs::hitkey genHitKey(uint16_t strip );

  //! get strip from hit key
  uint16_t getStrip(TrkrDefs::hitkey);

  /**
   * @brief Generate a cluster key using a hitsetkey and cluster id
   * @param[in] hskey hitsetkey
   * @param[in] clusid Cluster id
   * @param[out] cluskey
   */
  TrkrDefs::cluskey genClusterKey(TrkrDefs::hitsetkey hskey, uint32_t clusid);
}

#endif
