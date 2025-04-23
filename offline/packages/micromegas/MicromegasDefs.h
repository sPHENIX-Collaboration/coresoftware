// Tell emacs that this is a C++ source
// This file is really -*- C++ -*-.

#ifndef MICROMEGAS_MICROMEGASDEFS_H
#define MICROMEGAS_MICROMEGASDEFS_H

/*!
 * \file MicromegasDefs.h
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include <trackbase/TrkrDefs.h>

#include <array>
#include <cstdint>

namespace MicromegasDefs
{

  //! tells the direction along which a given cylinder is segmented
  enum class SegmentationType: uint8_t
  {
    SEGMENTATION_Z,
    SEGMENTATION_PHI
  };

  //! tells the drift direction for a given micromegas layer
  /*! this is needed for properly implementing transverse diffusion in the layer */
  enum class DriftDirection: uint8_t
  {
    INWARD,
    OUTWARD
  };

  /*!
   * @brief Generate a hitsetkey for the micromegas
   * @param[in] layer Layer index
   * @param[in] tile tile index
   * @param[out] hitsetkey
   *
   * Generate a hitsetkey for the mvtx. The tracker id is known
   * implicitly and used in the function.
   */
  TrkrDefs::hitsetkey genHitSetKey(uint8_t layer, SegmentationType segmentation, uint8_t tile );

  /*!
   * @brief Get the segmentation type from hitsetkey
   * @param[in] hitsetkey
   * @param[out] segmentation
   s*/
  SegmentationType getSegmentationType(TrkrDefs::hitsetkey);

  /*!
   * @brief Get the tile id from hitsetkey
   * @param[in] hitsetkey
   * @param[out] tile id
   s*/
  uint8_t getTileId(TrkrDefs::hitsetkey);

  /*!
   * @brief Generate a hitkey from strip index inside tile
   * @param[in] strip strip index
   */
  TrkrDefs::hitkey genHitKey(uint16_t strip );

  //! get strip from hit key
  uint16_t getStrip(TrkrDefs::hitkey);

  /*!
   * @brief Get the segmentation type from cluster key
   * @param[in] cluskey
   * @param[out] segmentation
   s*/
  SegmentationType getSegmentationType(TrkrDefs::cluskey);

  /*!
   * @brief Get the tile id from cluster key
   * @param[in] cluskey
   * @param[out] tile id
   s*/
  uint8_t getTileId(TrkrDefs::cluskey);

  //! number of TPOT tiles
  static constexpr int m_ntiles = 8;

  //! TPOT packet ids
  /**
   * note: TPOT only uses 2 packets.
   * For early runs (before 07/28/2023) they are 5000 and 5001
   * For later run, (after 07/28/2023) and because, as instructed by Martin, they are 5001 and 5002
   * We keep all 3 values here in order to be able to read both type of runs
   * One might need to change this in the future when 5000 becomes used by a different subsystem
   */
  static constexpr int m_npackets = 3;
  static constexpr int m_npackets_active = 2;
  static constexpr std::array<unsigned int,m_npackets> m_packet_ids = {5000, 5001, 5002};

  //! number of channels per fee board
  static constexpr int m_nchannels_fee = 256;

  //! number of sampa chips per fee board
  static constexpr int m_nsampa_fee = 8;

  //! number of fee boards
  static constexpr int m_nfee = 16;

  //! total number of channels
  static constexpr int m_nchannels_total = m_nfee*m_nchannels_fee;

  //! maximum valid ADC
  static constexpr uint16_t m_adc_max = 1024;

  //! mark invalid ADC values
  static constexpr uint16_t m_adc_invalid = 65000;

  /* see: https://git.racf.bnl.gov/gitea/Instrumentation/sampa_data/src/branch/fmtv2/README.md */
  // TODO: should move to online_distribution
  enum SampaDataType
  {
    HEARTBEAT_T = 0b000,
    TRUNCATED_DATA_T = 0b001,
    TRUNCATED_TRIG_EARLY_DATA_T = 0b011,
    NORMAL_DATA_T = 0b100,
    LARGE_DATA_T = 0b101,
    TRIG_EARLY_DATA_T = 0b110,
    TRIG_EARLY_LARGE_DATA_T = 0b111,
  };

}

#endif
