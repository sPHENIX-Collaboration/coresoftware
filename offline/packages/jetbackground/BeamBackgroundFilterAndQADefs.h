/// ===========================================================================
/*! \file    BeamBackgroundFilterAndQADefs.h
 *  \authors Derek Anderson
 *  \date    11.01.2024
 *
 *  A namespace to hold various definitions
 *  useful across the BeamBackgroundFilterAndQA
 *  module and its various filters.
 */
/// ===========================================================================

#ifndef BEAMBACKGROUNDFILTERANDQADEFS_H
#define BEAMBACKGROUNDFILTERANDQADEFS_H

// calobase libraries
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>

// c++ utilities
#include <array>

// ============================================================================
//! Misc beam background filter and QA definitions
// ============================================================================
namespace BeamBackgroundFilterAndQADefs
{

  // ==========================================================================
  //! Event status codes
  // ==========================================================================
  /*! This enumerates the outcomes of applying a filter:
   *    Evt     = null, no filter applied yet
   *    NoBkgd  = filter finds there is no beam background
   *    HasBkgd = filter finds there is beam background
   */
  enum Status
  {
    Evt,
    NoBkgd,
    HasBkgd
  };

  // ==========================================================================
  // Helper struct to scrape info from TowerInfo
  // ==========================================================================
  /*! This is a lightweight struct to scrape only the info relevant to various
   *  filters from a TowerInfo object.
   *
   *  An example usage is in the StreakSidebandFilter algorithm, where its used
   *  to build a 2D array of tower energies & status (see below) for quick
   *  lookup.
   */
  struct Tower
  {
    // members
    uint8_t status = 0;
    double energy = -1.;
    bool isGood = false;
    //! grab info from a TowerInfo object
    void SetInfo(TowerInfo* info)
    {
      status = info->get_status();
      energy = info->get_energy();
      isGood = info->get_isGood();
      return;
    }

    //! reset values
    void Reset()
    {
      status = 0;
      energy = -1.;
      return;
    }

  };  // end Tower

  // ==========================================================================
  //! Helper type for building (eta, phi) maps of towers
  // ==========================================================================
  /*! This is a lightweight struct to organize the above Tower structs into a
   *  handy (eta, phi) map. See the StreakSidebandFilter algorithm for an
   *  example usage.
   */
  template <std::size_t H, std::size_t F>
  struct TowerMap
  {
    // members
    std::array<std::array<Tower, F>, H> towers;

    //! build array
    void Build(TowerInfoContainer* container)
    {
      for (std::size_t iTwr = 0; iTwr < container->size(); ++iTwr)
      {
        const int32_t key = container->encode_key(iTwr);
        const int32_t iEta = container->getTowerEtaBin(key);
        const int32_t iPhi = container->getTowerPhiBin(key);
        towers.at(iEta).at(iPhi).SetInfo(container->get_tower_at_channel(iTwr));
      }
      return;
    }

    //! reset
    void Reset()
    {
      for (auto row : towers)
      {
        for (auto tower : row)
        {
          tower.Reset();
        }
      }
      return;
    }

  };  // end TowerArray

  // --------------------------------------------------------------------------
  //! Maps for specific calorimeters
  // --------------------------------------------------------------------------
  typedef TowerMap<96, 256> EMCalMap;
  typedef TowerMap<24, 64> IHCalMap;
  typedef TowerMap<24, 64> OHCalMap;

}  // namespace BeamBackgroundFilterAndQADefs

#endif

// end ========================================================================
