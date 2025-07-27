/// ===========================================================================
/*! \file    StreakSidebandFilter.h
 *  \authors Hanpu Jiang, Derek Anderson
 *  \date    11.01.2024
 *
 *  Part of the BeamBackgroundFilterAndQA module, this
 *  filter returns true if event contains a "streaK"
 *  as defined by Hanpu Jiang's streak identification
 *  algorithm.
 */
/// ===========================================================================

#ifndef STREAKSIDEBANDFILTER_H
#define STREAKSIDEBANDFILTER_H

// module components
#include "BaseBeamBackgroundFilter.h"
#include "BeamBackgroundFilterAndQADefs.h"

// c++ utilities
#include <array>
#include <cstddef>
#include <cstdint>
#include <string>

// forward declarations
class PHCompositeNode;
class TowerInfoContainer;

// ============================================================================
//! Identify streaks via sidebands
// ============================================================================
/*! A beam background filter which identifies streaks
 *  in the OHCal by comparing streak candidates vs.
 *  their sidebands, i.e. adjacent phi slices.
 */
class StreakSidebandFilter : public BaseBeamBackgroundFilter
{
 public:
  // ========================================================================
  //! User options for filter
  // ========================================================================
  struct Config
  {
    int verbosity{0};
    bool debug{true};
    float minStreakTwrEne{0.6};
    float maxAdjacentTwrEne{0.06};
    uint32_t minNumTwrsInStreak{5};
    std::string inNodeName{"TOWERINFO_CALIB_HCALOUT"};
  };

  // ctor/dtor
  StreakSidebandFilter(const std::string& name = "StreakSideband");
  StreakSidebandFilter(const Config& config, const std::string& name = "StreakSideband");
  virtual ~StreakSidebandFilter() = default;

  // inherited methods
  bool ApplyFilter(PHCompositeNode* topNode) override;
  void BuildHistograms(const std::string& module, const std::string& tag = "") override;

 private:
  // inherited methods
  void GrabNodes(PHCompositeNode* topNode) override;

  // filter-specific methods
  bool IsTowerNotStreaky(const BeamBackgroundFilterAndQADefs::Tower& tower) const;
  bool IsNeighborNotStreaky(const BeamBackgroundFilterAndQADefs::Tower& tower) const;

  ///! input node
  TowerInfoContainer* m_ohContainer{nullptr};

  ///! no. of streaks per event in ohcal
  std::array<std::size_t, 64> m_ohNumStreak{};

  ///! tower info (eta, phi) map
  BeamBackgroundFilterAndQADefs::OHCalMap m_ohMap;

  ///! configuration
  Config m_config;

};  // end StreakSidebandFilter

#endif

// end ========================================================================
