/// ===========================================================================
/*! \file    StreakSidebandFilter.cc
 *  \authors Hanpu Jiang, Derek Anderson
 *  \date    11.01.2024
 *
 *  Part of the BeamBackgroundFilterAndQA module, this
 *  filter returns true if event contains a "streaK"
 *  as defined by Hanpu Jiang's streak identification
 *  algorithm.
 */
/// ===========================================================================

// module components
#include "StreakSidebandFilter.h"

// calo base
#include <calobase/TowerInfoContainer.h>

// phool libraries
#include <phool/getClass.h>

// root libraries
#include <TH1.h>
#include <TH2.h>

// c++ utiilites
#include <algorithm>
#include <iostream>
#include <map>
#include <memory>
#include <vector>

// ctor/dtor ==================================================================

// ----------------------------------------------------------------------------
//! Default ctor
// ----------------------------------------------------------------------------
StreakSidebandFilter::StreakSidebandFilter(const std::string& name)
  : BaseBeamBackgroundFilter(name)

{
}  // end ctor()

// ----------------------------------------------------------------------------
//! ctor accepting config struct
// ----------------------------------------------------------------------------
StreakSidebandFilter::StreakSidebandFilter(const Config& config, const std::string& name)
  : BaseBeamBackgroundFilter(name)
  , m_config(config)
{
}  // end ctor(Config&)

// public methods =============================================================

// ----------------------------------------------------------------------------
// Apply filter to check for beam background or not
// ----------------------------------------------------------------------------
bool StreakSidebandFilter::ApplyFilter(PHCompositeNode* topNode)
{
  // print debug message
  if (m_config.debug && (m_config.verbosity > 2))
  {
    std::cout << "StreakSidebandFilter::ApplyFilter() Checking if streak found in OHCal via their sidebands" << std::endl;
  }

  // grab input node
  GrabNodes(topNode);

  // build tower map
  m_ohMap.Reset();
  m_ohMap.Build(m_ohContainer);

  // reset number of streaky towers per eta bin
  std::fill(m_ohNumStreak.begin(), m_ohNumStreak.end(), 0);

  // lambdas to get phi +- 1 neighbors
  auto getAdjacentUp = [this](const std::size_t phi)
  { return (phi + 1) % m_ohMap.towers.front().size(); };
  auto getAdjacentDown = [this](const std::size_t phi)
  { return (phi == 0) ? m_ohMap.towers.front().size() : (phi - 1); };

  // loop over tower (eta, phi) map to find streaks
  for (std::size_t iPhi = 0; iPhi < m_ohMap.towers.front().size(); ++iPhi)
  {
    for (std::size_t iEta = 0; iEta < m_ohMap.towers.size(); ++iEta)
    {
      // check if tower is a candidate for being in a streak
      const bool isNotStreak = IsTowerNotStreaky(m_ohMap.towers[iEta][iPhi]);
      if (isNotStreak)
      {
        continue;
      }

      // grab adjacent towers
      const std::size_t iUp = getAdjacentUp(iPhi);
      const std::size_t iDown = getAdjacentDown(iPhi);

      // and check if adjacent towers consistent w/ a streak
      const bool isUpNotStreak = IsNeighborNotStreaky(m_ohMap.towers[iEta][iUp]);
      const bool isDownNotStreak = IsNeighborNotStreaky(m_ohMap.towers[iEta][iDown]);
      if (isUpNotStreak || isDownNotStreak)
      {
        continue;
      }

      // finally, increment no. of streaky towers for this phi
      // and this phi + 1
      ++m_ohNumStreak[iPhi];
      ++m_ohNumStreak[iUp];

      // fill histograms
      m_hists["nstreakperphi"]->Fill(iPhi);
      m_hists["nstreakperphi"]->Fill(iPhi);
      m_hists["nstreaktwretavsphi"]->Fill(iEta, iPhi);
      m_hists["nstreaktwretavsphi"]->Fill(iEta, iUp);

    }  // end eta loop
  }    // end phi loop

  // now find longest streak
  const uint32_t nMaxStreak = *std::max_element(m_ohNumStreak.begin(), m_ohNumStreak.end());
  m_hists["nmaxstreak"]->Fill(nMaxStreak);

  // return if streak length above threshold
  return (nMaxStreak > m_config.minNumTwrsInStreak);

}  // end 'ApplyFilter()'

// ----------------------------------------------------------------------------
//! Construct histograms
// ----------------------------------------------------------------------------
void StreakSidebandFilter::BuildHistograms(const std::string& module, const std::string& tag)
{
  // print debug message
  if (m_config.debug && (m_config.verbosity > 2))
  {
    std::cout << "StreakSidebandFilter::BuildHistograms(std::string) Constructing histograms" << std::endl;
  }

  // make sure module name is lower case
  std::string moduleAndFilterName = module + "_" + m_name;

  // names of variables to be histogramed
  const std::vector<std::string> varNames = {
      "nmaxstreak",
      "nstreakperphi",
      "nstreaktwretavsphi"};

  // make qa-compliant hist names
  std::vector<std::string> histNames = BeamBackgroundFilterAndQADefs::MakeQAHistNames(varNames, moduleAndFilterName, tag);

  // construct histograms
  //   - n.b. reminder that there are
  //       24 ohcal towers in eta
  //       64 ohcal towers in phi
  m_hists[varNames[0]] = new TH1D(histNames[0].data(), "", 25, -0.5, 24.5);
  m_hists[varNames[1]] = new TH1D(histNames[1].data(), "", 65, -0.5, 64.5);
  m_hists[varNames[2]] = new TH2D(histNames[2].data(), "", 25, -0.5, 24.5, 65, -0.5, 64.5);
  return;

}  // end 'BuildHistograms(std::string&, std::string&)'

// private methods ============================================================

// ----------------------------------------------------------------------------
//! Grab input nodes
// ----------------------------------------------------------------------------
void StreakSidebandFilter::GrabNodes(PHCompositeNode* topNode)
{
  // print debug message
  if (m_config.debug && (m_config.verbosity > 2))
  {
    std::cout << "StreakSidebandFilter::GrabNodes(PHCompositeNode*) Grabbing input nodes" << std::endl;
  }

  m_ohContainer = findNode::getClass<TowerInfoContainer>(topNode, m_config.inNodeName);
  return;

}  // end 'GrabNodes(PHCompositeNode*)'

// ----------------------------------------------------------------------------
//! Check if tower not consistent w/ being in a streak
// ----------------------------------------------------------------------------
bool StreakSidebandFilter::IsTowerNotStreaky(const BeamBackgroundFilterAndQADefs::Tower& tower)
{
  // print debug message
  if (m_config.debug && (m_config.verbosity > 2))
  {
    std::cout << "StreakSidebandFilter::IsTowerNotStreaky() Checking if tower not consistent w/ streak" << std::endl;
  }

  const bool isBadStatus = !(tower.isGood);
  const bool isBelowEneCut = (tower.energy < m_config.minStreakTwrEne);
  return (isBadStatus || isBelowEneCut);

}  // end 'IsTowerNotStreaky(Tower& tower)'

// ----------------------------------------------------------------------------
//! Check if a neighboring tower consistent w/ a streak
// ----------------------------------------------------------------------------
bool StreakSidebandFilter::IsNeighborNotStreaky(const BeamBackgroundFilterAndQADefs::Tower& tower)
{
  // print debug message
  if (m_config.debug && (m_config.verbosity > 2))
  {
    std::cout << "StreakSidebandFilter::IsNeighborNotStreaky() Checking if neighboring tower not consistent w/ streak" << std::endl;
  }

  const bool isBadStatus = !(tower.isGood);
  const bool isAboveEneCut = (tower.energy > m_config.maxAdjacentTwrEne);
  return (isBadStatus || isAboveEneCut);

}  // end 'IsNeighborNotStreaky(Tower& tower)'

// end ========================================================================
