/// ===========================================================================
/*! \file    NullFilter.cc
 *  \authors Derek Anderson
 *  \date    11.20.2024
 *
 *  Part of the BeamBackgroundFilterAndQA module, this
 *  filter does nothing, but provides a template for
 *  others to copy and fill in.
 */
/// ===========================================================================

// module components
#include "NullFilter.h"

#include "BeamBackgroundFilterAndQADefs.h"

// root includes
#include <TH1.h>

// c++ includes
#include <iostream>
#include <map>
#include <memory>
#include <vector>

// ctor/dtor ==================================================================

// ----------------------------------------------------------------------------
//! Default ctor
// ----------------------------------------------------------------------------
NullFilter::NullFilter(const std::string& name)
  : BaseBeamBackgroundFilter(name)
{
}  // end ctor()

// ----------------------------------------------------------------------------
//! ctor accepting config struct
// ----------------------------------------------------------------------------
NullFilter::NullFilter(const Config& config, const std::string& name)
  : BaseBeamBackgroundFilter(name)
  , m_config(config)
{
}  // end ctor(Config&)

// public methods =============================================================

// ----------------------------------------------------------------------------
// Apply filter to check for beam background or not
// ----------------------------------------------------------------------------
bool NullFilter::ApplyFilter(PHCompositeNode* topNode)
{
  // print debug message
  if (m_config.debug && (m_config.verbosity > 2))
  {
    std::cout << "NullFilter::ApplyFilter() Checking if streak found in OHCal via their sidebands" << std::endl;
  }

  // grab input node(s)
  GrabNodes(topNode);

  //... actual filter algorithm goes here ...//

  // other histograms filled similarly
  m_hists["test"]->Fill(1);

  // should return
  //   true,  if found beam background
  //   false, if no beam background found
  return false;

}  // end 'ApplyFilter()'

// ----------------------------------------------------------------------------
//! Construct histograms
// ----------------------------------------------------------------------------
void NullFilter::BuildHistograms(const std::string& module, const std::string& tag)
{
  // print debug message
  if (m_config.debug && (m_config.verbosity > 2))
  {
    std::cout << "NullFilter::BuildHistograms(std::string) Constructing histograms" << std::endl;
  }

  // make sure module name is lower case
  std::string moduleAndFilterName = module + "_" + m_name;

  // names of variables to be histogramed
  const std::vector<std::string> varNames = {
      "test"
      //... variable names like NStreakTwr should go here ...//
  };

  // make qa-compliant hist names
  std::vector<std::string> histNames = BeamBackgroundFilterAndQADefs::MakeQAHistNames(varNames, moduleAndFilterName, tag);

  // construct histograms
  m_hists[varNames[0]] = new TH1D(histNames[0].data(), "", 2, -0.5, 1.5);
  return;

}  // end 'BuildHistograms(std::string&, std::string&)'

// private methods ============================================================

// ----------------------------------------------------------------------------
//! Grab input nodes
// ----------------------------------------------------------------------------
void NullFilter::GrabNodes(PHCompositeNode* /*topNode*/)
{
  // print debug message
  if (m_config.debug && (m_config.verbosity > 2))
  {
    std::cout << "NullFilter::GrabNodes(PHCompositeNode*) Grabbing input nodes" << std::endl;
  }

  //... grab input nodes off the node tree here ...//
  return;

}  // end 'GrabNodes(PHCompositeNode*)'

// end ========================================================================
