/// ===========================================================================
/*! \file    BeamBackgroundFilterAndQA.h
 *  \authors Hanpu Jiang, Derek Anderson
 *  \date    10.21.2024
 *
 *  A F4A module to filter out events with significant
 *  beam background (the so-called "streaky events")
 *  and produce some relevant QA histograms.
 */
/// ===========================================================================

#define BEAMBACKGROUNDFILTERANDQA_CC

// module components
#include "BeamBackgroundFilterAndQA.h"
#include "BeamBackgroundFilterAndQADefs.h"

// calo base
#include <calobase/TowerInfoContainer.h>

// f4a libraries
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllHistoManager.h>
#include <ffaobjects/FlagSavev1.h>

// phool libraries
#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/PHCompositeNode.h>

// qa utilities
#include <qautils/QAHistManagerDef.h>

// root libraries
#include <TH1.h>

// c++ utiilites
#include <cassert>
#include <iostream>



// ctor/dtor ==================================================================

// ----------------------------------------------------------------------------
//! Default module constructor
// ----------------------------------------------------------------------------
BeamBackgroundFilterAndQA::BeamBackgroundFilterAndQA(const std::string& name, const bool debug)
  : SubsysReco(name)
  , m_manager(nullptr)
  , m_flags(name)
{

  // print debug message
  if (debug && (Verbosity() > 0))
  {
    std::cout << "BeamBackgroundFilterAndQA::BeamBackgroundFilterAndQA(const std::string &name) Calling ctor" << std::endl;
  }

}  // end ctor(std::string&, bool)'



// ----------------------------------------------------------------------------
//! Module constructor accepting a configuration
// ----------------------------------------------------------------------------
BeamBackgroundFilterAndQA::BeamBackgroundFilterAndQA(const Config& config)
  : SubsysReco(config.moduleName)
  , m_manager(nullptr)
  , m_flags(config.moduleName)
  , m_config(config)
{

  // print debug message
  if (m_config.debug && (Verbosity() > 0))
  {
    std::cout << "BeamBackgroundFilterAndQA::BeamBackgroundFilterAndQA(const std::string &name) Calling ctor" << std::endl;
  }

}  // end ctor(BeamBackgroundFilterAndQAConfig&)'



// ----------------------------------------------------------------------------
//! Module destructor
// ----------------------------------------------------------------------------
BeamBackgroundFilterAndQA::~BeamBackgroundFilterAndQA()
{

  // print debug message
  if (m_config.debug && (Verbosity() > 0))
  {
    std::cout << "BeamBackgroundFilterAndQA::~BeamBackgroundFilterAndQA() Calling dtor" << std::endl;
  }

  /* nothing to do */

}  // end dtor()



// fun4all methods ============================================================

// ----------------------------------------------------------------------------
//! Initialize module
// ----------------------------------------------------------------------------
int BeamBackgroundFilterAndQA::Init(PHCompositeNode* topNode)
{

  if (m_config.debug)
  {
    std::cout << "BeamBackgroundFilterAndQA::Init(PHCompositeNode *topNode) Initializing" << std::endl;
  }

  // initialize relevant filters, histograms
  InitFilters();
  InitFlags(topNode);
  BuildHistograms();

  // if needed, initialize histograms + manager
  if (m_config.doQA)
  {
    InitHistManager();
    RegisterHistograms();
  }
  return Fun4AllReturnCodes::EVENT_OK;

}  // end 'Init(PHCompositeNode*)'



// ----------------------------------------------------------------------------
//! Grab inputs, check for beam background, and fill histograms
// ----------------------------------------------------------------------------
int BeamBackgroundFilterAndQA::process_event(PHCompositeNode* topNode)
{

  if (m_config.debug)
  {
    std::cout << "BeamBackgroundFilterAndQA::process_event(PHCompositeNode *topNode) Processing event" << std::endl;
  }

  // reset flags
  SetDefaultFlags();

  // check for beam background
  const bool hasBeamBkgd = ApplyFilters(topNode);

  // update flags on node tree
  UpdateFlags(topNode);

  // if debugging, print out flags
  if (m_config.debug)
  {
    m_flags.printint();
  }

  // if it does, abort event
  if (hasBeamBkgd && m_config.doEvtAbort)
  {
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  else
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

}  // end 'process_event(PHCompositeNode*)'



// ----------------------------------------------------------------------------
//! Run final calculations
// ----------------------------------------------------------------------------
int BeamBackgroundFilterAndQA::End(PHCompositeNode* /*topNode*/)
{

  if (m_config.debug)
  {
    std::cout << "BeamBackgroundFilterAndQA::End(PHCompositeNode *topNode) This is the end..." << std::endl;
  }

  //... nothing to do ...//
  return Fun4AllReturnCodes::EVENT_OK;

}  // end 'End(PHCompositeNode*)'



// private methods ============================================================

// ----------------------------------------------------------------------------
//! Initialize filters
// ---------------------------------------------------------------------------
void BeamBackgroundFilterAndQA::InitFilters()
{

  // print debug message
  if (m_config.debug && (Verbosity() > 1))
  {
    std::cout << "BeamBackgroundFilterAndQA::InitFilters() Initializing background filters" << std::endl;
  }

  m_filters["Null"] = std::make_unique<NullFilter>(m_config.null, "Null");
  m_filters["StreakSideband"] = std::make_unique<StreakSidebandFilter>(m_config.sideband, "StreakSideband");
  //... other filters added here ...//
  return;

}  // end 'InitFilters()'



// ----------------------------------------------------------------------------
//! Initialize flags
// ----------------------------------------------------------------------------
void BeamBackgroundFilterAndQA::InitFlags(PHCompositeNode* topNode)
{

  // print debug message
  if (m_config.debug && (Verbosity() > 1))
  {
    std::cout << "BeamBackgroundFilterAndQA::InitFlags() Initializing flags" << std::endl;
  }

  // add node for flags
  PHNodeIterator itNode(topNode);
  PHCompositeNode* parNode = dynamic_cast<PHCompositeNode*>(itNode.findFirst("PHCompositeNode", "PAR"));
  if (!parNode)
  {
    std::cerr << PHWHERE << " WARNING: No PAR node found! Cannot add node for background flags to node tree!" << std::endl;
  }
  else
  {
    m_flags.SaveToNodeTree(parNode, m_config.flagPrefix);
  }

  // initialize flags
  SetDefaultFlags();
  return;

}  // end 'InitFlags()'



// ----------------------------------------------------------------------------
//! Initialize histogram manager
// ----------------------------------------------------------------------------
void BeamBackgroundFilterAndQA::InitHistManager()
{

  // print debug message
  if (m_config.debug && (Verbosity() > 0))
  {
    std::cout << "BeamBackgroundFilterAndQA::InitHistManager() Initializing histogram manager" << std::endl;
  }

  m_manager = QAHistManagerDef::getHistoManager();
  if (!m_manager)
  {
    std::cerr << PHWHERE << ": PANIC! Couldn't grab histogram manager!" << std::endl;
    assert(m_manager);
  }
  return;

}  // end 'InitHistManager()'



// ----------------------------------------------------------------------------
//! Build histograms
// ----------------------------------------------------------------------------
void BeamBackgroundFilterAndQA::BuildHistograms()
{

  // print debug message
  if (m_config.debug && (Verbosity() > 0))
  {
    std::cout << "BeamBackgroundFilterAndQA::BuildHistograms() Creating histograms" << std::endl;
  }

  // construct module-wide variable names
  std::vector<std::string> varNames = {"nevts_overall"};
  for (const std::string& filterToApply : m_config.filtersToApply)
  {
    varNames.push_back("nevts_" + filterToApply);
  }

  // get module-wide histogram names
  std::vector<std::string> histNames = BeamBackgroundFilterAndQADefs::MakeQAHistNames(varNames, m_config.moduleName, m_config.histTag);

  // create module-wide histograms
  for (std::size_t iVar = 0; iVar < varNames.size(); ++iVar)
  {
    m_hists[varNames[iVar]] = new TH1D(histNames[iVar].data(), "", 3, -0.5, 2.5);
    m_hists[varNames[iVar]]->GetXaxis()->SetBinLabel(1, "All");
    m_hists[varNames[iVar]]->GetXaxis()->SetBinLabel(2, "No beam bkgd.");
    m_hists[varNames[iVar]]->GetXaxis()->SetBinLabel(3, "Beam bkgd.");
  }

  // build filter-specific histograms
  for (const std::string& filterToApply : m_config.filtersToApply)
  {
    m_filters.at(filterToApply)->BuildHistograms(m_config.moduleName, m_config.histTag);
  }
  return;

}  // end 'BuildHistograms()'



// ----------------------------------------------------------------------------
//! Register histograms
// ----------------------------------------------------------------------------
void BeamBackgroundFilterAndQA::RegisterHistograms()
{

  // print debug message
  if (m_config.debug && (Verbosity() > 0))
  {
    std::cout << "BeamBackgroundFilterAndQA::RegisterHistograms() Registering histograms w/ manager" << std::endl;
  }

  // register module-wide histograms
  for (auto& hist : m_hists)
  {
    m_manager->registerHisto(hist.second);
  }

  // register filter-specific histograms
  for (const std::string& filterToApply : m_config.filtersToApply)
  {
    m_filters.at(filterToApply)->RegisterHistograms(m_manager);
  }
  return;

}  // end 'RegisterHistograms()'



// ----------------------------------------------------------------------------
//! Set default values of flags
// ----------------------------------------------------------------------------
void BeamBackgroundFilterAndQA::SetDefaultFlags()
{

  // print debug message
  if (m_config.debug && (Verbosity() > 1))
  {
    std::cout << "BeamBackgroundFilterAndQA::SetDefaultFlags() Setting defuault flag values" << std::endl;
  }

  // set default values
  for (const std::string& filterToApply : m_config.filtersToApply)
  {
    m_flags.set_int_param(MakeFlagName(filterToApply), 0);
  }
  m_flags.set_int_param(MakeFlagName(), 0);
  return;

}  // end 'SetDefaultFlags()'



// ----------------------------------------------------------------------------
//! Update flags on the node tree
// ----------------------------------------------------------------------------
void BeamBackgroundFilterAndQA::UpdateFlags(PHCompositeNode* topNode)
{

  // print debug message
  if (m_config.debug && (Verbosity() > 0))
  {
    std::cout << "BeamBackgroundFilterAndQA::UpdateFlags(PHCompositeNode*) Updating flags on the node tree" << std::endl;
  }

  PHNodeIterator itNode(topNode);
  PHCompositeNode* parNode = dynamic_cast<PHCompositeNode*>(itNode.findFirst("PHCompositeNode", "PAR"));
  if (parNode)
  {
    m_flags.UpdateNodeTree(parNode, m_config.flagPrefix);
  }
  return;

}  // end 'UpdateFlags(PHCompositeNode*)'



// ----------------------------------------------------------------------------
//! Apply relevant filters
// ----------------------------------------------------------------------------
bool BeamBackgroundFilterAndQA::ApplyFilters(PHCompositeNode* topNode)
{

  // print debug message
  if (m_config.debug && (Verbosity() > 0))
  {
    std::cout << "BeamBackgroundFilterAndQA::ApplyFilters(PHCompositeNode*) Creating histograms" << std::endl;
  }

  // apply individual filters
  bool hasBkgd = false;
  for (const std::string& filterToApply : m_config.filtersToApply)
  {
    const bool filterFoundBkgd = m_filters.at(filterToApply)->ApplyFilter(topNode);
    if (filterFoundBkgd)
    {
      m_hists["nevts_" + filterToApply]->Fill(BeamBackgroundFilterAndQADefs::Status::HasBkgd);
      m_flags.set_int_param(MakeFlagName(filterToApply), 1);
    }
    else
    {
      m_hists["nevts_" + filterToApply]->Fill(BeamBackgroundFilterAndQADefs::Status::NoBkgd);
    }
    m_hists["nevts_" + filterToApply]->Fill(BeamBackgroundFilterAndQADefs::Status::Evt);
    hasBkgd += filterFoundBkgd;
  }

  // fill overall histograms and return
  m_hists["nevts_overall"]->Fill(BeamBackgroundFilterAndQADefs::Status::Evt);
  if (hasBkgd)
  {
    m_hists["nevts_overall"]->Fill(BeamBackgroundFilterAndQADefs::Status::HasBkgd);
    m_flags.set_int_param(MakeFlagName(), 1);
  }
  else
  {
    m_hists["nevts_overall"]->Fill(BeamBackgroundFilterAndQADefs::Status::NoBkgd);
  }
  return hasBkgd;

}  // end 'ApplyFilters(PHCompositeNode*)'



// ----------------------------------------------------------------------------
//! Create flag name
// ----------------------------------------------------------------------------
/*! Helper method to create a flag (PHParameters) name from a
 *  provided filter name.  If no filter is provided, method
 *  will return just the flag prefix.
 */
std::string BeamBackgroundFilterAndQA::MakeFlagName(const std::string& filter)
{

  // print debug message
  if (m_config.debug && (Verbosity() > 2))
  {
    std::cout << "BeamBackgroundFilterAndQA::MakeFlagName(std::string&) Creating flag name" << std::endl;
  }

  // by default, return flag prefix; otherwise,
  // combine prefix and filter name
  if (filter.empty())
  {
    return m_config.flagPrefix;
  }
  else
  {
    return m_config.flagPrefix + "_" + filter + "Filter";
  }

}  // end 'MakeFlagName(std::string&)'

// end ========================================================================
