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

// module components
#include "BeamBackgroundFilterAndQA.h"
#include "BeamBackgroundFilterAndQADefs.h"

// f4a includes
#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>

// phool includes
#include <phool/phool.h>
#include <phool/recoConsts.h>

// qa utilities
#include <qautils/QAHistManagerDef.h>

// root includes
#include <TH1.h>

// c++ includes
#include <cassert>
#include <cstddef>
#include <iostream>
#include <utility>

// ctor/dtor ==================================================================

// ----------------------------------------------------------------------------
//! Default module constructor
// ----------------------------------------------------------------------------
BeamBackgroundFilterAndQA::BeamBackgroundFilterAndQA(const std::string& name)
  : SubsysReco(name)
{
}  // end ctor(std::string&)'

// ----------------------------------------------------------------------------
//! Module constructor accepting a configuration
// ----------------------------------------------------------------------------
BeamBackgroundFilterAndQA::BeamBackgroundFilterAndQA(const Config& config)
  : SubsysReco(config.moduleName)
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
int BeamBackgroundFilterAndQA::Init(PHCompositeNode* /*topNode*/)
{
  if (m_config.debug)
  {
    std::cout << "BeamBackgroundFilterAndQA::Init(PHCompositeNode *topNode) Initializing" << std::endl;
  }

  // initialize relevant filters, histograms
  InitFilters();
  InitFlags();
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

  // check for beam background
  const bool hasBeamBkgd = ApplyFilters(topNode);

  // if debugging, print out flags
  if (m_config.debug)
  {
    m_consts->PrintIntFlags();
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
void BeamBackgroundFilterAndQA::InitFlags()
{
  // print debug message
  if (m_config.debug && (Verbosity() > 1))
  {
    std::cout << "BeamBackgroundFilterAndQA::InitFlags() Initializing reco consts and flags" << std::endl;
  }

  m_consts = recoConsts::instance();
  for (const std::string& filterToApply : m_config.filtersToApply)
  {
    m_consts->set_IntFlag("HasBeamBackground_" + filterToApply + "Filter", 0);
  }
  m_consts->set_IntFlag("HasBeamBackground", 0);
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
      m_consts->set_IntFlag("HasBeamBackground_" + filterToApply + "Filter", 1);
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
    m_consts->set_IntFlag("HasBeamBackground", 1);
  }
  else
  {
    m_hists["nevts_overall"]->Fill(BeamBackgroundFilterAndQADefs::Status::NoBkgd);
  }
  return hasBkgd;

}  // end 'ApplyFilters(PHCompositeNode*)'

// end ========================================================================
