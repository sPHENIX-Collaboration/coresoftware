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

#ifndef BEAMBACKGROUNDFILTERANDQA_H
#define BEAMBACKGROUNDFILTERANDQA_H

// module components
#include "NullFilter.h"
#include "StreakSidebandFilter.h"

// f4a libraries
#include <fun4all/SubsysReco.h>

// c++ utilities
#include <map>
#include <memory>
#include <string>
#include <vector>

// forward declarations
class BaseBeamBackgroundFilter;
class Fun4AllHistoManager;
class PHCompositeNode;
class recoConsts;
class TH1;

// ============================================================================
//! Filter beam background events and create QA
// ============================================================================
/*! A F4A module to filter out events with significant
 *  beam background and produce some relevant QA
 *  histograms.
 */
class BeamBackgroundFilterAndQA : public SubsysReco
{
 public:
  // ========================================================================
  //! User options for module
  // =======================================================================
  struct Config
  {
    // turn modes on/off
    bool debug = true;
    bool doQA = true;
    bool doEvtAbort = false;

    ///! module name
    std::string moduleName = "BeamBackgroundFilterAndQA";

    ///! histogram tags
    std::string histTag = "";

    ///! which filters to apply
    std::vector<std::string> filtersToApply = {"Null", "StreakSideband"};

    ///! filter configurations
    NullFilter::Config null;
    StreakSidebandFilter::Config sideband;
    //... add other configurations here ...//
  };

  // ctor/dtor
  BeamBackgroundFilterAndQA(const std::string& name = "BeamBackgroundFilterAndQA");
  BeamBackgroundFilterAndQA(const Config& config);
  ~BeamBackgroundFilterAndQA() override;

  // setters
  void SetConfig(const Config& config) { m_config = config; }

  // getters
  Config GetConfig() const { return m_config; }

  // f4a methods
  int Init(PHCompositeNode* /*topNode*/) override;
  int process_event(PHCompositeNode* topNode) override;
  int End(PHCompositeNode* /*topNode*/) override;

 private:
  // private methods
  void InitFilters();
  void InitFlags();
  void InitHistManager();
  void BuildHistograms();
  void RegisterHistograms();
  bool ApplyFilters(PHCompositeNode* topNode);

  ///! histogram manager
  Fun4AllHistoManager* m_manager{nullptr};

  ///! reco consts (for flags)
  recoConsts* m_consts{nullptr};

  ///! module-wide histograms
  std::map<std::string, TH1*> m_hists;

  ///! module configuration
  Config m_config;

  ///! filters
  std::map<std::string, std::unique_ptr<BaseBeamBackgroundFilter>> m_filters;

};  // end BeamBackgroundFilterAndQA

#endif

// end ========================================================================
