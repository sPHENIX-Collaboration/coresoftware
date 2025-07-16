/// ===========================================================================
/*! \file   TrksInJetQA.h
 *  \author Derek Anderson
 *  \date   03.25.2024
 *
 *  A "small" Fun4All module to produce QA plots for tracks,
 *  hits, and more.
 */
/// ============================================================================

#ifndef TRKSINJETQA_H
#define TRKSINJETQA_H

// module utilities
#include "JetQADefs.h"
#include "TrksInJetQAConfig.h"
#include "TrksInJetQAHist.h"
#include "TrksInJetQAInJetFiller.h"
#include "TrksInJetQAInclusiveFiller.h"

// calo trigger includes
#include <calotrigger/TriggerAnalyzer.h>

// f4a includes
#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

// phool includes
#include <phool/PHCompositeNode.h>
#include <phool/phool.h>

// qa utilities
#include <qautils/QAHistManagerDef.h>

// root includes
#include <TFile.h>

// c++ utilities
#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <optional>
#include <string>
#include <utility>
#include <vector>

// ============================================================================
//! Tracks-in-jet QA
// ============================================================================
/*! This Fun4All module generates histograms to QA low-level
 *  info of tracks found inside a jet, as well as some
 *  kinematic information.
 */ 
class TrksInJetQA : public SubsysReco
{
 public:

  ///! enumerates possibe output modes
  enum OutMode
  {
    File,
    QA
  };

  // ctor/dtor
  TrksInJetQA(const std::string& name = "TrksInJetQA");
  ~TrksInJetQA() override;

  // setters
  void SetOutFileName(const std::string& name) { m_outFileName = name; }
  void SetHistPrefix(const std::string& prefix) { m_histPrefix = prefix; }
  void SetHistSuffix(const std::string& suffix) { m_histSuffix = suffix; }
  void SetTrgToSelect(const uint32_t trig = JetQADefs::GL1::MBDNSJet1)
  {
    m_doTrgSelect = true;
    m_trgToSelect = trig;
  }

  // public methods
  void Configure(const TrksInJetQAConfig& config,
                 std::optional<TrksInJetQAHist> hist = std::nullopt);

  // f4a methods
  int Init(PHCompositeNode* /*topNode*/) override;
  int process_event(PHCompositeNode* topNode) override;
  int End(PHCompositeNode* /*topNode*/) override;

 private:

  // private methods
  void InitOutput();
  void InitHistograms();
  void RegisterHistograms();

  // io members
  TFile* m_outFile {nullptr};
  std::string m_moduleName;
  std::string m_outFileName {"tracksInJetsQA.root"};
  Fun4AllHistoManager* m_manager {nullptr};
  TriggerAnalyzer* m_analyzer {nullptr};

  // trigger selection
  bool m_doTrgSelect {false};
  uint32_t m_trgToSelect {JetQADefs::GL1::MBDNSJet1};

  // optional prefix, suffix for histograms
  std::optional<std::string> m_histPrefix {std::nullopt};
  std::optional<std::string> m_histSuffix {std::nullopt};

  // submodules to run
  std::unique_ptr<TrksInJetQAInJetFiller> m_inJet;
  std::unique_ptr<TrksInJetQAInclusiveFiller> m_inclusive;

  // module utilities
  TrksInJetQAConfig m_config;
  TrksInJetQAHist m_hist;

};  // end TrksInJetQA

#endif

// end ========================================================================
