// ----------------------------------------------------------------------------
// 'TrksInJetQA.h'
// Derek Anderson
// 03.25.2024
//
// A "small" Fun4All module to produce QA plots for tracks,
// hits, and more.
// ----------------------------------------------------------------------------

#ifndef TRKSINJETQA_H
#define TRKSINJETQA_H

// c++ utilities
#include <string>
#include <vector>
#include <cassert>
#include <utility>
#include <optional>
// root libraries
#include <TFile.h>
// f4a libraries
#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllHistoManager.h>
// phool libraries
#include <phool/phool.h>
#include <phool/PHCompositeNode.h>
// qa utilities
#include <qautils/QAHistManagerDef.h>
// module utilities
#include "TrksInJetQAHist.h"
#include "TrksInJetQAConfig.h"
// submodule definitions
#include "TrksInJetQAInJetFiller.h"
#include "TrksInJetQAInclusiveFiller.h"



// TrksInJetQA definition -----------------------------------------------------

class TrksInJetQA : public SubsysReco {

  public:

    //  output modes
    enum OutMode {File, QA};

    // ctor/dtor
    TrksInJetQA(const std::string& name);
    ~TrksInJetQA() override;

    // setters
    void SetOutFileName(const std::string& name)  {m_outFileName = name;}
    void SetHistSuffix(const std::string& suffix) {m_histSuffix  = suffix;}

    // public methods
    void Configure(
      TrksInJetQAConfig config,
      std::optional<TrksInJetQAHist> hist = std::nullopt
    );


    // f4a methods
    int Init(PHCompositeNode* topNode)          override;
    int process_event(PHCompositeNode* topNode) override;
    int End(PHCompositeNode* topNode)           override;

  private:

    // private methods
    void InitOutput();
    void InitHistograms();
    void RegisterHistograms();

    // io members
    //   - FIXME raw pointers should be smart ones!
    TFile*               m_outFile     = NULL;
    std::string          m_outFileName = "tracksInJetsQA.root";
    Fun4AllHistoManager* m_manager     = NULL;

    // optional suffix for histograms
    std::optional<std::string> m_histSuffix  = std::nullopt;

    // submodules to run
    std::unique_ptr<TrksInJetQAInJetFiller>     m_inJet;
    std::unique_ptr<TrksInJetQAInclusiveFiller> m_inclusive;

    // module utilities
    TrksInJetQAConfig m_config;
    TrksInJetQAHist   m_hist;

};  // end TrksInJetQA

#endif

// end ------------------------------------------------------------------------
