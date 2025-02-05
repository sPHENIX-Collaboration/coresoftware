// ----------------------------------------------------------------------------
// 'TrksInJetQAJetManager.h'
// Derek Anderson
// 03.26.2024
//
// A submodule for the TrksInJetQA module
// to generate QA plots for jets
// ----------------------------------------------------------------------------

#ifndef TRACKSINJETSQAMAKER_JETQAHISTMANAGER_H
#define TRACKSINJETSQAMAKER_JETQAHISTMANAGER_H

// submodule definitions
#include "TrksInJetQABaseManager.h"

// jet includes
#include <jetbase/Jet.h>

// tracking includes
#include <trackbase_historic/SvtxTrack.h>

// root includes
#include <TH1.h>
#include <TH2.h>

// c++ utilities
#include <limits>
#include <optional>
#include <utility>
#include <vector>

// TrksInJetQAJetManager definition -------------------------------------------

class TrksInJetQAJetManager : public TrksInJetQABaseManager
{
 public:
  // histogram accessors
  enum Type
  {
    All
  };
  enum H1D
  {
    Eta,
    Phi,
    Pt,
    NTrk,
    PtSum
  };
  enum H2D
  {
    PtVsEta,
    PtSumVsPt
  };

  // histogram content
  struct JetQAContent
  {
    double eta = std::numeric_limits<double>::max();
    double phi = std::numeric_limits<double>::max();
    double pt = std::numeric_limits<double>::max();
    double nTrk = std::numeric_limits<double>::max();
    double ptSum = std::numeric_limits<double>::max();
  };

  // ctor/dtor
  using TrksInJetQABaseManager::TrksInJetQABaseManager;
  virtual ~TrksInJetQAJetManager() override = default;

  // public methods
  void GetInfo(Jet* jet, std::optional<std::vector<SvtxTrack*>> tracks = std::nullopt);

 private:
  // private methods
  void FillHistograms(const int type, JetQAContent& content);

  // inherited interal methods
  void DefineHistograms() override;

};  // end TrksInJetQAJetManager

#endif

// end ------------------------------------------------------------------------
