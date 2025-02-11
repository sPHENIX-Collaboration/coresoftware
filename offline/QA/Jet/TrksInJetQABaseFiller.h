// ----------------------------------------------------------------------------
// 'TrksInJetQABaseFiller.h'
// Derek Anderson
// 04.11.2024
//
// A submodule for the TrksInJetQA F4A module to produce
// QA histograms for tracks and more in jets
// ----------------------------------------------------------------------------

#ifndef TRKSINJETQABASEFILLER_H
#define TRKSINJETQABASEFILLER_H

// submodule definitions
#include "TrksInJetQAClustManager.h"
#include "TrksInJetQAHitManager.h"
#include "TrksInJetQAJetManager.h"
#include "TrksInJetQATrkManager.h"

// module utilties
#include "TrksInJetQAConfig.h"
#include "TrksInJetQAHist.h"

// jet includes
#include <jetbase/JetContainer.h>

// phool libraries
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>

// tracking includes
#include <trackbase/ActsGeometry.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase_historic/SvtxTrackMap.h>

// root includes
#include <TFile.h>

// c++ includes
#include <string>

// TrksInJetQABaseFiller ------------------------------------------------------

class TrksInJetQABaseFiller
{
 public:
  // ctor/dtor
  TrksInJetQABaseFiller(const TrksInJetQAConfig& config, TrksInJetQAHist& hist);
  virtual ~TrksInJetQABaseFiller() = default;

  // public methods
  void MakeHistograms(const std::string& prefix = "", const std::string& suffix = "");
  void SaveHistograms(TFile* outFile, const std::string& outDirName);
  void GrabHistograms(std::vector<TH1D*>& vecOutHist1D, std::vector<TH2D*>& vecOutHist2D);

  // virtual public methods
  virtual void Fill(PHCompositeNode* topNode) = 0;

 protected:
  // private methods
  void GetNodes(PHCompositeNode* topNode);

  // necessary dst nodes
  ActsGeometry* m_actsGeom {nullptr};
  TrkrHitSetContainer* m_hitMap {nullptr};
  TrkrClusterContainer* m_clustMap {nullptr};
  SvtxTrackMap* m_trkMap {nullptr};
  JetContainer* m_jetMap {nullptr};

  // submodules to use
  std::unique_ptr<TrksInJetQAHitManager> m_hitManager {nullptr};
  std::unique_ptr<TrksInJetQAClustManager> m_clustManager {nullptr};
  std::unique_ptr<TrksInJetQATrkManager> m_trackManager {nullptr};
  std::unique_ptr<TrksInJetQAJetManager> m_jetManager {nullptr};

  // module utilities
  TrksInJetQAConfig m_config;
  TrksInJetQAHist m_hist;

};  // end TrksInJetQABaseFiller

#endif

// end ------------------------------------------------------------------------
