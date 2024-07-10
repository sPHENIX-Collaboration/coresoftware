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

// tracking includes
#include <trackbase/ActsGeometry.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase_historic/SvtxTrackMap.h>
// jet includes
#include <jetbase/JetContainer.h>

// phool libraries
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>

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
  //   - FIXME these should be smart pointers!
  ActsGeometry* m_actsGeom = NULL;
  TrkrHitSetContainer* m_hitMap = NULL;
  TrkrClusterContainer* m_clustMap = NULL;
  SvtxTrackMap* m_trkMap = NULL;
  JetContainer* m_jetMap = NULL;

  // submodules to use
  std::unique_ptr<TrksInJetQAHitManager> m_hitManager = NULL;
  std::unique_ptr<TrksInJetQAClustManager> m_clustManager = NULL;
  std::unique_ptr<TrksInJetQATrkManager> m_trackManager = NULL;
  std::unique_ptr<TrksInJetQAJetManager> m_jetManager = NULL;

  // module utilities
  TrksInJetQAConfig m_config;
  TrksInJetQAHist m_hist;

};  // end TrksInJetQABaseFiller

#endif

// end ------------------------------------------------------------------------
