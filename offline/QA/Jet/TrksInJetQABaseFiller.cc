// ----------------------------------------------------------------------------
// 'TrksInJetQABaseFiller.cc'
// Derek Anderson
// 04.11.2024
//
// A submodule for the TrksInJetQA F4A module to produce
// QA histograms for tracks and more in jets
// ----------------------------------------------------------------------------

#define TRKSINJETQABASEFILLER_CC

#include "TrksInJetQABaseFiller.h"



// ctor/dtor ------------------------------------------------------------------

TrksInJetQABaseFiller::TrksInJetQABaseFiller(
  TrksInJetQAConfig& config,
  TrksInJetQAHist& hist
) {

  // grab utilities
  m_config = config;
  m_hist   = hist;

  // initialize managers
  if (m_config.doHitQA)   m_hitManager   = std::make_unique<TrksInJetQAHitManager>(m_config, m_hist);
  if (m_config.doClustQA) m_clustManager = std::make_unique<TrksInJetQAClustManager>(m_config, m_hist);
  if (m_config.doTrackQA) m_trackManager = std::make_unique<TrksInJetQATrkManager>(m_config, m_hist);
  if (m_config.doJetQA)   m_jetManager   = std::make_unique<TrksInJetQAJetManager>(m_config, m_hist);

}  // end ctor()'



TrksInJetQABaseFiller::~TrksInJetQABaseFiller() {

  /* nothing to do */

}  // end dtor



// public methods -------------------------------------------------------------

void TrksInJetQABaseFiller::MakeHistograms(std::string label) {

  // initialize relevant submodules
  if (m_config.doHitQA)   m_hitManager   -> MakeHistograms(label);
  if (m_config.doClustQA) m_clustManager -> MakeHistograms(label);
  if (m_config.doTrackQA) m_trackManager -> MakeHistograms(label);
  if (m_config.doJetQA)   m_jetManager   -> MakeHistograms(label);
  return;

}  // end 'MakeHistograms(std::string)'



void TrksInJetQABaseFiller::SaveHistograms(TFile* outFile, std::string outDirName) {

  TDirectory* outDir = outFile -> mkdir(outDirName.data());
  if (!outDir) {
    std::cerr << PHWHERE << ": PANIC: unable to make output directory!" << std::endl;
    assert(outDir);
  }

  if (m_config.doHitQA)   m_hitManager   -> SaveHistograms(outDir, m_config.hitOutDir);
  if (m_config.doClustQA) m_clustManager -> SaveHistograms(outDir, m_config.clustOutDir);
  if (m_config.doTrackQA) m_trackManager -> SaveHistograms(outDir, m_config.trackOutDir);
  if (m_config.doJetQA)   m_jetManager   -> SaveHistograms(outDir, m_config.jetOutDir);
  return;

}  // end 'SaveHistograms(TFile*, std::string)'



void TrksInJetQABaseFiller::GrabHistograms(
  std::vector<TH1D*>& vecOutHist1D,
  std::vector<TH2D*>& vecOutHist2D
) {

  if (m_config.doHitQA)   m_hitManager   -> GrabHistograms(vecOutHist1D, vecOutHist2D);
  if (m_config.doClustQA) m_clustManager -> GrabHistograms(vecOutHist1D, vecOutHist2D);
  if (m_config.doTrackQA) m_trackManager -> GrabHistograms(vecOutHist1D, vecOutHist2D);
  if (m_config.doJetQA)   m_jetManager   -> GrabHistograms(vecOutHist1D, vecOutHist2D);
  return;

}  // end 'GrabHistograms(std::vector<TH1D*>&, std::vector<TH2D*>&)'



// private methods ------------------------------------------------------------

void TrksInJetQABaseFiller::GetNodes(PHCompositeNode* topNode) {

  // grab necessary jet nodes
  if (m_config.doJetQA) {
    m_jetMap = findNode::getClass<JetContainer>(topNode, m_config.jetInNode.data());
    if (!m_jetMap) {
      std::cerr << PHWHERE << ": PANIC: couldn't grab jet map from node tree!" << std::endl;
      assert(m_jetMap);
    }
  }

  // grab necessary track nodes
  if (m_config.doTrackQA) {
    m_trkMap = findNode::getClass<SvtxTrackMap>(topNode, m_config.trkInNode.data());
    if (!m_trkMap) {
      std::cerr << PHWHERE << ": PANIC: couldn't grab track map from node tree!" << std::endl;
      assert(m_trkMap);
    }
  }

  // grab necessary cluster nodes
  if (m_config.doClustQA) {

    m_actsGeom = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
    if (!m_actsGeom) {
      std::cerr << PHWHERE << ": PANIC: couldn't grab ACTS geometry from node tree!" << std::endl;
      assert(m_actsGeom);
    }

    m_clustMap = findNode::getClass<TrkrClusterContainer>(topNode, m_config.clustInNode.data());
    if (!m_clustMap) {
      std::cerr << PHWHERE << ": PANIC: couldn't grab cluster map from node tree!" << std::endl;
      assert(m_clustMap);
    }
  }

  // grab necessary hit nodes
  if (m_config.doHitQA) {
    m_hitMap = findNode::getClass<TrkrHitSetContainer>(topNode, m_config.hitInNode.data());
    if (!m_hitMap) {
      std::cerr << PHWHERE << ": PANIC: couldn't grab hit map from node tree!" << std::endl;
      assert(m_hitMap);
    }
  }
  return;

}  // end 'GetNodes(PHCompositeNode*)'

// end ------------------------------------------------------------------------
