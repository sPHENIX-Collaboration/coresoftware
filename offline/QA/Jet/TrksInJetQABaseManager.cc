// ----------------------------------------------------------------------------
// 'TrksInJetQABaseManager.cc'
// Derek Anderson
// 04.03.2024
//
// Base hist manager submodule for the TrksInJetQA module which
// consolidates methods/data common to all of the hist managers
// ----------------------------------------------------------------------------

#define TRKSINJETQABASEMANAGER_CC

// submodule definition
#include "TrksInJetQABaseManager.h"



// ctor/dtor ------------------------------------------------------------------

TrksInJetQABaseManager::TrksInJetQABaseManager(
  TrksInJetQAConfig& config,
  TrksInJetQAHist& hist
) {

  ResetVectors();

  // grab vectors
  m_config = config;
  m_hist   = hist;

}  // end ctor(TrksInJetQAConfig&, TrksInJetQAHist&)



TrksInJetQABaseManager::~TrksInJetQABaseManager() {

  ResetVectors();

}  // end dtor



// public methods -------------------------------------------------------------

void TrksInJetQABaseManager::MakeHistograms(std::string label) {

  DefineHistograms();
  BuildHistograms(label);
  return;

}  // end 'MakeHistograms(std::string)'



void TrksInJetQABaseManager::SaveHistograms(TDirectory* topDir, std::string outDirName) {

  TDirectory* outDir = topDir -> mkdir(outDirName.data());
  if (!outDir) {
    std::cerr << PHWHERE << ": PANIC: unable to make output directory!" << std::endl;
    assert(outDir);
  }

  outDir -> cd();
  for (auto hist1Ds : m_vecHist1D) {
    for (TH1D* hist1D : hist1Ds) {
      hist1D -> Write();
    }
  }
  for (auto hist2Ds : m_vecHist2D) {
    for (TH2D* hist2D : hist2Ds) {
      hist2D -> Write();
    }
  }
  return;

}  // end 'SaveHistograms(TDirectory*, std::string)'



void TrksInJetQABaseManager::GrabHistograms(
  std::vector<TH1D*>& vecOutHist1D,
  std::vector<TH2D*>& vecOutHist2D
) {

  for (auto hist1Ds : m_vecHist1D) {
    for (TH1D* hist1D : hist1Ds) {
      vecOutHist1D.push_back(hist1D);
    }
  }
  for (auto hist2Ds : m_vecHist2D) {
    for (TH2D* hist2D : hist2Ds) {
      vecOutHist2D.push_back(hist2D);
    }
  }
  return;

}  // end 'GrabHistograms(std::vector<TH1D*>&, std::vector<TH2D*>&)'



// private methods ------------------------------------------------------------

void TrksInJetQABaseManager::BuildHistograms(std::string label) {

  // build 1d histograms
  m_vecHist1D.resize( m_vecHistTypes.size() );
  for (size_t iType = 0; iType < m_vecHistTypes.size(); iType++) {
    for (HistDef1D histDef1D : m_vecHistDef1D) {

      // make name
      std::string sHistName("h");
      sHistName += m_vecHistTypes.at(iType);
      sHistName += std::get<0>(histDef1D);
      sHistName += "_";
      sHistName += label;

      // create histogram
      m_vecHist1D.at(iType).push_back(
        new TH1D(
          sHistName.data(),
          "",
          std::get<1>(histDef1D).first,
          std::get<1>(histDef1D).second.first,
          std::get<1>(histDef1D).second.second
        )
      );
    }  // end hist loop
  }  // end type loop

  // build 2d histograms
  m_vecHist2D.resize( m_vecHistTypes.size() );
  for (size_t iType = 0; iType < m_vecHistTypes.size(); iType++) {
    for (HistDef2D histDef2D : m_vecHistDef2D) {

      // make name
      std::string sHistName("h");
      sHistName += m_vecHistTypes.at(iType);
      sHistName += std::get<0>(histDef2D);
      sHistName += "_";
      sHistName += label;

      // create histogram
      m_vecHist2D.at(iType).push_back(
        new TH2D(
          sHistName.data(),
          "",
          std::get<1>(histDef2D).first,
          std::get<1>(histDef2D).second.first,
          std::get<1>(histDef2D).second.second,
          std::get<2>(histDef2D).first,
          std::get<2>(histDef2D).second.first,
          std::get<2>(histDef2D).second.second
        )
      );
    }  // end hist loop
  }  // end type loop
  return;

}  // end 'BuildHistograms(std::string)'



void TrksInJetQABaseManager::ResetVectors() {

  m_vecHist1D.clear();
  m_vecHist2D.clear();
  m_vecHistTypes.clear();
  m_vecHistDef1D.clear();
  m_vecHistDef2D.clear();
  return;

}  // end 'ResetVectors()'



// private helper methods -----------------------------------------------------

bool TrksInJetQABaseManager::IsInMvtx(const uint16_t layer) {

  return (layer < m_config.nMvtxLayer);

}  // end 'IsInMvtx(uint16_t)'



bool TrksInJetQABaseManager::IsInIntt(const uint16_t layer) {

  return (
    (layer >= m_config.nMvtxLayer) &&
    (layer <  m_config.nInttLayer + m_config.nMvtxLayer)
  );

}  // end 'IsInIntt(uint16_t)'



bool TrksInJetQABaseManager::IsInTpc(const uint16_t layer) {

  return (layer >= m_config.nMvtxLayer + m_config.nInttLayer);

}  // end 'IsInTpc(uint16_t)'

// end ------------------------------------------------------------------------
