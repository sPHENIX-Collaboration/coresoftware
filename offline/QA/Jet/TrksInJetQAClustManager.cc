// ----------------------------------------------------------------------------
// 'TrksInJetQAClustManager.cc'
// Derek Anderson
// 03.25.2024
//
// A submodule for the TrksInJetQA module to generate
// QA plots for track clusters
// ----------------------------------------------------------------------------

#define TRKSINJETQACLUSTMANAGER_CC

// submodule definition
#include "TrksInJetQAClustManager.h"



// public methods -------------------------------------------------------------

void TrksInJetQAClustManager::GetInfo(TrkrCluster* cluster, TrkrDefs::cluskey& clustKey, ActsGeometry* actsGeom) {

  // check which subsystem cluster is in
  const uint16_t layer  = TrkrDefs::getLayer(clustKey);
  const bool     isMvtx = IsInMvtx(layer);
  const bool     isIntt = IsInIntt(layer);
  const bool     isTpc  = IsInTpc(layer);

  // get cluster position
  Acts::Vector3 actsPos = actsGeom -> getGlobalPosition(clustKey, cluster);

  // collect cluster info
  ClustQAContent content {
    .x = actsPos(0),
    .y = actsPos(1),
    .z = actsPos(2),
    .r = std::hypot( actsPos(0), actsPos(1) )
  };

  // fill histograms
  FillHistograms(Type::All, content);
  if (isMvtx) {
    FillHistograms(Type::Mvtx, content);
  } else if (isIntt) {
    FillHistograms(Type::Intt, content);
  } else if (isTpc) {
    FillHistograms(Type::Tpc, content);
  }

}  // end GetInfo(TrkrCluster*, TrkrDefs::cluskey&, ActsGeometry*)'



// private methods ------------------------------------------------------------

void TrksInJetQAClustManager::FillHistograms(const int type, ClustQAContent& content) {

  // fill 1d histograms
  m_vecHist1D.at(type).at(H1D::PosX) -> Fill(content.x);
  m_vecHist1D.at(type).at(H1D::PosY) -> Fill(content.y);
  m_vecHist1D.at(type).at(H1D::PosZ) -> Fill(content.z);
  m_vecHist1D.at(type).at(H1D::PosR) -> Fill(content.r);

  // fill 2d histograms
  m_vecHist2D.at(type).at(H2D::PosYvsX) -> Fill(content.x, content.y);
  m_vecHist2D.at(type).at(H2D::PosRvsZ) -> Fill(content.z, content.r);
  return;

}  //  end 'FillHistograms(int, ClustQAContent&)'



void TrksInJetQAClustManager::DefineHistograms() {

  // grab binning schemes
  std::vector<BinDef> vecBins = m_hist.GetVecHistBins();

  // set histogram types
  m_vecHistTypes.push_back( "Mvtx" );
  m_vecHistTypes.push_back( "Intt" );
  m_vecHistTypes.push_back( "Tpc"  );
  m_vecHistTypes.push_back( "All"  );

  // define 1d histograms
  m_vecHistDef1D.push_back( std::make_tuple( "ClustPosX", vecBins.at(TrksInJetQAHist::Var::PosXY) ));
  m_vecHistDef1D.push_back( std::make_tuple( "ClustPosY", vecBins.at(TrksInJetQAHist::Var::PosXY) ));
  m_vecHistDef1D.push_back( std::make_tuple( "ClustPosZ", vecBins.at(TrksInJetQAHist::Var::PosZ)  ));
  m_vecHistDef1D.push_back( std::make_tuple( "ClustPosR", vecBins.at(TrksInJetQAHist::Var::PosR)  ));

  // define 2d histograms
  m_vecHistDef2D.push_back(
    std::make_tuple(
      "ClustPosYvsX",
      vecBins.at(TrksInJetQAHist::Var::PosXY),
      vecBins.at(TrksInJetQAHist::Var::PosXY)
    )
  );
  m_vecHistDef2D.push_back(
    std::make_tuple(
      "ClustPosRvsZ",
      vecBins.at(TrksInJetQAHist::Var::PosZ),
      vecBins.at(TrksInJetQAHist::Var::PosR)
    )
  );
  return;

}  // end 'BuildHistograms()'

// end ------------------------------------------------------------------------
