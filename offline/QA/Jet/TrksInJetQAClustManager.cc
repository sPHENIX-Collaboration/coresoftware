/// ===========================================================================
/*! \file   TrksInJetQAClustManager.cc
 *  \author Derek Anderson
 *  \date   03.25.2024
 *
 *  A submodule for the TrksInJetQA module to generate
 *  QA plots for track clusters
 */
/// ============================================================================

#define TRKSINJETQACLUSTMANAGER_CC

// submodule definition
#include "TrksInJetQAClustManager.h"

// public methods =============================================================

// ----------------------------------------------------------------------------
//! Grab information from a tracker cluster
// ----------------------------------------------------------------------------
void TrksInJetQAClustManager::GetInfo(TrkrCluster* cluster,
                                      TrkrDefs::cluskey& clustKey,
                                      ActsGeometry* actsGeom)
{
  // check which subsystem cluster is in
  const uint16_t layer = TrkrDefs::getLayer(clustKey);
  const bool isMvtx = IsInMvtx(layer);
  const bool isIntt = IsInIntt(layer);
  const bool isTpc = IsInTpc(layer);

  // get cluster position
  Acts::Vector3 actsPos = actsGeom->getGlobalPosition(clustKey, cluster);

  // collect cluster info
  ClustQAContent content{
      .x = actsPos(0),
      .y = actsPos(1),
      .z = actsPos(2),
      .r = std::hypot(actsPos(0), actsPos(1))};

  // fill histograms
  FillHistograms(Type::All, content);
  if (m_config.doSubsysHist)
  {
    if (isMvtx)
    {
      FillHistograms(Type::Mvtx, content);
    }
    else if (isIntt)
    {
      FillHistograms(Type::Intt, content);
    }
    else if (isTpc)
    {
      FillHistograms(Type::Tpc, content);
    }
  }
}  // end GetInfo(TrkrCluster*, TrkrDefs::cluskey&, ActsGeometry*)'

// private methods ============================================================

// ----------------------------------------------------------------------------
//! Fill tracker cluster histograms
// ----------------------------------------------------------------------------
void TrksInJetQAClustManager::FillHistograms(const int type, ClustQAContent& content)
{
  // fill 1d histograms
  m_mapHist1D[Index(type, H1D::PosX)]->Fill(content.x);
  m_mapHist1D[Index(type, H1D::PosY)]->Fill(content.y);
  m_mapHist1D[Index(type, H1D::PosZ)]->Fill(content.z);
  m_mapHist1D[Index(type, H1D::PosR)]->Fill(content.r);

  // fill 2d histograms
  m_mapHist2D[Index(type, H2D::PosYvsX)]->Fill(content.x, content.y);
  m_mapHist2D[Index(type, H2D::PosRvsZ)]->Fill(content.z, content.r);
}  //  end 'FillHistograms(int, ClustQAContent&)'

// ----------------------------------------------------------------------------
//! Define tracker cluster histograms
// ----------------------------------------------------------------------------
void TrksInJetQAClustManager::DefineHistograms()
{
  // grab binning schemes
  std::vector<TrksInJetQADefs::BinDef> vecBins = m_hist.GetVecHistBins();

  // set histogram types
  m_mapHistTypes[Type::All] = "All";
  if (m_config.doSubsysHist)
  {
    m_mapHistTypes[Type::Mvtx] = "Mvtx";
    m_mapHistTypes[Type::Intt] = "Intt";
    m_mapHistTypes[Type::Tpc] = "Tpc";
  }

  // define 1d histograms
  m_mapHistDef1D[H1D::PosX] = std::tuple("ClustPosX", vecBins.at(TrksInJetQAHist::Var::PosXY));
  m_mapHistDef1D[H1D::PosY] = std::tuple("ClustPosY", vecBins.at(TrksInJetQAHist::Var::PosXY));
  m_mapHistDef1D[H1D::PosZ] = std::tuple("ClustPosZ", vecBins.at(TrksInJetQAHist::Var::PosZ));
  m_mapHistDef1D[H1D::PosR] = std::tuple("ClustPosR", vecBins.at(TrksInJetQAHist::Var::PosR));

  // define 2d histograms
  m_mapHistDef2D[H2D::PosYvsX] = std::tuple("ClustPosYvsX",
                                            vecBins.at(TrksInJetQAHist::Var::PosXY),
                                            vecBins.at(TrksInJetQAHist::Var::PosXY));
  m_mapHistDef2D[H2D::PosRvsZ] = std::tuple("ClustPosRvsZ",
                                            vecBins.at(TrksInJetQAHist::Var::PosZ),
                                            vecBins.at(TrksInJetQAHist::Var::PosR));
}  // end 'BuildHistograms()'

// end ========================================================================
