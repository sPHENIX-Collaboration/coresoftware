/// ===========================================================================
/*! \file   TrksInJetQAHitManager.cc
 *  \author Derek Anderson
 *  \date   03.25.2024
 *
 *  A submodule for the TrksInJetQA module
 *  to generate QA plots for track hits
 */
/// ===========================================================================

#define TRKSINJETQAHITMANAGER_CC

// submodule definition
#include "TrksInJetQAHitManager.h"

// public methods =============================================================

// ----------------------------------------------------------------------------
//! Get information from a tracker hit
// ----------------------------------------------------------------------------
void TrksInJetQAHitManager::GetInfo(TrkrHit* hit,
                                    TrkrDefs::hitsetkey& setKey,
                                    TrkrDefs::hitkey& hitKey)
{
  // check which subsystem hit is in
  const uint16_t layer = TrkrDefs::getLayer(setKey);
  const bool isMvtx = IsInMvtx(layer);
  const bool isIntt = IsInIntt(layer);
  const bool isTpc = IsInTpc(layer);

  // get phi and z values
  //   - FIXME should be more explicit about
  //     row/column vs. z/phi...
  uint16_t phiBin = std::numeric_limits<uint16_t>::max();
  uint16_t zBin = std::numeric_limits<uint16_t>::max();
  if (isMvtx)
  {
    phiBin = MvtxDefs::getCol(hitKey);
    zBin = MvtxDefs::getRow(hitKey);
  }
  else if (isIntt)
  {
    phiBin = InttDefs::getCol(hitKey);
    zBin = InttDefs::getRow(hitKey);
  }
  else if (isTpc)
  {
    phiBin = TpcDefs::getPad(hitKey);
    /* TODO put in z calculation */
  }

  // collect hit info
  HitQAContent content{
      .ene = hit->getEnergy(),
      .adc = hit->getAdc(),
      .layer = layer,
      .phiBin = phiBin,
      .zBin = zBin};

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
}  // end 'GetInfo(TrkrHit*, TrkrDefs::hitsetkey&, TrkrDefs::hitkey&)'

// private methods ============================================================

// ----------------------------------------------------------------------------
//! Fill tracker hit histograms
// ----------------------------------------------------------------------------
void TrksInJetQAHitManager::FillHistograms(const int type, HitQAContent& content)
{
  // fill 1d histograms
  m_mapHist1D[Index(type, H1D::Ene)]->Fill(content.ene);
  m_mapHist1D[Index(type, H1D::ADC)]->Fill(content.adc);
  m_mapHist1D[Index(type, H1D::Layer)]->Fill(content.layer);
  m_mapHist1D[Index(type, H1D::PhiBin)]->Fill(content.phiBin);
  m_mapHist1D[Index(type, H1D::ZBin)]->Fill(content.zBin);

  // fill 2d histograms
  m_mapHist2D[Index(type, H2D::EneVsLayer)]->Fill(content.layer, content.ene);
  m_mapHist2D[Index(type, H2D::EneVsADC)]->Fill(content.adc, content.ene);
  m_mapHist2D[Index(type, H2D::PhiVsZBin)]->Fill(content.zBin, content.phiBin);
}  //  end 'FillHistograms(Type, HitQAContent&)'

// ----------------------------------------------------------------------------
//! Define tracker hit histograms
// ----------------------------------------------------------------------------
void TrksInJetQAHitManager::DefineHistograms()
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
  m_mapHistDef1D[H1D::Ene] = std::tuple("HitEne", vecBins.at(TrksInJetQAHist::Var::Ene));
  m_mapHistDef1D[H1D::ADC] = std::tuple("HitAdc", vecBins.at(TrksInJetQAHist::Var::Adc));
  m_mapHistDef1D[H1D::Layer] = std::tuple("HitLayer", vecBins.at(TrksInJetQAHist::Var::Layer));
  m_mapHistDef1D[H1D::PhiBin] = std::tuple("HitPhiBin", vecBins.at(TrksInJetQAHist::Var::PhiBin));
  m_mapHistDef1D[H1D::ZBin] = std::tuple("HitZBin", vecBins.at(TrksInJetQAHist::Var::ZBin));

  // define 2d histograms
  m_mapHistDef2D[H2D::EneVsLayer] = std::tuple("HitEneVsLayer",
                                               vecBins.at(TrksInJetQAHist::Var::Layer),
                                               vecBins.at(TrksInJetQAHist::Var::Ene));
  m_mapHistDef2D[H2D::EneVsADC] = std::tuple("HitEneVsADC",
                                             vecBins.at(TrksInJetQAHist::Var::Adc),
                                             vecBins.at(TrksInJetQAHist::Var::Ene));
  m_mapHistDef2D[H2D::PhiVsZBin] = std::tuple("HitPhiVsZBin",
                                              vecBins.at(TrksInJetQAHist::Var::ZBin),
                                              vecBins.at(TrksInJetQAHist::Var::PhiBin));
}  // end 'DefineHistograms()'

// end ========================================================================
