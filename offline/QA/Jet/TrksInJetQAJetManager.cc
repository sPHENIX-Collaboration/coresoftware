/// ===========================================================================
/*! \file   TrksInJetQAJetManager.cc
 *  \author Derek Anderson
 *  \date   03.25.2024
 *
 *  A submodule for the TrksInJetQA module
 *  to generate QA plots for jets
 */
/// ===========================================================================

#define TRKSINJETQAJETMANAGER_CC

// submodule definition
#include "TrksInJetQAJetManager.h"

// public methods =============================================================

// ----------------------------------------------------------------------------
//! Get information from a jet
// ----------------------------------------------------------------------------
void TrksInJetQAJetManager::GetInfo(Jet* jet,
                                    std::optional<std::vector<SvtxTrack*>> tracks)
{
  // collect jet info
  JetQAContent jetContent{
      .eta = jet->get_eta(),
      .phi = jet->get_phi(),
      .pt = jet->get_pt(),
      .nTrk = 0,
      .ptSum = 0.};

  // if tracks are available, do calculations
  if (tracks.has_value())
  {
    for (SvtxTrack* track : tracks.value())
    {
      jetContent.nTrk++;
      jetContent.ptSum += track->get_pt();

    }  // end track loop
  }  // end if tracks.has_value

  // fill histograms
  FillHistograms(Type::All, jetContent);
}  // end 'GetInfo(Jet*, std::optional<std::vector<SvtxTrack*>>)'

// private methods ============================================================

// ----------------------------------------------------------------------------
//! Fill jet histograms
// ----------------------------------------------------------------------------
void TrksInJetQAJetManager::FillHistograms(const int type, JetQAContent& content)
{
  // fill 1d histograms
  if (m_config.doOptHist)
  {
    m_mapHist1D[Index(type, H1D::Eta)]->Fill(content.eta);
    m_mapHist1D[Index(type, H1D::Pt)]->Fill(content.pt);
    m_mapHist1D[Index(type, H1D::PtSum)]->Fill(content.ptSum);
  }
  m_mapHist1D[Index(type, H1D::Phi)]->Fill(content.phi);
  m_mapHist1D[Index(type, H1D::NTrk)]->Fill(content.nTrk);

  // fill 2d histograms
  if (m_config.doOptHist)
  {
    m_mapHist2D[Index(type, H2D::PtSumVsPt)]->Fill(content.pt, content.ptSum);
  }
  m_mapHist2D[Index(type, H2D::PtVsEta)]->Fill(content.eta, content.pt);
  m_mapHist2D[Index(type, H2D::NTrkVsEta)]->Fill(content.eta, content.nTrk);
  m_mapHist2D[Index(type, H2D::NTrkVsPt)]->Fill(content.pt, content.nTrk);
}  //  end 'FillHistograms(Type, JetQAContent&)'

// ----------------------------------------------------------------------------
//! Define jet histograms
// ----------------------------------------------------------------------------
void TrksInJetQAJetManager::DefineHistograms()
{
  // grab binning schemes
  std::vector<TrksInJetQADefs::BinDef> vecBins = m_hist.GetVecHistBins();

  // histogram labels
  m_mapHistTypes[Type::All] = "All";

  // 1d histogram definitions
  if (m_config.doOptHist)
  {
    m_mapHistDef1D[H1D::Eta] = std::tuple("JetEta", vecBins.at(TrksInJetQAHist::Var::Eta));
    m_mapHistDef1D[H1D::Pt] = std::tuple("JetPt", vecBins.at(TrksInJetQAHist::Var::Ene));
    m_mapHistDef1D[H1D::PtSum] = std::tuple("SumTrkPt", vecBins.at(TrksInJetQAHist::Var::Ene));
  }
  m_mapHistDef1D[H1D::Phi] = std::tuple("JetPhi", vecBins.at(TrksInJetQAHist::Var::Phi));
  m_mapHistDef1D[H1D::NTrk] = std::tuple("JetNTrks", vecBins.at(TrksInJetQAHist::Var::Num));

  // 2d histogram definitions
  if (m_config.doOptHist)
  {
    m_mapHistDef2D[H2D::PtSumVsPt] = std::tuple("SumTrkVsJetPt",
                                                vecBins.at(TrksInJetQAHist::Var::Ene),
                                                vecBins.at(TrksInJetQAHist::Var::Ene));
  }
  m_mapHistDef2D[H2D::PtVsEta] = std::tuple("JetPtVsEta",
                                            vecBins.at(TrksInJetQAHist::Var::Eta),
                                            vecBins.at(TrksInJetQAHist::Var::Ene));
  m_mapHistDef2D[H2D::NTrkVsEta] = std::tuple("JetNTrkVsEta",
                                              vecBins.at(TrksInJetQAHist::Var::Eta),
                                              vecBins.at(TrksInJetQAHist::Var::Num));
  m_mapHistDef2D[H2D::NTrkVsPt] = std::tuple("JetNTrkVsPt",
                                             vecBins.at(TrksInJetQAHist::Var::Ene),
                                             vecBins.at(TrksInJetQAHist::Var::Num));
}  // end 'DefineHistograms()'

// end ========================================================================
