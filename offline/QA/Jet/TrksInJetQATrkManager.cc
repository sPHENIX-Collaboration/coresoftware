/// ===========================================================================
/*! \file   TrksInJetQATrkManager.cc
 *  \author Derek Anderson
 *  \date   03.25.2024
 *
 *  A submodule for the TrksInJetQA module to generate
 *  QA plots for tracks
 */
/// ===========================================================================


// submodule definition
#include "TrksInJetQATrkManager.h"

// root libraries
#include <TVector3.h>

// public methods =============================================================

// ----------------------------------------------------------------------------
//! Get information from a track
// ----------------------------------------------------------------------------
void TrksInJetQATrkManager::GetInfo(SvtxTrack* track,
                                    std::optional<Jet*> jet)
{
  // collect track-only info
  TrackQAContent content = {
      .eta = track->get_eta(),
      .phi = track->get_phi(),
      .pt = track->get_pt(),
      .qual = track->get_quality()};

  // if a jet was provided, calculate z and jt
  if (jet.has_value())
  {
    // grab track & jet momenta
    TVector3 pTrk(track->get_px(), track->get_py(), track->get_pz());
    TVector3 pJet(jet.value()->get_px(), jet.value()->get_py(), jet.value()->get_pz());

    // now calculate z & jt
    content.z = pJet.Dot(pTrk) / pJet.Mag2();
    content.jt = pJet.Cross(pTrk).Mag() / pJet.Mag2();
  }

  // fill histograms
  FillHistograms(Type::All, content);
}  // end 'GetInfo(SvtxTrack*, std::optional<Jet*>)'

// private methods ============================================================

// ----------------------------------------------------------------------------
//! Fill track histograms
// ---------------------------------------------------------------------------- 
void TrksInJetQATrkManager::FillHistograms(const int type, TrackQAContent& content)
{
  // fill 1d histograms
  m_mapHist1D[Index(type, H1D::Eta)]->Fill(content.eta);
  m_mapHist1D[Index(type, H1D::Phi)]->Fill(content.phi);
  m_mapHist1D[Index(type, H1D::Pt)]->Fill(content.pt);
  m_mapHist1D[Index(type, H1D::Z)]->Fill(content.z);
  m_mapHist1D[Index(type, H1D::Jt)]->Fill(content.jt);
  m_mapHist1D[Index(type, H1D::Qual)]->Fill(content.qual);

  // fill 2d histograms
  m_mapHist2D[Index(type, H2D::EtaVsPhi)]->Fill(content.phi, content.eta);
  m_mapHist2D[Index(type, H2D::PtVsQual)]->Fill(content.qual, content.pt);
}  //  end 'FillHistograms(Type, TrackQAContent&)'

// ----------------------------------------------------------------------------
//! Define track histograms
// ----------------------------------------------------------------------------
void TrksInJetQATrkManager::DefineHistograms()
{
  // grab binning schemes
  std::vector<TrksInJetQADefs::BinDef> vecBins = m_hist.GetVecHistBins();

  // set histogram types
  m_mapHistTypes[Type::All] = "All";

  // define 1d histograms
  m_mapHistDef1D[H1D::Eta] = std::tuple("TrackEta", vecBins.at(TrksInJetQAHist::Var::Eta));
  m_mapHistDef1D[H1D::Phi] = std::tuple("TrackPhi", vecBins.at(TrksInJetQAHist::Var::Phi));
  m_mapHistDef1D[H1D::Pt] = std::tuple("TrackPt", vecBins.at(TrksInJetQAHist::Var::Ene));
  m_mapHistDef1D[H1D::Z] = std::tuple("TrackZ", vecBins.at(TrksInJetQAHist::Var::Frac));
  m_mapHistDef1D[H1D::Jt] = std::tuple("TrackJt", vecBins.at(TrksInJetQAHist::Var::Frac));
  m_mapHistDef1D[H1D::Qual] = std::tuple("TrackQual", vecBins.at(TrksInJetQAHist::Var::Qual));

  // define 2d histograms
  m_mapHistDef2D[H2D::EtaVsPhi] = std::tuple("TrackEtaVsPhi",
                                             vecBins.at(TrksInJetQAHist::Var::Phi),
                                             vecBins.at(TrksInJetQAHist::Var::Eta));
  m_mapHistDef2D[H2D::PtVsQual] = std::tuple("TrackPtVsQual",
                                             vecBins.at(TrksInJetQAHist::Var::Qual),
                                             vecBins.at(TrksInJetQAHist::Var::Ene));
}  // end 'DefineHistograms()'

// end ========================================================================
