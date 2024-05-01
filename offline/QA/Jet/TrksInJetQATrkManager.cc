// ----------------------------------------------------------------------------
// 'TrksInJetQATrkManager.cc'
// Derek Anderson
// 03.25.2024
//
// A submodule for the TrksInJetQA module to generate
// QA plots for tracks
// ----------------------------------------------------------------------------

#define TRKSINJETQATRKMANAGER_CC

// submodule definition
#include "TrksInJetQATrkManager.h"



// public methods -------------------------------------------------------------

void TrksInJetQATrkManager::GetInfo(SvtxTrack* track) {

  // collect track info
  TrackQAContent content = {
    .eta  = track -> get_eta(),
    .phi  = track -> get_phi(),
    .pt   = track -> get_pt(),
    .qual = track -> get_quality()
  };

  // fill histograms
  FillHistograms(Type::All, content);
  return;

}  // end 'GetInfo(SvtxTrack*)'



// private methods ------------------------------------------------------------

void TrksInJetQATrkManager::FillHistograms(const int type, TrackQAContent& content) {

  // fill 1d histograms
  m_vecHist1D.at(type).at(H1D::Eta)  -> Fill(content.eta);
  m_vecHist1D.at(type).at(H1D::Phi)  -> Fill(content.phi);
  m_vecHist1D.at(type).at(H1D::Pt)   -> Fill(content.pt);
  m_vecHist1D.at(type).at(H1D::Qual) -> Fill(content.qual);

  // fill 2d histograms
  m_vecHist2D.at(type).at(H2D::EtaVsPhi) -> Fill(content.phi, content.eta);
  m_vecHist2D.at(type).at(H2D::PtVsQual) -> Fill(content.qual, content.pt);
  return;

}  //  end 'FillHistograms(Type, TrackQAContent&)'



void TrksInJetQATrkManager::DefineHistograms() {

  // grab binning schemes
  std::vector<BinDef> vecBins = m_hist.GetVecHistBins();

  // set histogram types
  m_vecHistTypes.push_back( "All" );

  // define 1d histograms
  m_vecHistDef1D.push_back( std::make_tuple( "TrackEta",  vecBins.at(TrksInJetQAHist::Var::Eta)  ));
  m_vecHistDef1D.push_back( std::make_tuple( "TrackPhi",  vecBins.at(TrksInJetQAHist::Var::Phi)  ));
  m_vecHistDef1D.push_back( std::make_tuple( "TrackPt",   vecBins.at(TrksInJetQAHist::Var::Ene)  ));
  m_vecHistDef1D.push_back( std::make_tuple( "TrackQual", vecBins.at(TrksInJetQAHist::Var::Qual) ));

  // define 2d histograms
  m_vecHistDef2D.push_back(
    std::make_tuple(
      "TrackEtaVsPhi",
      vecBins.at(TrksInJetQAHist::Var::Phi),
      vecBins.at(TrksInJetQAHist::Var::Eta)
    )
  );
  m_vecHistDef2D.push_back(
    std::make_tuple(
      "TrackPtVsQual",
      vecBins.at(TrksInJetQAHist::Var::Qual),
      vecBins.at(TrksInJetQAHist::Var::Ene)
    )
  );
  return;

}  // end 'BuildHistograms()'

// end ------------------------------------------------------------------------
