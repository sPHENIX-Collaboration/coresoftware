// ----------------------------------------------------------------------------
// 'TrksInJetQAInJetFiller.cc'
// Derek Anderson
// 04.03.2024
//
// A submodule for the TrksInJetsQAM F4A module to produce
// QA histograms for tracks and more in jets
// ----------------------------------------------------------------------------

#define TRKSINJETQAINJETFILLER_CC

// submodule definition
#include "TrksInJetQAInJetFiller.h"



// inherited public methods ---------------------------------------------------

void TrksInJetQAInJetFiller::Fill(PHCompositeNode* topNode) {

  GetNodes(topNode);

  FillJetAndTrackQAHists(topNode);
  return;

}  // end 'Fill(PHCompositeNode* topNode)'



// private methods ------------------------------------------------------------

void TrksInJetQAInJetFiller::GetNode(const int node, PHCompositeNode* topNode) {

  // jump to relevant node
  switch (node) {

    case Node::Flow:
      m_flowStore = findNode::getClass<ParticleFlowElementContainer>(topNode, "ParticleFlowElements");
      if (!m_flowStore) {
        std::cerr << PHWHERE << ": PANIC: Couldn't grab particle flow container from node tree!" << std::endl;
        assert(m_flowStore);
      }
      break;

    default:
      std::cerr << PHWHERE << ": WARNING: trying to grab unkown additional node..." << std::endl;
      break;
  }
  return;

}  // end 'GetNode(int, PHCompositeNode*)'



void TrksInJetQAInJetFiller::FillJetAndTrackQAHists(PHCompositeNode* topNode) {

  // loop over jets
  for (
    uint64_t iJet = 0;
    iJet < m_jetMap -> size();
    ++iJet
  ) {

    // grab jet and make sure track vector is clear
    Jet* jet = m_jetMap -> get_jet(iJet);
    m_trksInJet.clear();

    // get all tracks "in" a jet
    GetCstTracks(jet, topNode);
    GetNonCstTracks(jet, topNode);

    // grab jet info and fill histograms
    if (m_config.doJetQA) m_jetManager -> GetInfo(jet, m_trksInJet);

    // loop over tracks in the jet
    for (SvtxTrack* track : m_trksInJet) {

      // grab track info and fill histograms
      if (m_config.doTrackQA) m_trackManager -> GetInfo(track);

      // fill cluster and hit histograms as needed
      if (m_config.doClustQA || m_config.doHitQA) {
        FillClustAndHitQAHists(track, topNode);
      }
    }  // end track loop
  }  // end jet loop
  return;

}  // end 'FillJetAndTrackQAHists(PHCompositeNode*)'



void TrksInJetQAInJetFiller::FillClustAndHitQAHists(SvtxTrack* track, PHCompositeNode* topNode) {

  // get cluster keys
  for (auto clustKey : ClusKeyIter(track)) {

    // grab cluster and its info
    if (m_config.doClustQA) {
      m_clustManager -> GetInfo(
        m_clustMap -> findCluster(clustKey),
        clustKey,
        m_actsGeom
      );
    }

    // get hits if needed
    if (m_config.doHitQA) {

      // grab hit set and key associated with cluster key
      TrkrDefs::hitsetkey setKey = TrkrDefs::getHitSetKeyFromClusKey(clustKey);
      TrkrHitSet*         set    = m_hitMap -> findHitSet(setKey);
      if (!set || !(set -> size() > 0)) return;

      // loop over all hits in hit set
      TrkrHitSet::ConstRange hits = set -> getHits();
      for (
        TrkrHitSet::ConstIterator itHit = hits.first;
        itHit != hits.second;
        ++itHit
      ) {

        // grab hit
        TrkrDefs::hitkey hitKey = itHit -> first;
        TrkrHit*         hit    = itHit -> second;

        // grab info and fill histograms
        m_hitManager -> GetInfo(hit, setKey, hitKey);

      }  // end hit loop
    }
  }  // end cluster key loop
  return;

}  // end 'FillClustQAHists(SvtxTrack*, PHCompositeNode*)'



void TrksInJetQAInJetFiller::GetCstTracks(Jet* jet, PHCompositeNode* topNode) {

  // loop over consituents
  auto csts = jet -> get_comp_vec();
  for (
    auto cst = csts.begin();
    cst != csts.end();
    ++cst
  ) {

    // ignore cst if non-relevent type
    const uint32_t src = cst -> first;
    if ( IsCstNotRelevant(src) ) continue;

    // if cst is track, add to list
    if (src == Jet::SRC::TRACK) {
      m_trksInJet.push_back( m_trkMap -> get(cst -> second) );
    }

    // if pfo, grab track if needed
    if (src == Jet::SRC::PARTICLE) {
      PFObject*  pfo   = GetPFObject(cst -> second, topNode);
      SvtxTrack* track = GetTrkFromPFO(pfo);
      if (track) {
        m_trksInJet.push_back( track );
      }
    }
  }  // end cst loop
  return;

}  // end 'GetCstTracks(Jet* jet, PHCompositeNode* topNode)'



void TrksInJetQAInJetFiller::GetNonCstTracks(Jet* jet, PHCompositeNode* topNode) {

  // loop over tracks
  for (
    SvtxTrackMap::Iter itTrk = m_trkMap -> begin();
    itTrk != m_trkMap -> end();
    ++itTrk
  ) {

    // grab track
    SvtxTrack* track = itTrk -> second;

    // ignore tracks we've already added to the list
    if ( IsTrkInList(track -> get_id()) ) continue;

    // FIXME this can be improved!
    //   - jets don't necessarily have areas of
    //     pi*(Rjet)^2
    //   - it may be better to instead check
    //     if a track projection falls in
    //     a constituent tower/cluster
    //   - Also track projections to a calo
    //     would be better to use than just
    //     the track
    const double dr = GetTrackJetDist(track, jet);
    if (dr < m_config.rJet) {
      m_trksInJet.push_back( track );
    }
  }  // end track loop
  return;

}  // end 'GetNonCstTracks(Jet* jet, PHCompositeNode* topNode)'



bool TrksInJetQAInJetFiller::IsCstNotRelevant(const uint32_t type) {

  const bool isVoid   = (type == Jet::SRC::VOID);
  const bool isImport = (type == Jet::SRC::HEPMC_IMPORT);
  const bool isProbe  = (type == Jet::SRC::JET_PROBE);
  return (isVoid || isImport || isProbe);

}  // end 'IsCstNotRelevant(uint32_t)'



bool TrksInJetQAInJetFiller::IsTrkInList(const uint32_t id) {

  bool isAdded = false;
  for (SvtxTrack* trkInJet : m_trksInJet) {
    if (id == trkInJet -> get_id() ) {
      isAdded = true;
      break;
    }
  }
  return isAdded;

}  // end 'IsTrkInList(uint32_t)'



double TrksInJetQAInJetFiller::GetTrackJetDist(SvtxTrack* track, Jet* jet) {

  // get delta eta
  const double dEta = (track -> get_eta()) - (jet -> get_eta());

  // get delta phi
  double dPhi = (track -> get_phi()) - (jet -> get_phi());
  if (dPhi < (-1. * TMath::Pi())) dPhi += TMath::TwoPi();
  if (dPhi > (1. * TMath::Pi()))  dPhi -= TMath::TwoPi();

  // return distance
  return std::hypot(dEta, dPhi);

}  // end 'GetTrackJetDist(SvtxTrack*, Jet*)'



PFObject* TrksInJetQAInJetFiller::GetPFObject(const uint32_t id, PHCompositeNode* topNode) {

  // pointer to pfo 
  PFObject* pfoToFind = NULL;

  // grab pf node if needed
  if (!m_flowStore) GetNode(Node::Flow, topNode);

  // loop over pfos
  for (
      ParticleFlowElementContainer::ConstIterator itFlow = m_flowStore -> getParticleFlowElements().first;
      itFlow != m_flowStore -> getParticleFlowElements().second;
      ++itFlow
  ) {

    // get pfo
    PFObject* pfo = itFlow -> second;

    // if has provided id, set pointer and exit
    if (id == pfo -> get_id()) {
      pfoToFind = pfo;
      break;
    }
  }  // end pfo loop
  return pfoToFind;

}  // end 'GetPFObject(uint32_t, PHCompositeNode*)'



SvtxTrack* TrksInJetQAInJetFiller::GetTrkFromPFO(PFObject* pfo) {

  // pointer to track
  SvtxTrack* track = NULL;

  // if pfo has track, try to grab
  const auto type = pfo -> get_type();
  if (
    (type == ParticleFlowElement::PFLOWTYPE::MATCHED_CHARGED_HADRON) ||
    (type == ParticleFlowElement::PFLOWTYPE::UNMATCHED_CHARGED_HADRON)
  ) {
    track = pfo -> get_track();
  }
  return track;

}  // end 'GetTrkFromPFO(PFObject*)'

// end ------------------------------------------------------------------------
