/// ===========================================================================
/*! \file   TrksInJetQAInJetFiller.cc
 *  \author Derek Anderson
 *  \date   04.03.2024
 *
 *  A submodule for the TrksInJetsQAM F4A module to produce
 *  QA histograms for tracks and more in jets
 */
/// ===========================================================================


// submodule definition
#include "TrksInJetQAInJetFiller.h"

// inherited public methods ===================================================

// ----------------------------------------------------------------------------
//! Run fill for relevant histograms
// ----------------------------------------------------------------------------
void TrksInJetQAInJetFiller::Fill(PHCompositeNode* topNode)
{
  GetNodes(topNode);

  FillJetAndTrackQAHists(topNode);
}  // end 'Fill(PHCompositeNode* topNode)'

// private methods ============================================================

// ----------------------------------------------------------------------------
//! Grab additional input nodes
// ----------------------------------------------------------------------------
void TrksInJetQAInJetFiller::GetNode(const int node, PHCompositeNode* topNode)
{
  // jump to relevant node
  // NOLINTNEXTLINE(hicpp-multiway-paths-covered)
  switch (node)
  {
    case Node::Flow:
      m_flowStore = findNode::getClass<ParticleFlowElementContainer>(topNode, "ParticleFlowElements");
      if (!m_flowStore)
      {
        std::cerr << PHWHERE << ": PANIC: Couldn't grab particle flow container from node tree!" << std::endl;
        assert(m_flowStore);
      }
      break;
    default:
      std::cerr << PHWHERE << ": WARNING: trying to grab unkown additional node..." << std::endl;
      break;
  }
}  // end 'GetNode(int, PHCompositeNode*)'

// ----------------------------------------------------------------------------
//! Fill histograms for jet and in-jet track histograms
// ----------------------------------------------------------------------------
/*! This defines how to (1) loop over jets, and (2) extract
 *  tracks inside them. Some care is required here since
 *  jets can be made
 *    a. completely without tracks (e.g. from calo towers),
 *    b. completely with tracks (i.e. track jets),
 *    c. or with a mix of tracks and other objects (e.g.
 *       particle flow jets).
 *  Once tracks have been identified, the relevant hit and
 *  cluster populations are extracted from the tracks, and
 *  their histograms are filled.
 */
void TrksInJetQAInJetFiller::FillJetAndTrackQAHists(PHCompositeNode* topNode)
{
  // loop over jets
  for (
      uint64_t iJet = 0;
      iJet < m_jetMap->size();
      ++iJet)
  {
    // grab jet and make sure track vector is clear
    Jet* jet = m_jetMap->get_jet(iJet);
    m_trksInJet.clear();

    // get all tracks "in" a jet
    GetCstTracks(jet, topNode);
    GetNonCstTracks(jet);

    // grab jet info and fill histograms
    if (m_config.doJetQA)
    {
      m_jetManager->GetInfo(jet, m_trksInJet);
    }

    // loop over tracks in the jet
    for (SvtxTrack* track : m_trksInJet)
    {
      // grab track info and fill histograms
      if (m_config.doTrackQA)
      {
        m_trackManager->GetInfo(track, jet);
      }

      // fill cluster and hit histograms as needed
      if (m_config.doClustQA || m_config.doHitQA)
      {
        FillClustAndHitQAHists(track);
      }
    }  // end track loop
  }  // end jet loop
}  // end 'FillJetAndTrackQAHists(PHCompositeNode*)'

// ----------------------------------------------------------------------------
//! Fill histograms for in-jet track hits and clusters
// ----------------------------------------------------------------------------
void TrksInJetQAInJetFiller::FillClustAndHitQAHists(SvtxTrack* track)
{
  // get cluster keys
  for (auto clustKey : ClusKeyIter(track))
  {
    // grab cluster and its info
    if (m_config.doClustQA)
    {
      m_clustManager->GetInfo(
          m_clustMap->findCluster(clustKey),
          clustKey,
          m_actsGeom);
    }

    // get hits if needed
    if (m_config.doHitQA)
    {
      // grab hit set and key associated with cluster key
      TrkrDefs::hitsetkey setKey = TrkrDefs::getHitSetKeyFromClusKey(clustKey);
      TrkrHitSet* set = m_hitMap->findHitSet(setKey);
      if (!set || !(set->size() > 0))
      {
        return;
      }

      // loop over all hits in hit set
      TrkrHitSet::ConstRange hits = set->getHits();
      for (
          TrkrHitSet::ConstIterator itHit = hits.first;
          itHit != hits.second;
          ++itHit)
      {
        // grab hit
        TrkrDefs::hitkey hitKey = itHit->first;
        TrkrHit* hit = itHit->second;

        // grab info and fill histograms
        m_hitManager->GetInfo(hit, setKey, hitKey);

      }  // end hit loop
    }
  }  // end cluster key loop
}  // end 'FillClustQAHists(SvtxTrack*)'

// ----------------------------------------------------------------------------
//! Identify constituent tracks in a jet
// ----------------------------------------------------------------------------
void TrksInJetQAInJetFiller::GetCstTracks(Jet* jet, PHCompositeNode* topNode)
{
  // loop over consituents
  auto csts = jet->get_comp_vec();
  for (auto& cst : csts)
  {
    // ignore cst if non-relevent type
    const uint32_t src = cst.first;
    if (IsCstNotRelevant(src))
    {
      continue;
    }

    // if cst is track, add to list
    if (src == Jet::SRC::TRACK)
    {
      m_trksInJet.push_back(m_trkMap->get(cst.second));
    }

    // if pfo, grab track if needed
    if (src == Jet::SRC::PARTICLE)
    {
      TrksInJetQADefs::PFObject* pfo = GetPFObject(cst.second, topNode);
      SvtxTrack* track = GetTrkFromPFO(pfo);
      if (track)
      {
        m_trksInJet.push_back(track);
      }
    }
  }  // end cst loop
}  // end 'GetCstTracks(Jet* jet, PHCompositeNode* topNode)'

// ----------------------------------------------------------------------------
//! Identify tracks in jet cone that are NOT constituents
// ----------------------------------------------------------------------------
/*! This method comes into play if, e.g., the jets being
 *  analyzed were made with calorimeter information only.
 *  In this case you can still have tracks lying inside
 *  the jet cone despite them not being consituents of
 *  the constructed jet.
 */
void TrksInJetQAInJetFiller::GetNonCstTracks(Jet* jet)
{
  // loop over tracks
  for (auto& itTrk : *m_trkMap)
  {
    // grab track
    SvtxTrack* track = itTrk.second;

    // ignore tracks we've already added to the list
    if (IsTrkInList(track->get_id()))
    {
      continue;
    }

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
    if (dr < m_config.rJet)
    {
      m_trksInJet.push_back(track);
    }
  }  // end track loop
}  // end 'GetNonCstTracks(Jet* jet)'

// ----------------------------------------------------------------------------
//! Check if a consituent is not a relevant type
// ----------------------------------------------------------------------------
bool TrksInJetQAInJetFiller::IsCstNotRelevant(const uint32_t type)
{
  const bool isVoid = (type == Jet::SRC::VOID);
  const bool isImport = (type == Jet::SRC::HEPMC_IMPORT);
  const bool isProbe = (type == Jet::SRC::JET_PROBE);
  return (isVoid || isImport || isProbe);
}  // end 'IsCstNotRelevant(uint32_t)'

// ----------------------------------------------------------------------------
//! Check if a track has already been added to the list of tracks in the jet
// ----------------------------------------------------------------------------
bool TrksInJetQAInJetFiller::IsTrkInList(const uint32_t id)
{
  bool isAdded = false;
  for (SvtxTrack* trkInJet : m_trksInJet)
  {
    if (id == trkInJet->get_id())
    {
      isAdded = true;
      break;
    }
  }
  return isAdded;
}  // end 'IsTrkInList(uint32_t)'

// ----------------------------------------------------------------------------
//! Get eta-phi distance between track & jet axis
// ----------------------------------------------------------------------------
double TrksInJetQAInJetFiller::GetTrackJetDist(SvtxTrack* track, Jet* jet)
{
  // get delta eta
  const double dEta = (track->get_eta()) - (jet->get_eta());

  // get delta phi
  double dPhi = (track->get_phi()) - (jet->get_phi());
  if (dPhi < (-1. * M_PI))
  {
    dPhi += 2. * M_PI;
  }
  if (dPhi > (1. * M_PI))
  {
    dPhi -= 2. * M_PI;
  }

  // return distance
  return std::hypot(dEta, dPhi);

}  // end 'GetTrackJetDist(SvtxTrack*, Jet*)'

// ----------------------------------------------------------------------------
//! Find a particle flow element in node
// ----------------------------------------------------------------------------
TrksInJetQADefs::PFObject* TrksInJetQAInJetFiller::GetPFObject(const uint32_t id,
                                                               PHCompositeNode* topNode)
{
  // pointer to pfo
  TrksInJetQADefs::PFObject* pfoToFind = nullptr;

  // grab pf node if needed
  if (!m_flowStore)
  {
    GetNode(Node::Flow, topNode);
  }

  // loop over pfos
  for (
      ParticleFlowElementContainer::ConstIterator itFlow = m_flowStore->getParticleFlowElements().first;
      itFlow != m_flowStore->getParticleFlowElements().second;
      ++itFlow)
  {
    // get pfo
    TrksInJetQADefs::PFObject* pfo = itFlow->second;

    // if has provided id, set pointer and exit
    if (id == pfo->get_id())
    {
      pfoToFind = pfo;
      break;
    }
  }  // end pfo loop
  return pfoToFind;
}  // end 'GetPFObject(uint32_t, PHCompositeNode*)'

// ----------------------------------------------------------------------------
//! Get track from a particle flow element
// ----------------------------------------------------------------------------
SvtxTrack* TrksInJetQAInJetFiller::GetTrkFromPFO(TrksInJetQADefs::PFObject* pfo)
{
  // pointer to track
  SvtxTrack* track = nullptr;

  // if pfo has track, try to grab
  const auto type = pfo->get_type();
  if (
      (type == ParticleFlowElement::PFLOWTYPE::MATCHED_CHARGED_HADRON) ||
      (type == ParticleFlowElement::PFLOWTYPE::UNMATCHED_CHARGED_HADRON))
  {
    track = pfo->get_track();
  }
  return track;
}  // end 'GetTrkFromPFO(TrksInJetQADefs::PFObject*)'

// end ========================================================================
