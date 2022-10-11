#include "TruthRecoTrackMatching.h"

#include <trackbase/TrkrTruthTrackContainer.h>
#include <trackbase/TrkrTruthTrackContainerv1.h>
#include <trackbase/TrkrTruthTrackv1.h>
#include <trackbase_historic/SvtxTrackMap.h>

#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterv4.h>

#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterContainerv4.h>

// not actually sure if I need all of these
#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>  // for PHDataNode
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>  // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>  // for PHObject
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <g4main/PHG4TruthInfoContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <TSystem.h>

/* #include <vector> */

int TruthRecoTrackMatching::InitRun(PHCompositeNode *topNode)
{ 
  if (createNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
  {
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}


int TruthRecoTrackMatching::process_event(PHCompositeNode* topnode) 
{
  if (topnode == nullptr) {
    std::cout << " topnode is null " << std::endl;
  }

  // FIXME orange 0
  std::cout << " FIXME orange: Thing in the top node: " << std::endl;
  topnode->print();
  std::cout << " processing event in TruthRecoTrackMatching " << std::endl;

  // make a matching to the "best track" number
  // There are 
  // now just loop through the TrkrTruthTrackContainer
  std::vector<bool> matched_truthtracks (m_TrkrTruthTrackContainer->getTruthTracks().size(), false);
  std::vector<bool> matched_recotracks  (m_SvtxTrackMap->size(), false);
  
  std::cout << "id-C size: " << m_TrkrTruthTrackContainer->getTruthTracks().size() << std::endl;
  for (auto _track : m_TrkrTruthTrackContainer->getTruthTracks()) {
    auto track = static_cast<TrkrTruthTrackv1*>(_track);
    /* std::cout << " id-C: " << track->getTrackid() << std::endl; */
    std::cout << Form(" id-C (%2i)   phi(%5.2f) eta(%5.2f) pt(%5.2f)",
        track->getTrackid(), track->getPhi(), track->getPseudoRapidity(), track->getPt()) << std::endl;
  }

  std::cout << "id size: " << m_SvtxTrackMap->size() << std::endl;
  for (auto reco = m_SvtxTrackMap->begin(); reco != m_SvtxTrackMap->end(); ++reco) {
    auto _ = reco->second;
    std::cout << Form(" id(%2i-%2i)   phi(%5.2f) eta(%5.2f) pt(%5.2f)",
        reco->first, _->get_id(), _->get_phi(), _->get_eta(), _->get_pt()) << std::endl;
  }


  /* std::cout << " FIXME loop 3 " << std::endl; */
  /* auto particle_map = m_PHG4TruthInfoContainer->GetMap(); */

  /* std::cout << " FIXME loop 0 " << std::endl; */
  /* for (const auto& pair : particle_map) { */
  /*   int id = pair.first; */
  /*   if (!m_PHG4TruthInfoContainer->isEmbeded(id)) continue; */

  /*   // do the matching -> find the best reconstructed track */
  /*   std::cout << " id-A: " << id << "  is embedded? " << m_PHG4TruthInfoContainer->isEmbeded(id) << std::endl; */
  /*   // get the particle eta, phi, and pt; must match the to the reconstructed tracks */
  /* } */
  /* std::cout << " FIXME loop 1 " << std::endl; */
  /* // test loop directly through the embedded track ids */
  /* auto embTrkIds = m_PHG4TruthInfoContainer->GetEmbeddedTrkIds(); */
  /* for (auto ids = embTrkIds.first; ids != embTrkIds.second; ++ids) { */
  /*   std::cout << " id-B: " << ids->first << "  is embedded? " << m_PHG4TruthInfoContainer->isEmbeded(ids->first) << std::endl; */
  /* } */
  /* std::cout << " FIXME loop 2 " << std::endl; */


  return Fun4AllReturnCodes::EVENT_OK;
}

int TruthRecoTrackMatching::End(PHCompositeNode * /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}


int TruthRecoTrackMatching::createNodes(PHCompositeNode* topNode)
{
  // Initiailize the modules
  //
  //
  //
  //
  m_PHG4TruthInfoContainer = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!m_PHG4TruthInfoContainer)
  {
    std::cout << "Could not locate G4TruthInfo node when running \"TruthRecoTrackMatching\" module." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_SvtxTrackMap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if (!m_SvtxTrackMap)
  {
    std::cout << "Could not locate SvtxTrackMap node when running \"TruthRecoTrackMatching\" module." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }


  m_TruthClusterContainer = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_TRUTHCLUSTERCONTAINER");
  if (!m_TruthClusterContainer)
  {
    std::cout << "Could not locate TRKR_TRUTHCLUSTERCONTAINER node when running \"TruthRecoTrackMatching\" module." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_RecoClusterContainer = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!m_RecoClusterContainer)
  {
    std::cout << "Could not locate TRKR_CLUSTER node when running \"TruthRecoTrackMatching\" module." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_TrkrTruthTrackContainer = findNode::getClass<TrkrTruthTrackContainer>(topNode, "TRKR_TRUTHTRACKCONTAINER");
  if (!m_TrkrTruthTrackContainer)
  {
    std::cout << "Could not locate TRKR_TRUTHTRACKCONTAINER node when running \"TruthRecoTrackMatching\" module." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
