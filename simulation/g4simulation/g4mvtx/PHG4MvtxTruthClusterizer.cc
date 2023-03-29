#include "PHG4MvtxTruthClusterizer.h"
#include "PHG4MvtxDigitizer.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <g4tracking/TrkrTruthTrack.h>
#include <g4tracking/TrkrTruthTrackContainer.h>
#include <iostream>
#include <mvtx/MvtxClusterizer.h>
#include <mvtx/MvtxHitPruner.h>
#include <phool/PHIODataNode.h>
#include <trackbase/TrkrClusterContainerv4.h>
#include <trackbase/TrkrHitSetContainerv1.h>

PHG4MvtxTruthClusterizer::PHG4MvtxTruthClusterizer (
        PHG4MvtxDigitizer* _digitiser
      , MvtxHitPruner*     _pruner
      , MvtxClusterizer*   _clusterizer
      , int                _verbosity ) 
  : TruthClusterizerBase (_verbosity)
  , m_digitiser   { _digitiser   }
  , m_pruner      { _pruner      }
  , m_clusterizer { _clusterizer }
{}

void PHG4MvtxTruthClusterizer::init_run(PHCompositeNode*& _topNode) {
  m_topNode = _topNode;
  init_nodes(_topNode); // base class TruthClusterizerBase gets required nodes
}

int PHG4MvtxTruthClusterizer::clusterize_hits(TrkrClusterContainer* clusters)
{
  m_digitiser ->DigitizeMvtxLadderCells(m_topNode, m_hits, true, m_verbosity);

  int prune = m_pruner ->process_TrkrHitSetContainer(m_hits, m_verbosity);
  if (prune != Fun4AllReturnCodes::EVENT_OK) return prune;

  int cluster = m_clusterizer ->ClusterMvtx(m_topNode, m_hits, clusters, true, m_verbosity);
  if (cluster != Fun4AllReturnCodes::EVENT_OK) return prune;

  return Fun4AllReturnCodes::EVENT_OK;
}

void PHG4MvtxTruthClusterizer::check_g4hit(PHG4Hit* hit) {
  if (m_verbosity>10) std::cout << " -> Checking PHG4Hit" << std::endl;
  check_g4hit_status(hit);
  if (m_was_emb) {
    if (m_verbosity>3) {
      std::cout << PHWHERE << std::endl 
      << " -> Pre clustering " << (int) m_hits->size() << " hits" << std::endl;
    }
    TrkrClusterContainerv4 clusters{};
    clusterize_hits   (&clusters);
    transfer_clusters (&clusters);
    if (m_verbosity>3) {
      std::cout << PHWHERE << std::endl 
      << " -> Clustered " << (int) clusters.size() << " clusters" << std::endl;
    }
  }
  if (m_is_new_track) update_track();
}

void PHG4MvtxTruthClusterizer::end_of_event() {
  check_g4hit(nullptr); // flush out last data if ended in truth track
  m_hitsetkey_cnt.clear();
  if (m_verbosity>2) { 
      std::cout << PHWHERE << " :: tracks with clusters after clustering in MVTX" << std::endl;
      for (auto& track : m_truthtracks->getMap()) {
        std::cout << "  track("<< track.first <<") nclusters: " << track.second->getClusters().size();
        for (auto& cluster : track.second->getClusters()) std::cout << " " << (int) TrkrDefs::getLayer(cluster);
        std::cout << std::endl;
      }
  }
}

