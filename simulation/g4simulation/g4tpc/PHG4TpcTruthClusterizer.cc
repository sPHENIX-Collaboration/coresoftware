#include "PHG4TpcTruthClusterizer.h"
#include "PHG4TpcDigitizer.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <g4tracking/TrkrTruthTrack.h>
#include <g4tracking/TrkrTruthTrackContainer.h>
#include <iostream>
#include <tpc/TpcClusterizer.h>
#include <phool/PHIODataNode.h>
#include <trackbase/TrkrClusterContainerv4.h>
#include <trackbase/TrkrHitSetContainerv1.h>

PHG4TpcTruthClusterizer::PHG4TpcTruthClusterizer (
        PHG4TpcDigitizer* _digitiser
      , TpcClusterizer*   _clusterizer
      , int                _verbosity ) 
  : TruthClusterizerBase (_verbosity)
  , m_digitiser          { _digitiser      }
  , m_clusterizer        { _clusterizer    }
{}

void PHG4TpcTruthClusterizer::init_run(PHCompositeNode*& _topNode) {
  m_topNode = _topNode;
  init_nodes(_topNode); // base class TruthClusterizerBase gets required nodes
}

int PHG4TpcTruthClusterizer::clusterize_hits(TrkrClusterContainer* clusters)
{
  std::cout << (clusters == nulltpr) << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

void PHG4TpcTruthClusterizer::check_g4hit(PHG4Hit* hit) {
  if (m_verbosity>10) std::cout << " -> Checking PHG4Hit in Tpc" << std::endl;
  check_g4hit_status(hit);
  if (m_was_emb) {
    if (m_verbosity>3) {
      std::cout << PHWHERE << std::endl 
      << " -> Pre Tpc Truth("<<m_current_track->getTrackid() << ") clustering " << (int) m_hits->size() << " hits" << std::endl;
      m_hits->identify();
    }
    TrkrClusterContainerv4 clusters{};
    clusterize_hits   (&clusters);
    transfer_clusters (&clusters);
    if (m_verbosity>3) {
      std::cout << PHWHERE << std::endl 
      << " -> Clustered " << (int) clusters.size() << " Tpc Truth clusters" << std::endl;
    }
  }
  if (m_is_new_track) update_track();
}

void PHG4TpcTruthClusterizer::end_of_event() {
  check_g4hit(nullptr); // flush out last data if ended in truth track
  m_hitsetkey_cnt.clear();
  if (m_verbosity>2) { 
      std::cout << PHWHERE << std::endl << " :: tracks with clusters after clustering in TPC" << std::endl;
      for (auto& track : m_truthtracks->getMap()) {
        std::cout << " track("<< track.first <<") nclusters: " << track.second->getClusters().size();
        for (auto& cluster : track.second->getClusters()) std::cout << " " << (int) TrkrDefs::getLayer(cluster);
        std::cout << std::endl;
      }
  }
}

