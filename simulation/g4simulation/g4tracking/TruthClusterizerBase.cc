#include "TruthClusterizerBase.h"

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4tracking/TrkrTruthTrack.h>
#include <g4tracking/TrkrTruthTrackContainer.h>
#include <g4tracking/TrkrTruthTrackContainerv1.h>
#include <intt/CylinderGeomIntt.h>

/* #include <trackbase/InttDefs.h> */
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrClusterContainerv4.h>
#include <trackbase/TrkrClusterv4.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainerv1.h>
#include <trackbase/TrkrHitv2.h>  // for TrkrHit

#include <g4detectors/PHG4CylinderGeom.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <g4mvtx/PHG4MvtxDigitizer.h>
#include <mvtx/MvtxHitPruner.h>
#include <mvtx/MvtxClusterizer.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>                         // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>

#include <TMatrixFfwd.h>                            // for TMatrixF
#include <TMatrixT.h>                               // for TMatrixT, operator*
#include <TMatrixTUtils.h>                          // for TMatrixTRow


#include <array>
#include <cmath>
#include <iostream>
#include <set>
#include <vector> 

TruthClusterizerBase::TruthClusterizerBase ( )
    : m_hits      { new TrkrHitSetContainerv1 }
{ }

void TruthClusterizerBase::init_clusterizer_base( PHCompositeNode*& _topNode, int _verbosity ) {
  m_topNode = _topNode;
  m_verbosity = _verbosity;
  PHNodeIterator iter(m_topNode);
  auto dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  assert(dstNode);
  PHNodeIterator dstiter(dstNode);
  auto DetNode = dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", "TRKR"));
  if (!DetNode)
  {
    DetNode = new PHCompositeNode("TRKR");
    dstNode->addNode(DetNode);
  }

  m_truthtracks = findNode::getClass<TrkrTruthTrackContainer>(m_topNode, "TRKR_TRUTHTRACKCONTAINER");
  if (!m_truthtracks)
  {
    PHNodeIterator dstiter(dstNode);
    m_truthtracks = new TrkrTruthTrackContainerv1();
    auto newNode = new PHIODataNode<PHObject>(m_truthtracks, "TRKR_TRUTHTRACKCONTAINER", "PHObject");
    DetNode->addNode(newNode);
  }

  m_clusters = findNode::getClass<TrkrClusterContainer>(m_topNode, "TRKR_TRUTHCLUSTERCONTAINER");
  if (!m_clusters)
  {
    m_clusters = new TrkrClusterContainerv4;
    auto newNode = new PHIODataNode<PHObject>(m_clusters, "TRKR_TRUTHCLUSTERCONTAINER", "PHObject");
    DetNode->addNode(newNode);
  }

  m_truthinfo = findNode::getClass<PHG4TruthInfoContainer>(m_topNode, "G4TruthInfo");
  if (!m_truthinfo)
  {
    std::cout << PHWHERE << " PHG4TruthInfoContainer node not found on node tree" << std::endl;
    assert(m_truthinfo);
  }
}      
      
TruthClusterizerBase::~TruthClusterizerBase() {
  delete m_hits;
}

void TruthClusterizerBase::check_g4hit_status(PHG4Hit* hit) {
  int new_trkid = (hit==nullptr) ? -1 : hit->get_trkid();
  m_is_new_track = (new_trkid != m_trkid);
  if (m_verbosity>5) std::cout << PHWHERE << std::endl << " -> Checking status of PHG4Hit. Track id("<<new_trkid<<")" << std::endl;
  if (!m_is_new_track) return;

  m_trkid = new_trkid; // although m_current_track won't catch up until update_track is called
  m_was_emb = m_is_emb;
  m_is_emb = m_truthinfo->isEmbeded(m_trkid);
}

// to call if m_was_emb=true and after clustering
void TruthClusterizerBase::transfer_clusters(TrkrClusterContainer* pass_clusters) {
  m_hits->Reset(); // clear out the old  hits
  for (auto hitsetkey : pass_clusters->getHitSetKeys()) {
    m_hitsetkey_cnt.try_emplace(hitsetkey,0);
    unsigned int& cnt = m_hitsetkey_cnt[hitsetkey];
    auto range = pass_clusters->getClusters(hitsetkey);
    for (auto cluster = range.first; cluster != range.second; ++cluster) {
      auto ckey = TrkrDefs::genClusKey(hitsetkey, cnt);
      m_clusters->addClusterSpecifyKey(ckey, cluster->second);
      m_current_track->addCluster(ckey);
      ++cnt;
    }
  }
  m_was_emb = false;
}

// if m_is_new_track
void TruthClusterizerBase::update_track() {
  m_current_track = m_is_emb ? m_truthtracks->getTruthTrack(m_trkid, m_truthinfo) : nullptr;
}

void TruthClusterizerBase::addhitset(
    TrkrDefs::hitsetkey hitsetkey, 
    TrkrDefs::hitkey hitkey, 
    float neffelectrons) 
{
  if (!m_is_emb) return;
  TrkrHitSetContainer::Iterator hitsetit = m_hits->findOrAddHitSet(hitsetkey);
  // See if this hit already exists
  TrkrHit *hit = nullptr;
  hit = hitsetit->second->getHit(hitkey);
  if (!hit)
  {
    // create a new one
    hit = new TrkrHitv2();
    hitsetit->second->addHitSpecificKey(hitkey, hit);
  }
  // Either way, add the energy to it  -- adc values will be added at digitization
  hit->addEnergy(neffelectrons);
}

void TruthClusterizerBase::print_clusters(int nclusprint) {
  std::cout << PHWHERE << ": content of clusters " << std::endl;
  auto& tmap = m_truthtracks->getMap();
  std::cout << " Number of tracks: " << tmap.size() << std::endl;
  for (auto& _pair : tmap) {
    auto& track = _pair.second;

    printf("id(%2i) phi:eta:pt(", (int)track->getTrackid());
    std::cout << "phi:eta:pt(";
    printf("%5.2f:%5.2f:%5.2f", track->getPhi(), track->getPseudoRapidity(), track->getPt());
      /* Form("%5.2:%5.2:%5.2", track->getPhi(), track->getPseudoRapidity(), track->getPt()) */
      //<<track->getPhi()<<":"<<track->getPseudoRapidity()<<":"<<track->getPt() 
    std::cout << ") nclusters(" << track->getClusters().size() <<") ";
    if (m_verbosity <= 10) { std::cout << std::endl; }
    else {
      int nclus = 0;
      for (auto cluskey : track->getClusters()) {
        std::cout << " " 
          << ((int) TrkrDefs::getHitSetKeyFromClusKey(cluskey)) <<":index(" <<
          ((int)  TrkrDefs::getClusIndex(cluskey)) << ")";
        ++nclus;
        if (nclusprint > 0 && nclus >= nclusprint) {
          std::cout << " ... "; 
          break;
        }
      }
    }
  }
  std::cout << PHWHERE << " ----- end of clusters " << std::endl;
}

