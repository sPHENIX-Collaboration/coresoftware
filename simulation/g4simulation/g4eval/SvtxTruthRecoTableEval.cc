
#include "SvtxTruthRecoTableEval.h"
#include "SvtxEvalStack.h"
#include "SvtxTrackEval.h"

#include "SvtxClusterEval.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>

#include <trackbase_historic/PHG4ParticleSvtxMap_v1.h>
#include <trackbase_historic/SvtxPHG4ParticleMap_v1.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrack_FastSim.h>
#include <trackbase_historic/SvtxTrackMap.h>

#include <cassert>
#include <cstddef>
#include <iostream>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

//____________________________________________________________________________..
SvtxTruthRecoTableEval::SvtxTruthRecoTableEval(const std::string &name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
SvtxTruthRecoTableEval::~SvtxTruthRecoTableEval() = default;

//____________________________________________________________________________..
int SvtxTruthRecoTableEval::Init(PHCompositeNode * /*unused*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int SvtxTruthRecoTableEval::InitRun(PHCompositeNode *topNode)
{
  if (createNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
  {
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int SvtxTruthRecoTableEval::process_event(PHCompositeNode *topNode)
{
  const int verbosity = Verbosity();

  if (!m_svtxevalstack)
  {
    m_svtxevalstack = std::make_unique<SvtxEvalStack>(topNode);
    m_svtxevalstack->set_strict(false);
    m_svtxevalstack->set_verbosity(verbosity);
    m_svtxevalstack->set_use_initial_vertex(true);
    m_svtxevalstack->set_use_genfit_vertex(false);
    m_svtxevalstack->next_event(topNode);
  }
  else
  {
    m_svtxevalstack->next_event(topNode);
  }

  SvtxTrackEval *trackeval = m_svtxevalstack->get_track_eval();
  assert(trackeval);
  trackeval->set_verbosity(verbosity);

  if (verbosity > 1)
  {
    std::cout << "Fill truth/reco maps " << std::endl;
  }
  fillTruthRecoMaps(topNode, trackeval, verbosity);

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int SvtxTruthRecoTableEval::ResetEvent(PHCompositeNode * /*unused*/)
{
  if (Verbosity() > 0)
  {
    std::cout << "Truth track map " << std::endl;
    m_truthMap->identify();
    std::cout << std::endl
              << "Reco track map " << std::endl;
    m_recoMap->identify();
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int SvtxTruthRecoTableEval::End(PHCompositeNode * /*unused*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

void SvtxTruthRecoTableEval::fillTruthRecoMaps(PHCompositeNode *topNode, SvtxTrackEval *trackeval, const int verbosity)
{
  PHG4TruthInfoContainer *truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  assert(truthinfo);

  SvtxTrackMap *trackMap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  assert(trackMap);

  PHG4TruthInfoContainer::ConstRange range = truthinfo->GetParticleRange();
  if (m_scanForPrimaries)
  {
    range = truthinfo->GetPrimaryParticleRange();
  }

  std::vector<int> selectedTruthIds;
  std::unordered_set<int> selectedTruthIdSet;
  const double minMomentumTruthMap2 = m_minMomentumTruthMap * m_minMomentumTruthMap;

  for (auto iter = range.first; iter != range.second; ++iter)
  {
    PHG4Particle *g4particle = iter->second;

    const double px = g4particle->get_px();
    const double py = g4particle->get_py();
    const double pz = g4particle->get_pz();
    const double momentum2 = px * px + py * py + pz * pz;

    // only record particle above minimal momentum (square) requirement.
    // doing this saves us a slow sqrt operation to calculate the momentum itself
    if (momentum2 < minMomentumTruthMap2)
    {
      continue;
    }

    const int gtrackID = g4particle->get_track_id();
    selectedTruthIds.push_back(gtrackID);
    selectedTruthIdSet.insert(gtrackID);
  }

  SvtxClusterEval *clustereval = trackeval->get_cluster_eval();
  std::map<int, PHG4ParticleSvtxMap::WeightedRecoTrackMap> truthMaps;

  for (const auto &[key, track] : *trackMap)
  {
    TrackSeed *siliconSeed = track->get_silicon_seed();
    TrackSeed *tpcSeed = track->get_tpc_seed();

    std::size_t nclusterKeys = 0;
    if (siliconSeed)
    {
      nclusterKeys += siliconSeed->size_cluster_keys();
    }
    if (tpcSeed)
    {
      nclusterKeys += tpcSeed->size_cluster_keys();
    }

    std::unordered_map<int, unsigned int> nclustersByTruthId;
    nclustersByTruthId.reserve(nclusterKeys);

    const auto add_cluster_contributions = [&](TrackSeed *seed)
    {
      if (!seed)
      {
        return;
      }

      for (auto clusterIter = seed->begin_cluster_keys();
           clusterIter != seed->end_cluster_keys();
           ++clusterIter)
      {
        const std::set<PHG4Particle *> particles = clustereval->all_truth_particles(*clusterIter);
        for (PHG4Particle *g4particle : particles)
        {
          ++nclustersByTruthId[g4particle->get_track_id()];
        }
      }
    };

    // Match SvtxTrackEval::get_track_ckeys ordering.
    add_cluster_contributions(siliconSeed);
    add_cluster_contributions(tpcSeed);

    SvtxPHG4ParticleMap::WeightedTruthTrackMap truthmap;
    SvtxTrack_FastSim *fastsim_track = dynamic_cast<SvtxTrack_FastSim *>(track);

    const unsigned int trackID = track->get_id();
    for (const auto &[gtrackID, nclusters] : nclustersByTruthId)
    {
      const float clusCont = static_cast<float>(nclusters);
      if (selectedTruthIdSet.contains(gtrackID))
      {
        truthMaps[gtrackID][clusCont].insert(trackID);
      }
      if (!fastsim_track)
      {
        truthmap[clusCont].insert(gtrackID);
      }
    }

    if (fastsim_track)
    {
      // Preserve SvtxTrackEval::all_truth_particles fast-sim special case for reco->truth maps only.
      PHG4Particle *g4particle = truthinfo->GetParticle(fastsim_track->get_truth_track_id());
      const float clusCont = trackeval->get_nclusters_contribution(track, g4particle);
      truthmap[clusCont].insert(g4particle->get_track_id());
    }

    if (verbosity > 1)
    {
      std::cout << " Inserting track id " << key << " with truth map size " << truthmap.size() << std::endl;
    }
    m_recoMap->insert(key, std::move(truthmap));
  }

  for (const int gtrackID : selectedTruthIds)
  {
    auto truthMapIter = truthMaps.find(gtrackID);
    if (truthMapIter == truthMaps.end() || truthMapIter->second.empty())
    {
      continue;
    }

    if (verbosity > 1)
    {
      std::cout << " Inserting gtrack id " << gtrackID << " with map size " << truthMapIter->second.size() << std::endl;
    }

    m_truthMap->insert(gtrackID, std::move(truthMapIter->second));
  }

  m_truthMap->setProcessed(true);
  m_recoMap->setProcessed(true);
}

int SvtxTruthRecoTableEval::createNodes(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));

  if (!dstNode)
  {
    std::cerr << "DST node is missing, quitting" << std::endl;
    throw std::runtime_error("Failed to find DST node in SvtxTruthRecoTableEval::createNodes");
  }

  PHCompositeNode *svtxNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "SVTX"));

  if (!svtxNode)
  {
    svtxNode = new PHCompositeNode("SVTX");
    dstNode->addNode(svtxNode);
  }

  m_truthMap = findNode::getClass<PHG4ParticleSvtxMap>(topNode, "PHG4ParticleSvtxMap");
  if (!m_truthMap)
  {
    m_truthMap = new PHG4ParticleSvtxMap_v1;
    PHIODataNode<PHObject> *truthNode =
        new PHIODataNode<PHObject>(m_truthMap, "PHG4ParticleSvtxMap", "PHObject");
    svtxNode->addNode(truthNode);
  }

  m_recoMap = findNode::getClass<SvtxPHG4ParticleMap>(topNode, "SvtxPHG4ParticleMap");
  if (!m_recoMap)
  {
    m_recoMap = new SvtxPHG4ParticleMap_v1;
    PHIODataNode<PHObject> *recoNode =
        new PHIODataNode<PHObject>(m_recoMap, "SvtxPHG4ParticleMap", "PHObject");
    svtxNode->addNode(recoNode);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
