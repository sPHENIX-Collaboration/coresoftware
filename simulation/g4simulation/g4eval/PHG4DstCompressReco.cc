#include "PHG4DstCompressReco.h"

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4Shower.h>
#include <g4main/PHG4TruthInfoContainer.h>

#include <g4detectors/PHG4CellContainer.h>

#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <trackbase_historic/PHG4ParticleSvtxMap.h>
#include <trackbase_historic/SvtxPHG4ParticleMap.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHPointerListIterator.h>
#include <phool/getClass.h>

#include <iostream>
#include <utility>

using namespace std;

PHG4DstCompressReco::PHG4DstCompressReco(const string& name)
  : SubsysReco(name)
  , _truth_info(nullptr)
  , _compress_g4hit_names()
  , _compress_g4cell_names()
  , _g4cells()
  , _g4hits()
  , _keep_g4hits()
{
}

int PHG4DstCompressReco::InitRun(PHCompositeNode* topNode)
{
  _truth_info = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!_truth_info)
  {
    cout << "PHG4DstCompressReco::InitRun(): Can't find G4TruthInfo" << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  SearchG4HitNodes(topNode);

  for (std::set<std::string>::iterator iter = _compress_g4cell_names.begin();
       iter != _compress_g4cell_names.end(); ++iter)
  {
    std::string name = *iter;

    PHG4CellContainer* g4cells = findNode::getClass<PHG4CellContainer>(topNode, name.c_str());
    if (g4cells)
    {
      _g4cells.insert(g4cells);
    }
  }

  for (std::set<std::string>::iterator iter = _compress_tower_names.begin();
       iter != _compress_tower_names.end(); ++iter)
  {
    std::string name = *iter;

    RawTowerContainer* towers = findNode::getClass<RawTowerContainer>(topNode, name.c_str());
    if (towers)
    {
      _towers.insert(towers);
    }
  }

  if (m_keepRecoTrackMatchedParticles)
  {
    _recoTruthMap = findNode::getClass<SvtxPHG4ParticleMap>(topNode, "SvtxPHG4ParticleMap");
    if (_recoTruthMap == nullptr)
    {
      cout << __PRETTY_FUNCTION__ << " Fatal error: missing SvtxPHG4ParticleMap while m_keepRecoTrackMatchedParticles is set. "
           << "Was PHG4DstCompressReco called before this module?"
           << endl;
      exit(1);
    }

    _truthRecoMap = findNode::getClass<PHG4ParticleSvtxMap>(topNode, "PHG4ParticleSvtxMap");
    if (_truthRecoMap == nullptr)
    {
      cout << __PRETTY_FUNCTION__ << " Fatal error: missing PHG4ParticleSvtxMap while m_keepRecoTrackMatchedParticles is set. "
           << "Was PHG4DstCompressReco called before this module?"
           << endl;
      exit(1);
    }
  }  //  if (m_keepRecoTrackMatchedParticles)

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4DstCompressReco::process_event(PHCompositeNode* /*topNode*/)
{
  if (_g4hits.empty() && _g4cells.empty() && _towers.empty()) return Fun4AllReturnCodes::EVENT_OK;

  //---cells--------------------------------------------------------------------

  for (std::set<PHG4CellContainer*>::iterator iter = _g4cells.begin();
       iter != _g4cells.end();
       ++iter)
  {
    PHG4CellContainer* cells = *iter;
    cells->Reset();  // DROP ALL COMPRESSED G4CELLS
  }

  //---hits---------------------------------------------------------------------

  for (std::set<PHG4HitContainer*>::iterator iter = _g4hits.begin();
       iter != _g4hits.end();
       ++iter)
  {
    PHG4HitContainer* hits = *iter;
    hits->Reset();  // DROP ALL COMPRESSED G4HITS
  }

  //---secondary particles and vertexes-----------------------------------------

  std::set<int> keep_particle_ids;
  for (std::set<PHG4HitContainer*>::iterator iter = _keep_g4hits.begin();
       iter != _keep_g4hits.end();
       ++iter)
  {
    PHG4HitContainer* hits = *iter;

    for (PHG4HitContainer::ConstIterator jter = hits->getHits().first;
         jter != hits->getHits().second;
         ++jter)
    {
      PHG4Hit* hit = jter->second;
      keep_particle_ids.insert(hit->get_trkid());
      // this will need to include all parents too in a trace back to
      // the primary, but let's start here for now
    }
  }

  //---tracker truth map, if set-----------------------------------------
  if (_truthRecoMap)
  {
    for (const auto& [particle_id, map] : *_truthRecoMap)
    {
      keep_particle_ids.insert(particle_id);
    }
  }
  if (_recoTruthMap)
  {
    for (const auto& [track_id, weighted_truth_track_map] : *_recoTruthMap)
    {
      for (const auto& [weight, particle_set] : weighted_truth_track_map)
      {
        if (weight > 0)
        {
          for (const auto& particle_id : particle_set)
          {
            keep_particle_ids.insert(particle_id);
          }
        }
      }
    }  //    for (const auto& [track_id, weighted_truth_track_map] : *_recoTruthMap)

  }  //  if (_recoTruthMap)

  std::set<int> keep_vertex_ids;
  PHG4TruthInfoContainer::Range range = _truth_info->GetSecondaryParticleRange();
  for (PHG4TruthInfoContainer::Iterator iter = range.first;
       iter != range.second;)
  {
    int id = iter->first;
    PHG4Particle* particle = iter->second;

    if (keep_particle_ids.find(id) != keep_particle_ids.end())
    {
      ++iter;
      keep_vertex_ids.insert(particle->get_vtx_id());
      continue;
    }
    else
    {
      _truth_info->delete_particle(iter++);  // DROP PARTICLES NOT ASSOCIATED TO A PRESERVED HIT
    }
  }

  PHG4TruthInfoContainer::VtxRange vrange = _truth_info->GetSecondaryVtxRange();
  for (PHG4TruthInfoContainer::VtxIterator iter = vrange.first;
       iter != vrange.second;)
  {
    int id = iter->first;

    if (keep_vertex_ids.find(id) != keep_vertex_ids.end())
    {
      ++iter;
      continue;
    }
    else
    {
      _truth_info->delete_vtx(iter++);  // DROP VERTEXES NOT ASSOCIATED TO A PRESERVED HIT
    }
  }

  //---shower entries-----------------------------------------------------------

  PHG4TruthInfoContainer::ShowerRange srange = _truth_info->GetShowerRange();
  for (PHG4TruthInfoContainer::ShowerIterator iter = srange.first;
       iter != srange.second;
       ++iter)
  {
    PHG4Shower* shower = iter->second;

    shower->clear_g4particle_id();
    shower->clear_g4vertex_id();
    shower->clear_g4hit_id();
  }

  //---tower cell entries-------------------------------------------------------
  for (std::set<RawTowerContainer*>::iterator iter = _towers.begin();
       iter != _towers.end();
       ++iter)
  {
    RawTowerContainer* towers = *iter;

    // loop over all the towers
    for (RawTowerContainer::Iterator jter = towers->getTowers().first;
         jter != towers->getTowers().second;
         ++jter)
    {
      RawTower* tower = jter->second;
      tower->clear_g4cells();
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void PHG4DstCompressReco::SearchG4HitNodes(PHCompositeNode* top)
{
  // fill a lookup map between the g4hit container ids and the containers
  // themselves
  // without knowing what the container names are in advance, only that they
  // begin G4HIT_*

  // separate the names into those in the compression list and those not in the
  // compression list

  PHNodeIterator nodeiter(top);
  PHPointerListIterator<PHNode> iter(nodeiter.ls());
  PHNode* thisNode;
  while ((thisNode = iter()))
  {
    if (thisNode->getType() == "PHCompositeNode")
    {
      SearchG4HitNodes(static_cast<PHCompositeNode*>(thisNode));
    }
    else if (thisNode->getType() == "PHIODataNode")
    {
      if (thisNode->getName().find("G4HIT_") == 0)
      {
        PHIODataNode<PHG4HitContainer>* DNode =
            static_cast<PHIODataNode<PHG4HitContainer>*>(thisNode);
        if (DNode)
        {
          PHG4HitContainer* object =
              dynamic_cast<PHG4HitContainer*>(DNode->getData());
          if (object)
          {
            if (_compress_g4hit_names.find(thisNode->getName()) !=
                _compress_g4hit_names.end())
            {
              _g4hits.insert(object);
            }
            else
            {
              _keep_g4hits.insert(object);
            }
          }
        }
      }
    }
  }
}
