#include "PHTruthTrackSeeding.h"

#include "AssocInfoContainer.h"

#include <trackbase_historic/SvtxTrack.h>     // for SvtxTrack, SvtxTra...
#include <trackbase_historic/SvtxTrackMap.h>  // for SvtxTrackMap, Svtx...
#include <trackbase_historic/SvtxTrack_FastSim.h>

#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterHitAssoc.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHitTruthAssoc.h>

#include <g4eval/SvtxClusterEval.h>
#include <g4eval/SvtxEvalStack.h>
#include <g4eval/SvtxEvaluator.h>
#include <g4eval/SvtxTrackEval.h>

#include <g4main/PHG4Hit.h>  // for PHG4Hit
#include <g4main/PHG4Particle.h>  // for PHG4Particle
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4HitDefs.h>  // for keytype
#include <g4main/PHG4TruthInfoContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/phool.h>

#include <cstdlib>   // for exit
#include <iostream>  // for operator<<, endl
#include <map>       // for multimap, map<>::c...
#include <memory>
#include <utility>  // for pair
#include <cassert>
#include <set>

#define LogDebug(exp) std::cout << "DEBUG: " << __FILE__ << ": " << __LINE__ << ": " << exp << std::endl
#define LogError(exp) std::cout << "ERROR: " << __FILE__ << ": " << __LINE__ << ": " << exp << std::endl
#define LogWarning(exp) std::cout << "WARNING: " << __FILE__ << ": " << __LINE__ << ": " << exp << std::endl

class PHCompositeNode;

using namespace std;

PHTruthTrackSeeding::PHTruthTrackSeeding(const std::string& name)
  : PHTrackSeeding(name)
{}

int PHTruthTrackSeeding::Setup(PHCompositeNode* topNode)
{
  cout << "Enter PHTruthTrackSeeding:: Setup" << endl;

  int ret = PHTrackSeeding::Setup(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHTruthTrackSeeding::Process(PHCompositeNode* topNode)
{
  typedef std::map<int, std::set<TrkrCluster*> > TrkClustersMap;
  TrkClustersMap m_trackID_clusters;

  // loop over all clusters
  TrkrClusterContainer::ConstRange clusrange = _cluster_map->getClusters();
  for (TrkrClusterContainer::ConstIterator clusiter = clusrange.first; clusiter != clusrange.second; ++clusiter)
  {
    TrkrCluster* cluster = clusiter->second;
    TrkrDefs::cluskey cluskey = clusiter->first;
    unsigned int trkrid = TrkrDefs::getTrkrId(cluskey);
    unsigned int layer = TrkrDefs::getLayer(cluskey);
    if(layer<_min_layer) continue;
    if(layer>=_max_layer) continue;

    if (Verbosity() >= 3)
    {
      cout <<__PRETTY_FUNCTION__<<" process cluster ";
      cluster->identify();
    }

    // get the hits for this cluster
    TrkrClusterHitAssoc::ConstRange hitrange = clusterhitassoc->getHits(cluskey);  // returns range of pairs {cluster key, hit key} for this cluskey
    for (TrkrClusterHitAssoc::ConstIterator clushititer = hitrange.first; clushititer != hitrange.second; ++clushititer)
    {
      TrkrDefs::hitkey hitkey = clushititer->second;
      // TrkrHitTruthAssoc uses a map with (hitsetkey, std::pair(hitkey, g4hitkey)) - get the hitsetkey from the cluskey
      TrkrDefs::hitsetkey hitsetkey = TrkrDefs::getHitSetKeyFromClusKey(cluskey);

      // get all of the g4hits for this hitkey
      std::multimap<TrkrDefs::hitsetkey, std::pair<TrkrDefs::hitkey, PHG4HitDefs::keytype> > temp_map;
      hittruthassoc->getG4Hits(hitsetkey, hitkey, temp_map);  // returns pairs (hitsetkey, std::pair(hitkey, g4hitkey)) for this hitkey only
      for (std::multimap<TrkrDefs::hitsetkey, std::pair<TrkrDefs::hitkey, PHG4HitDefs::keytype> >::iterator htiter = temp_map.begin(); htiter != temp_map.end(); ++htiter)
      {
        // extract the g4 hit key here and add the hits to the set
        PHG4HitDefs::keytype g4hitkey = htiter->second.second;
        PHG4Hit* phg4hit = nullptr;
        switch( trkrid )
        {
          case TrkrDefs::mvtxId:
          if (phg4hits_mvtx) phg4hit = phg4hits_mvtx->findHit( g4hitkey );
          break;

          case TrkrDefs::inttId:
          if (phg4hits_intt) phg4hit = phg4hits_intt->findHit( g4hitkey );
          break;

          case TrkrDefs::tpcId:
          if (phg4hits_tpc) phg4hit = phg4hits_tpc->findHit( g4hitkey );
          break;

          case TrkrDefs::micromegasId:
          if (phg4hits_micromegas) phg4hit = phg4hits_micromegas->findHit( g4hitkey );
          break;
        }
   
        if( !phg4hit )
        {
          std::cout<<PHWHERE<<" unable to find g4hit from key " << g4hitkey << std::endl;
          continue;
        }
     
        int particle_id = phg4hit->get_trkid();

        // monentum cut-off
        if (_min_momentum>0)
        {
          PHG4Particle* particle = _g4truth_container->GetParticle(particle_id);
          if (!particle)
          {
            cout <<__PRETTY_FUNCTION__<<" - validity check failed: missing truth particle with ID of "<<particle_id<<". Exiting..."<<endl;
            exit(1);
          }
          const double monentum2 =
              particle->get_px() * particle->get_px()
              +
              particle->get_py() * particle->get_py()
              +
              particle->get_pz() * particle->get_pz()
              ;

          if (Verbosity() >= 10)
          {
            cout <<__PRETTY_FUNCTION__<<" check momentum for particle"<<particle_id<<" -> cluster "<<cluskey
                <<" = "<<sqrt(monentum2)<<endl;;
            particle->identify();
          }

          if (monentum2 < _min_momentum * _min_momentum)
          {
            if (Verbosity() >= 3)
            {
              cout <<__PRETTY_FUNCTION__<<" ignore low momentum particle"<<particle_id<<" -> cluster "<<cluskey<<endl;;
              particle->identify();
            }
            continue;
          }
        }


        TrkClustersMap::iterator it = m_trackID_clusters.find(particle_id);

        if (it != m_trackID_clusters.end())
        {
          it->second.insert(cluster);
          if (Verbosity() >= 3)
          {
            cout <<__PRETTY_FUNCTION__<<" append particle"<<particle_id<<" -> cluster "<<cluskey<<endl;;
            cluster->identify();
          }
        }
        else
        {
          std::set<TrkrCluster*> clusters;
          clusters.insert(cluster);
          m_trackID_clusters.insert(std::pair<int, std::set<TrkrCluster*> >(particle_id, clusters));


          if (Verbosity() >= 3)
          {
            cout <<__PRETTY_FUNCTION__<<" new particle"<<particle_id<<" -> cluster "<<cluskey<<endl;;
            cluster->identify();
          }

        }
      }  // loop over g4hits associated with hit
    }    // loop over hits associated with cluster
  }      // loop over clusters

  //==================================

  if (Verbosity() >= 2)
  {
    cout <<__PRETTY_FUNCTION__
        <<" _track_map->size = "<<_track_map->size()<<endl;
  }

  // Build track
  for (TrkClustersMap::const_iterator trk_clusters_itr = m_trackID_clusters.begin();
       trk_clusters_itr != m_trackID_clusters.end(); ++trk_clusters_itr)
  {
    if (trk_clusters_itr->second.size() <  _min_clusters_per_track)
      continue;

    // check number of layers also pass the _min_clusters_per_track cut to avoid tight loopers
    set<uint8_t> layers;
    for (TrkrCluster* cluster : trk_clusters_itr->second)
    {
      assert(cluster);
      const uint8_t layer = TrkrDefs::getLayer(cluster->getClusKey());
      layers.insert(layer);
    }
    if (Verbosity()>2)
    {
      cout <<__PRETTY_FUNCTION__<<" particle "<<trk_clusters_itr->first<<" -> "
          <<trk_clusters_itr->second.size()<<" clusters covering "<<layers.size()<<" layers."
          <<" Layer/clusters cuts are > "<<_min_clusters_per_track
          <<endl;
    }

    if (layers.size() >=  _min_clusters_per_track)
    {

      std::unique_ptr<SvtxTrack_FastSim> svtx_track(new SvtxTrack_FastSim());

      svtx_track->set_id(_track_map->size());
      svtx_track->set_truth_track_id(trk_clusters_itr->first);

      // dummy values, set px to make it through the minimum pT cut
      svtx_track->set_px(10.);
      svtx_track->set_py(0.);
      svtx_track->set_pz(0.);
      for (TrkrCluster* cluster : trk_clusters_itr->second)
      {
        svtx_track->insert_cluster_key(cluster->getClusKey());
        _assoc_container->SetClusterTrackAssoc(cluster->getClusKey(), svtx_track->get_id());
      }

      _track_map->insert(svtx_track.get());

      if (Verbosity() >= 2)
      {
        cout <<__PRETTY_FUNCTION__<<" particle "<<trk_clusters_itr->first<<" -> "
            <<trk_clusters_itr->second.size()<<" clusters"
            <<" _track_map->size = "<< (_track_map->size()) <<": ";
        _track_map->identify();
      }
    }
  }

  if (Verbosity() >= 2)
  {
    cout << "Loop over SvtxTrackMap entries " << endl;
    for (SvtxTrackMap::Iter iter = _track_map->begin();
         iter != _track_map->end(); ++iter)
    {
      SvtxTrack* svtx_track = iter->second;

      //svtx_track->identify();
      //continue;

      cout << "Track ID: " << svtx_track->get_id() << ", Dummy Track pT: "
           << svtx_track->get_pt() << ", Truth Track/Particle ID: "
           << svtx_track->get_truth_track_id() << endl;

      //Print associated clusters;
      for (SvtxTrack::ConstClusterKeyIter iter =
               svtx_track->begin_cluster_keys();
           iter != svtx_track->end_cluster_keys(); ++iter)
      {
        TrkrDefs::cluskey cluster_key = *iter;
        TrkrCluster* cluster = _cluster_map->findCluster(cluster_key);
        float radius = sqrt(
            cluster->getX() * cluster->getX() + cluster->getY() * cluster->getY());
        cout << "       cluster ID: "
             << cluster->getClusKey() << ", cluster radius: " << radius
             << endl;
      }
    }
  }
  //==================================

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHTruthTrackSeeding::GetNodes(PHCompositeNode* topNode)
{
  _g4truth_container = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!_g4truth_container)
  {
    cerr << PHWHERE << " ERROR: Can't find node G4TruthInfo" << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  clusterhitassoc = findNode::getClass<TrkrClusterHitAssoc>(topNode, "TRKR_CLUSTERHITASSOC");
  if (!clusterhitassoc)
  {
    cout << PHWHERE << "Failed to find TRKR_CLUSTERHITASSOC node, quit!" << endl;
    exit(1);
  }

  hittruthassoc = findNode::getClass<TrkrHitTruthAssoc>(topNode, "TRKR_HITTRUTHASSOC");
  if (!hittruthassoc)
  {
    cout << PHWHERE << "Failed to find TRKR_HITTRUTHASSOC node, quit!" << endl;
    exit(1);
  }

  using nodePair = std::pair<std::string, PHG4HitContainer*&>;
  std::initializer_list<nodePair> nodes =
  {
    { "G4HIT_TPC", phg4hits_tpc },
    { "G4HIT_INTT", phg4hits_intt },
    { "G4HIT_MVTX", phg4hits_mvtx },
    { "G4HIT_MICROMEGAS", phg4hits_micromegas }
  };
  
  for( auto&& node: nodes )
  {
    if( !( node.second = findNode::getClass<PHG4HitContainer>( topNode, node.first ) ) )
    { std::cerr << PHWHERE << " PHG4HitContainer " << node.first << " not found" << std::endl; }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHTruthTrackSeeding::End()
{
  return 0;
}
