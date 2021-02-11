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
#include <trackbase_historic/SvtxVertexMap.h>
#include <trackbase_historic/SvtxVertex.h>

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
  , _clustereval(nullptr)
{}

int PHTruthTrackSeeding::Setup(PHCompositeNode* topNode)
{
  cout << "Enter PHTruthTrackSeeding:: Setup" << endl;

  int ret = PHTrackSeeding::Setup(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;
  _clustereval = new  SvtxClusterEval(topNode);
  _clustereval->do_caching(true);
  // _clustereval.set_strict(strict);
  //  _clustereval.set_verbosity(verbosity);
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHTruthTrackSeeding::Process(PHCompositeNode* topNode)
{
  _clustereval->next_event(topNode);
  _track_map->Clear();

  typedef std::map<int, std::set<TrkrCluster*> > TrkClustersMap;
  TrkClustersMap m_trackID_clusters;

  vector<TrkrDefs::cluskey> ClusterKeyList; 

  PHG4TruthInfoContainer::ConstRange range = _g4truth_container->GetPrimaryParticleRange();
  //  Float_t gntracks = (Float_t) truthinfo->GetNumPrimaryVertexParticles();
  for (PHG4TruthInfoContainer::ConstIterator iter = range.first;
       iter != range.second;
       ++iter){
    ClusterKeyList.clear();
    PHG4Particle* g4particle = iter->second;

    if (g4particle==NULL){
      cout <<__PRETTY_FUNCTION__<<" - validity check failed: missing truth particle" <<endl;
      exit(1);
    }

    float gtrackID = g4particle->get_track_id();

    // monentum cut-off
    if (_min_momentum>0){
      const double monentum2 =
	g4particle->get_px() * g4particle->get_px()+
	g4particle->get_py() * g4particle->get_py()+
	g4particle->get_pz() * g4particle->get_pz();
      
      if (monentum2 < _min_momentum * _min_momentum){
	if (Verbosity() >= 3){
	  cout <<__PRETTY_FUNCTION__<<" ignore low momentum particle "<< gtrackID <<endl;
	  g4particle->identify();
	}
	continue;
      }
    }

    for(Int_t layer = _min_layer;layer < _max_layer;layer++){
      TrkrDefs::cluskey cluskey = _clustereval->best_cluster_by_nhit(gtrackID, layer);
      if(cluskey!=0)
	ClusterKeyList.push_back(cluskey);
    }
    if(ClusterKeyList.size()< _min_clusters_per_track)
      continue;

    std::unique_ptr<SvtxTrack_FastSim> svtx_track(new SvtxTrack_FastSim());
    svtx_track->set_id(_track_map->size());
    svtx_track->set_truth_track_id(gtrackID);
   
    // Smear the truth values out by 5% so that the seed momentum and
    // position aren't completely exact
    
    double random = 0;//((double) rand() / (RAND_MAX)) * 0.05;
    // make it negative sometimes
    if(rand() % 2)
      random *= -1;
    svtx_track->set_px(g4particle->get_px() * (1 + random));
    svtx_track->set_py(g4particle->get_py() * (1 + random));
    svtx_track->set_pz(g4particle->get_pz() * (1 + random));
      
    for (TrkrDefs::cluskey cluskey : ClusterKeyList){
      svtx_track->insert_cluster_key(cluskey);
      _assoc_container->SetClusterTrackAssoc(cluskey,svtx_track->get_id());
    }
    svtx_track->set_ndf(ClusterKeyList.size()*3-5);
    svtx_track->set_chisq(1.5*ClusterKeyList.size()*3-5);
    _track_map->insert(svtx_track.get());
  }
  
  if (Verbosity() >= 5)
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
      cout << " nhits: " << svtx_track->size_cluster_keys()<< endl;
      //Print associated clusters;
      for (SvtxTrack::ConstClusterKeyIter iter_clus =
               svtx_track->begin_cluster_keys();
           iter_clus != svtx_track->end_cluster_keys(); ++iter_clus)
      {
        TrkrDefs::cluskey cluster_key = *iter_clus;
	cout << "Key: "  << cluster_key<< endl;
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
