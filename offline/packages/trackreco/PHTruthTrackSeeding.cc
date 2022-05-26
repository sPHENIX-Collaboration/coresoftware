#include "PHTruthTrackSeeding.h"

#include <trackbase_historic/TrackSeed.h>     
#include <trackbase_historic/TrackSeedContainer.h> 
#include <trackbase_historic/TrackSeed_FastSim_v1.h>

#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>


#include <trackbase/TrkrHitTruthAssoc.h>
#include <trackbase_historic/SvtxVertexMap.h>
#include <trackbase_historic/SvtxVertex.h>
#include <trackbase_historic/ActsTransformations.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <g4eval/SvtxClusterEval.h>
#include <g4eval/SvtxEvalStack.h>
#include <g4eval/SvtxEvaluator.h>
#include <g4eval/SvtxTrackEval.h>

#include <g4main/PHG4Hit.h>  // for PHG4Hit
#include <g4main/PHG4Particle.h>  // for PHG4Particle
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4HitDefs.h>  // for keytype
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>                 
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>  
#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/PHRandomSeed.h>

#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterCrossingAssoc.h> 
#include <trackbase/TrkrHitTruthAssoc.h>

#include <trackbase_historic/ActsTransformations.h>
#include <trackbase_historic/TrackSeed_v1.h>
#include <trackbase_historic/SvtxTrackSeed_v1.h>
#include <trackbase_historic/TrackSeedContainer_v1.h>
#include <trackbase_historic/SvtxVertexMap.h>
#include <trackbase_historic/SvtxVertex.h>

#include <TDatabasePDG.h>
#include <TParticlePDG.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>                            // for gsl_rng_alloc

#include <cassert>
#include <cstdlib>   // for exit
#include <iostream>  // for operator<<, std::endl 
#include <map>       // for multimap, map<>::c...
#include <memory>
#include <set>
#include <utility>  // for pair

class PHCompositeNode;

namespace
{ template< class T> inline constexpr T square( const T& x ) { return x*x; } }

PHTruthTrackSeeding::PHTruthTrackSeeding(const std::string& name)
  : PHTrackSeeding(name)
  , _clustereval(nullptr)
{
  // initialize random generator
  const uint seed = PHRandomSeed();
  m_rng.reset( gsl_rng_alloc(gsl_rng_mt19937) );
  gsl_rng_set( m_rng.get(), seed );
}

int PHTruthTrackSeeding::Setup(PHCompositeNode* topNode)
{
  std::cout << "Enter PHTruthTrackSeeding:: Setup" << std::endl ;

  int ret = PHTrackSeeding::Setup(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  if(_track_map_name.find("Svtx") != std::string::npos)
    { ret = CreateNodes(topNode); }
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  _clustereval = new  SvtxClusterEval(topNode);
  _clustereval->do_caching(true);
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHTruthTrackSeeding::Process(PHCompositeNode* topNode)
{
  _clustereval->next_event(topNode);
  _track_map->Clear();

  using TrkClustersMap = std::map<int, std::set<TrkrCluster*> >;
  TrkClustersMap m_trackID_clusters;

  std::vector<TrkrDefs::cluskey> ClusterKeyList; 

  PHG4TruthInfoContainer::ConstRange range = m_g4truth_container->GetPrimaryParticleRange();
  for (PHG4TruthInfoContainer::ConstIterator iter = range.first;
       iter != range.second;
       ++iter){
    ClusterKeyList.clear();
    PHG4Particle* g4particle = iter->second;

    if (g4particle==NULL){
      std::cout <<__PRETTY_FUNCTION__<<" - validity check failed: missing truth particle" << std::endl;
      exit(1);
    }

    const float gtrackID = g4particle->get_track_id();
    
    // monentum cut-off
    if (_min_momentum>0){
      const double monentum2 =
	g4particle->get_px() * g4particle->get_px()+
	g4particle->get_py() * g4particle->get_py()+
	g4particle->get_pz() * g4particle->get_pz();
      
      if (monentum2 < _min_momentum * _min_momentum){
	if (Verbosity() >= 3){
	  std::cout <<__PRETTY_FUNCTION__<<" ignore low momentum particle "<< gtrackID << std::endl;
	  g4particle->identify();
	}
	continue;
      }
    }

    for(unsigned int layer = _min_layer;layer < _max_layer;layer++){
      TrkrDefs::cluskey cluskey = _clustereval->best_cluster_by_nhit(gtrackID, layer);
      if(cluskey!=0) 
	ClusterKeyList.push_back(cluskey);
    }
    if(ClusterKeyList.size()< _min_clusters_per_track)
      continue;

    if(_track_map_name.find("Svtx") != std::string::npos)
      {
	buildFullTrack(ClusterKeyList, g4particle);
      }
    else
      {
	buildTrackSeed(ClusterKeyList, g4particle, _track_map);
      }
  }

  if (Verbosity() >= 5)
  {
    std::cout << "Loop over TrackMap " << _track_map_name << " entries " << _track_map->size() << std::endl;
    for (TrackSeedContainer::Iter iter = _track_map->begin();
         iter != _track_map->end(); ++iter)
    {
      TrackSeed* svtx_track = *iter;
      unsigned int trackid = _track_map->index(iter);
      svtx_track->identify();
    

      std::cout << "Track ID: " << trackid << ", Dummy Track pT: "
           << svtx_track->get_pt() << ", Truth Track/Particle ID: "
           << svtx_track->get_truth_track_id() 
	   << " (X,Y,Z) " << svtx_track->get_x() << ", " << svtx_track->get_y() << ", " << svtx_track->get_z()  
	   << std::endl ;
      std::cout << " nhits: " << svtx_track->size_cluster_keys()<< std::endl ;
      //Print associated clusters;
      ActsTransformations transformer;
      for (TrackSeed::ConstClusterKeyIter iter_clus =
               svtx_track->begin_cluster_keys();
           iter_clus != svtx_track->end_cluster_keys(); ++iter_clus)
      {
        TrkrDefs::cluskey cluster_key = *iter_clus;
	std::cout << "Key: "  << cluster_key<< std::endl;
        TrkrCluster* cluster = _cluster_map->findCluster(cluster_key);

        Acts::Vector3 global = transformer.getGlobalPosition(
          cluster_key, cluster,
          surfmaps,
          tgeometry);
        float radius = std::sqrt(square(global(0)) + square(global(1)));

        std::cout << "       cluster ID: "
		  << cluster_key << ", cluster radius: " << radius
		  << std::endl ;
      }

      trackid++;
    }
  }

  //==================================

  return Fun4AllReturnCodes::EVENT_OK;
}
void PHTruthTrackSeeding::buildFullTrack(std::vector<TrkrDefs::cluskey>& clusters,
					 PHG4Particle *g4particle)
{
  auto track = std::make_unique<SvtxTrackSeed_v1>();

  buildTrackSeed(clusters, g4particle, _tpc_seeds);

  /// The ids will by definition be the last entry in the container because
  /// the seeds were just added
  track->set_tpc_seed_index(_tpc_seeds->size()-1); 
 
  if(Verbosity() > 2)
    {
      std::cout << "adding svtxtrackseed " << std::endl;
      track->identify();
      auto tpcseed = _tpc_seeds->get(track->get_tpc_seed_index());
      tpcseed->identify();
 
    }

  _track_map->insert(track.get());

}
void PHTruthTrackSeeding::buildTrackSeed(std::vector<TrkrDefs::cluskey> clusters, PHG4Particle *g4particle, TrackSeedContainer* container)
{
  auto track = std::make_unique<TrackSeed_FastSim_v1>();
  bool silicon = false;
  for (const auto& cluskey : clusters){
    if( TrkrDefs::getTrkrId(cluskey) == TrkrDefs::TrkrId::mvtxId || 
	TrkrDefs::getTrkrId(cluskey) == TrkrDefs::TrkrId::inttId)
      { silicon = true; }
    track->insert_cluster_key(cluskey);
  }
  
  auto random = gsl_ran_flat(m_rng.get(), 0.95, 1.05);
  
  const auto particle = TDatabasePDG::Instance()->GetParticle(g4particle->get_pid());
  int charge = 1;
  if(particle) 
    { 
      if(particle->Charge() < 0)
	{ charge = -1; }
    }
  
  float px = g4particle->get_px() * random;
  float py = g4particle->get_py() * random;
  float pz = g4particle->get_pz() * random;
  const auto g4vertex = m_g4truth_container->GetVtx(g4particle->get_vtx_id());
  float x = g4vertex->get_x() * random; 
  float y = g4vertex->get_y() * random;
  float z = g4vertex->get_z() * random;

  float pt = sqrt(px*px+py*py);
  float phi = atan2(py,px);
  float R = 100 * pt / (0.3*1.4);
  float theta = atan2(pt,pz);
  if(theta < 0)
    { theta += M_PI; }
  if(theta > M_PI)
    { theta -= M_PI; }
  
  float eta = -log(tan(theta/2.));

  // We have two equations, phi = atan2(-(X0-x),y-Y0) and 
  //R^2 = (x-X0)^2 + (y-Y0)^2. Solve for X0 and Y0 knowing R and phi
  float tanphisq = square(tan(phi));
  float a = tanphisq + 1;
  float b =-2*y*(tanphisq+1);
  float c = (tanphisq+1)*square(y)-square(R);
  
  float Y0_1 = (-b + sqrt(square(b)-4*a*c)) / (2.*a);
  float Y0_2 = (-b - sqrt(square(b)-4*a*c)) / (2.*a);
  float X0_1 = sqrt(pow(R, 2) - pow(Y0_1 - y, 2)) + x;
  float X0_2 = -sqrt(pow(R, 2) - pow(Y0_2 - y, 2)) + x;
  track->set_X0(X0_1);
  track->set_Y0(Y0_1);
  track->set_qOverR(charge / R);
  track->set_slope(1. / tan(theta));
  track->set_Z0(z);
  
  /// Need to find the right one for the bend angle
  
  float newphi = track->get_phi(m_clusterMap, surfmaps, tgeometry);
  /// We have to pick the right one based on the bend angle, so iterate
  /// through until you find the closest phi match
  if( fabs(newphi-phi) > 0.03)
    {
      track->set_X0(X0_2);
      newphi = track->get_phi(m_clusterMap, surfmaps, tgeometry);
  
      if( fabs(newphi-phi) > 0.03)
	{
	  track->set_Y0(Y0_2);
	  newphi = track->get_phi(m_clusterMap, surfmaps, tgeometry);

	  if( fabs(newphi-phi) > 0.03)
	    {
	      track->set_X0(X0_1);
	      newphi = track->get_phi(m_clusterMap, surfmaps, tgeometry);
	    }
	}
    }
  
  if(Verbosity() > 2)
    {
      std::cout << "Charge is " << charge << std::endl;
      std::cout << "truth/reco px " << px << ", " << track->get_px(m_clusterMap, surfmaps, tgeometry) << std::endl;
      std::cout << "truth/reco py " << py << ", " << track->get_py(m_clusterMap, surfmaps, tgeometry) << std::endl;
      std::cout << "truth/reco pz " << pz << ", " << track->get_pz() << std::endl;
      std::cout << "truth/reco pt " << pt << ", " << track->get_pt() << std::endl;
      std::cout << "truth/reco phi " << phi << ", " << track->get_phi(m_clusterMap, surfmaps, tgeometry) << std::endl;
      std::cout << "truth/reco eta " << eta << ", " << track->get_eta() << std::endl;
      std::cout << "truth/reco x " << x << ", " << track->get_x() << std::endl;
      std::cout << "truth/reco y " << y << ", " << track->get_y() << std::endl;
      std::cout << "truth/reco z " << z << ", " << track->get_z() << std::endl;
	
    }
  
  // set intt crossing
  if(silicon)
    {
      // silicon tracklet
      /* inspired from PHtruthSiliconAssociation */
      const auto intt_crossings = getInttCrossings(track.get());
      if(intt_crossings.empty()) 
	{
	  if(Verbosity() > 1)  std::cout << "PHTruthTrackSeeding::Process - Silicon track " << container->size() - 1 << " has no INTT clusters" << std::endl;
	  return ;
	} else if( intt_crossings.size() > 1 ) {
	if(Verbosity() > 1) 
	  { std::cout << "PHTruthTrackSeeding::Process - INTT crossings not all the same for track " << container->size() - 1 << " crossing_keep - dropping this match " << std::endl; }
	
      } else {
	const auto& crossing = *intt_crossings.begin();
	track->set_crossing(crossing);
	if(Verbosity() > 1)
	  std::cout << "PHTruthTrackSeeding::Process - Combined track " << container->size() - 1  << " bunch crossing " << crossing << std::endl;           
      }
    }  // end if _min_layer
  else
    {
      // no INTT layers, crossing is unknown
      track->set_crossing(SHRT_MAX);	
    }

  container->insert(track.get());
}
int PHTruthTrackSeeding::CreateNodes(PHCompositeNode* topNode)
{
   // create nodes...
  PHNodeIterator iter(topNode);

  PHCompositeNode* dstNode = static_cast<PHCompositeNode*>(iter.findFirst(
      "PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cerr << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  PHNodeIterator iter_dst(dstNode);

  // Create the SVTX node
  PHCompositeNode* tb_node =
      dynamic_cast<PHCompositeNode*>(iter_dst.findFirst("PHCompositeNode",
                                                        "SVTX"));
  if (!tb_node)
  {
    tb_node = new PHCompositeNode("SVTX");
    dstNode->addNode(tb_node);
    if (Verbosity() > 0)
      std::cout << PHWHERE << "SVTX node added" << std::endl;
  }

  _tpc_seeds = findNode::getClass<TrackSeedContainer>(topNode,"TpcTrackSeedContainer");
  if(!_tpc_seeds)
    {
      _tpc_seeds = new TrackSeedContainer_v1;
      PHIODataNode<PHObject>* tracks_node = 
	new PHIODataNode<PHObject>(_tpc_seeds, "TpcTrackSeedContainer", "PHObject");
      tb_node->addNode(tracks_node);
    }

  _silicon_seeds = findNode::getClass<TrackSeedContainer>(topNode, "SiliconTrackSeedContainer");
   if(!_silicon_seeds)
    {
      _silicon_seeds = new TrackSeedContainer_v1;
      PHIODataNode<PHObject>* tracks_node = 
	new PHIODataNode<PHObject>(_silicon_seeds, "SiliconTrackSeedContainer", "PHObject");
      tb_node->addNode(tracks_node);
    }

   return Fun4AllReturnCodes::EVENT_OK;
}
int PHTruthTrackSeeding::GetNodes(PHCompositeNode* topNode)
{

  tgeometry = findNode::getClass<ActsTrackingGeometry>(topNode, "ActsTrackingGeometry");
  surfmaps = findNode::getClass<ActsSurfaceMaps>(topNode, "ActsSurfaceMaps");
  
  if(!tgeometry or !surfmaps) 
    {
      std::cerr << PHWHERE << "Error, can' find needed Acts nodes " << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  m_clusterMap = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  
  if(!m_clusterMap)
    {
      std::cerr << PHWHERE << "Error: Can't find node TRKR_CLUSTER" << std::endl ;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  m_cluster_crossing_map = findNode::getClass<TrkrClusterCrossingAssoc>(topNode, "TRKR_CLUSTERCROSSINGASSOC");
  if (!m_cluster_crossing_map)
  {
    std::cerr << PHWHERE << " ERROR: Can't find TRKR_CLUSTERCROSSINGASSOC " << std::endl ;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
    
  m_g4truth_container = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!m_g4truth_container)
  {
    std::cerr << PHWHERE << " ERROR: Can't find node G4TruthInfo" << std::endl ;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  hittruthassoc = findNode::getClass<TrkrHitTruthAssoc>(topNode, "TRKR_HITTRUTHASSOC");
  if (!hittruthassoc)
  {
    std::cout << PHWHERE << "Failed to find TRKR_HITTRUTHASSOC node, quit!" << std::endl ;
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
{ return 0; }


//_____________________________________________________________________________________________
std::set<short int> PHTruthTrackSeeding::getInttCrossings(TrackSeed *si_track) const
{
  std::set<short int> intt_crossings;

  // If the Si track contains an INTT hit, use it to get the bunch crossing offset
  // loop over associated clusters to get keys for silicon cluster
  
  for( auto iter = si_track->begin_cluster_keys(); 
       iter != si_track->end_cluster_keys(); ++iter)
  {
    
    const TrkrDefs::cluskey& cluster_key = *iter;
    const unsigned int trkrid = TrkrDefs::getTrkrId(cluster_key);
    if(trkrid == TrkrDefs::inttId)
    {
      
      // get layer from cluster key
      const unsigned int layer = TrkrDefs::getLayer(cluster_key);
      
      // get the bunch crossings for all hits in this cluster
      const auto crossings = m_cluster_crossing_map->getCrossings(cluster_key);
      for(auto iter = crossings.first; iter != crossings.second; ++iter)
      {
        const auto& [key, crossing] = *iter;
        if( Verbosity() )
        { std::cout << "PHTruthTrackSeeding::getInttCrossings - si Track cluster " << key << " layer " << layer << " crossing " << crossing  << std::endl; }
        intt_crossings.insert(crossing);
      }
    }
  }
  
  return intt_crossings;
}
