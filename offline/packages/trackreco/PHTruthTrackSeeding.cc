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

#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/PHRandomSeed.h>

#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterCrossingAssoc.h> 
#include <trackbase/TrkrHitTruthAssoc.h>
#include <trackbase_historic/ActsTransformations.h>
#include <trackbase_historic/SvtxTrack.h>     // for SvtxTrack, SvtxTra...
#include <trackbase_historic/SvtxTrackMap.h>  // for SvtxTrackMap, Svtx...
#include <trackbase_historic/SvtxTrack_FastSim_v3.h>
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

    auto svtx_track = std::make_unique<TrackSeed_FastSim_v1>();
    svtx_track->set_truth_track_id(gtrackID);
      
    // Smear the truth values out by 5% so that the seed momentum and
    // position aren't completely exact
    auto random = gsl_ran_flat(m_rng.get(), 0.95, 1.05);
    // make it negative sometimes
    if(rand() % 2)
      random *= -1;

    for (const auto& cluskey : ClusterKeyList){
      svtx_track->insert_cluster_key(cluskey);
    }

    svtx_track->circleFitByTaubin(m_clusterMap, surfmaps, tgeometry,
				  _min_layer, _max_layer);
    svtx_track->lineFit(m_clusterMap, surfmaps, tgeometry,
			_min_layer, _max_layer);

    // set intt crossing
    if(_min_layer < 7)
      {
	// silicon tracklet
	/* inspired from PHtruthSiliconAssociation */
	const auto intt_crossings = getInttCrossings(svtx_track.get());
	if(intt_crossings.empty()) 
	  {
	    if(Verbosity() > 1)  std::cout << "PHTruthTrackSeeding::Process - Silicon track " << gtrackID << " has no INTT clusters" << std::endl;
	    continue ;
	  } else if( intt_crossings.size() > 1 ) {
	  
	  if(Verbosity() > 1) 
	    { std::cout << "PHTruthTrackSeeding::Process - INTT crossings not all the same for track " << gtrackID << " crossing_keep - dropping this match " << std::endl; }
	  
	} else {
	  
	  const auto& crossing = *intt_crossings.begin();
	  svtx_track->set_crossing(crossing);
	  if(Verbosity() > 1)
	    std::cout << "PHTruthTrackSeeding::Process - Combined track " << gtrackID  << " bunch crossing " << crossing << std::endl;           
	}
      }  // end if _min_layer
    else
      {
	// no INTT layers, crossing is unknown
	svtx_track->set_crossing(SHRT_MAX);	
      }
 
    _track_map->insert(svtx_track.get());
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
