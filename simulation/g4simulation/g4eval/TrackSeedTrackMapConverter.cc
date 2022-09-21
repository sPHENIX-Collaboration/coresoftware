
#include "TrackSeedTrackMapConverter.h"

#include <trackbase/ActsGeometry.h>
#include <trackbase/TrkrClusterContainer.h>

#include <trackbase_historic/SvtxTrackMap_v1.h>
#include <trackbase_historic/SvtxTrack_v4.h>
#include <trackbase_historic/TrackSeedContainer.h>
#include <trackbase_historic/TrackSeed.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>                        // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>                      // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h> 

//____________________________________________________________________________..
TrackSeedTrackMapConverter::TrackSeedTrackMapConverter(const std::string &name):
 SubsysReco(name)
{
}

//____________________________________________________________________________..
TrackSeedTrackMapConverter::~TrackSeedTrackMapConverter()
{
}

//____________________________________________________________________________..
int TrackSeedTrackMapConverter::InitRun(PHCompositeNode *topNode)
{
  int ret = getNodes(topNode);
  return ret;
}

//____________________________________________________________________________..
int TrackSeedTrackMapConverter::process_event(PHCompositeNode*)
{
  if(Verbosity() > 1)
    {
      std::cout <<"silicon seed map size " << m_siContainer->size() << std::endl;
      for(auto iter = m_siContainer->begin(); iter != m_siContainer->end();
	  ++iter)
	{
	  auto seed = *iter;
	  if(!seed){
	    std::cout << "no seed at index "<< m_siContainer->index(iter) 
		      << std::endl;
	    continue;
	  }
	  seed->identify();
	}
      if(m_tpcContainer)
	{
	  std::cout << "tpc seed map size " << m_tpcContainer->size() << std::endl;
	  for(auto iter = m_tpcContainer->begin(); iter != m_tpcContainer->end();
	      ++iter)
	    {
	      auto seed = *iter;
	      if(!seed) {
		std::cout << "no tpc seed at entry " << m_tpcContainer->index(iter) 
			  << std::endl;
		continue;
	      }
	      seed->identify();
	    }
	}
    }

  unsigned int trackid = 0;
  for(const auto& trackSeed : *m_seedContainer)
    {
      /// If the seed was removed, skip it
      if(!trackSeed)
	{ continue; }

      if(m_trackSeedName.find("SvtxTrackSeed") != std::string::npos)
	{
	  /// Catches entries in the vector removed by ghost finder
	  unsigned int tpcindex = trackSeed->get_tpc_seed_index();
	  TrackSeed* seed = m_tpcContainer->get(tpcindex);
	  if(!seed)
	    { continue; }
	}

      auto svtxtrack = std::make_unique<SvtxTrack_v4>();

      if(Verbosity() > 0)
	{
	  std::cout << "iterating track seed " << trackid << std::endl;
	}

      svtxtrack->set_id(trackid);
      trackid++;

      /// If we've run the track matching
      if(m_trackSeedName.find("SvtxTrackSeed") != std::string::npos)
	{
	  if(Verbosity() > 0)
	    {
	      std::cout << "tpc seed id " << trackSeed->get_tpc_seed_index() <<std::endl;
	      std::cout << "si seed id " << trackSeed->get_silicon_seed_index() << std::endl;
	    }

	  unsigned int seedindex = trackSeed->get_tpc_seed_index();
	  TrackSeed *tpcseed = m_tpcContainer->get(seedindex);
	  if(trackSeed->get_silicon_seed_index() == std::numeric_limits<unsigned int>::max())
	    {      
	      /// Didn't find a match, so just use the tpc seed
	      svtxtrack->set_x(tpcseed->get_x());
	      svtxtrack->set_y(tpcseed->get_y());
	      svtxtrack->set_z(tpcseed->get_z()); 
	    }
	  else
	    {
	      TrackSeed *siseed = m_siContainer->get(trackSeed->get_silicon_seed_index());
	      svtxtrack->set_x(siseed->get_x());
	      svtxtrack->set_y(siseed->get_y());
	      svtxtrack->set_z(siseed->get_z());
	      addKeys(svtxtrack, siseed);
	      svtxtrack->set_silicon_seed(siseed);
	    }


	  svtxtrack->set_charge( tpcseed->get_qOverR() > 0 ? 1 : -1);
	  svtxtrack->set_px(tpcseed->get_px(m_clusters,m_tGeometry));
	  svtxtrack->set_py(tpcseed->get_py(m_clusters,m_tGeometry));
	  svtxtrack->set_pz(tpcseed->get_pz());
	  
	  addKeys(svtxtrack, tpcseed);
	  svtxtrack->set_tpc_seed(tpcseed);
	

	}
      else
	{
	  /// Otherwise we are using an individual subdetectors container
	  svtxtrack->set_x(trackSeed->get_x());
	  svtxtrack->set_y(trackSeed->get_y());
	  svtxtrack->set_z(trackSeed->get_z());
	  svtxtrack->set_charge( trackSeed->get_qOverR() > 0 ? 1 : -1);
	  svtxtrack->set_px(trackSeed->get_px(m_clusters,m_tGeometry));
	  svtxtrack->set_py(trackSeed->get_py(m_clusters,m_tGeometry));
	  svtxtrack->set_pz(trackSeed->get_pz());
	  
	  addKeys(svtxtrack, trackSeed);
	  if(m_trackSeedName.find("SiliconTrackSeed") != std::string::npos)
	    {
	      svtxtrack->set_silicon_seed(trackSeed);
	      svtxtrack->set_tpc_seed(nullptr);
	    }
	  else if(m_trackSeedName.find("TpcTrackSeed") != std::string::npos)
	    {
	      svtxtrack->set_tpc_seed(trackSeed);
	      svtxtrack->set_silicon_seed(nullptr);
	    }
	}

      if(Verbosity() > 0)
	{
	  std::cout << "Inserting svtxtrack into map " << std::endl;
	  svtxtrack->identify();
	}
      
      m_trackMap->insert(svtxtrack.get());

    }

  return Fun4AllReturnCodes::EVENT_OK;
}


//____________________________________________________________________________..
int TrackSeedTrackMapConverter::End(PHCompositeNode*)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

void TrackSeedTrackMapConverter::addKeys(std::unique_ptr<SvtxTrack_v4>& track, TrackSeed *seed)
{
  for(TrackSeed::ConstClusterKeyIter iter = seed->begin_cluster_keys();
      iter != seed->end_cluster_keys();
      ++iter)
    {
      track->insert_cluster_key(*iter);
    }
}


int TrackSeedTrackMapConverter::getNodes(PHCompositeNode *topNode)
{
  m_trackMap = findNode::getClass<SvtxTrackMap>(topNode,"SvtxTrackMap");
  if(!m_trackMap)
    {
      // create it
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
	    { std::cout << PHWHERE << "SVTX node added" << std::endl; }
	}
      
      
      m_trackMap = findNode::getClass<SvtxTrackMap>(topNode, m_trackMapName);
      if (!m_trackMap)
	{
	  m_trackMap = new SvtxTrackMap_v1;
	  PHIODataNode<PHObject>* tracks_node = 
	    new PHIODataNode<PHObject>(m_trackMap, m_trackMapName, "PHObject");
	  tb_node->addNode(tracks_node);
	  if (Verbosity() > 0){
	    std::cout << PHWHERE << "Svtx/" << m_trackMapName  << " node added" << std::endl;
	  }
	}
    }
  
  m_seedContainer = findNode::getClass<TrackSeedContainer>(topNode, m_trackSeedName);
  if(!m_seedContainer)
    {
      std::cout << PHWHERE << " Can't find track seed container " << m_trackSeedName << ", can't continue."
		<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  m_tpcContainer = findNode::getClass<TrackSeedContainer>(topNode, "TpcTrackSeedContainer");
  if(!m_tpcContainer)
    {
      std::cout << PHWHERE << "WARNING, TrackSeedTrackMapConverter may seg fault depending on what seeding algorithm this is run after" << std::endl;
    }

  m_siContainer = findNode::getClass<TrackSeedContainer>(topNode, "SiliconTrackSeedContainer");
  if(!m_siContainer)
    {
        std::cout << PHWHERE << "WARNING, TrackSeedTrackMapConverter may seg fault depending on what seeding algorithm this is run after" << std::endl;
    }

  m_clusters = findNode::getClass<TrkrClusterContainer>(topNode,"TRKR_CLUSTER");
  if(!m_clusters)
    {
      std::cout << PHWHERE << " Can't find cluster container, can't continue."
		<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  
  m_tGeometry = findNode::getClass<ActsGeometry>(topNode,"ActsGeometry");
  if(!m_tGeometry)
    {
      std::cout << PHWHERE << " Can't find ActsGeometry, can't continue."
		<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  return Fun4AllReturnCodes::EVENT_OK;
}
