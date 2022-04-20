
#include "TrackSeedTrackMapConverter.h"

#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackMap_v1.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrack_v3.h>
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
  unsigned int trackid = 0;
  for(const auto& trackSeed : *m_seedContainer)
    {
      auto svtxtrack = std::make_unique<SvtxTrack_v3>();

      if(Verbosity() > 0)
	{
	  std::cout << "iterating track seed " << trackid << std::endl;
	  trackSeed->identify();
	}

      trackid++;

      svtxtrack->set_x(trackSeed->get_x());
      svtxtrack->set_y(trackSeed->get_y());
      svtxtrack->set_z(trackSeed->get_z());
      svtxtrack->set_charge( trackSeed->get_qOverR() > 0 ? 1 : -1);
      svtxtrack->set_px(trackSeed->get_px());
      svtxtrack->set_py(trackSeed->get_py());
      svtxtrack->set_pz(trackSeed->get_pz());

      for(TrackSeed::ConstClusterKeyIter iter = trackSeed->begin_cluster_keys();
	  iter != trackSeed->end_cluster_keys();
	  ++iter)
	{
	  svtxtrack->insert_cluster_key(*iter);
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

  return Fun4AllReturnCodes::EVENT_OK;
}
