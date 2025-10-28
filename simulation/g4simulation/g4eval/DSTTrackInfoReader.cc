/*!
 * \file DSTTrackInfoReader.cc
 * \author Alex Patton <aopatton@mit.edu>
 */

#include "DSTTrackInfoReader.h"

#include <trackbase_historic/TrackInfoContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>


//_____________________________________________________________________

//_____________________________________________________________________
DSTTrackInfoReader::DSTTrackInfoReader(const std::string& name)
  : SubsysReco(name)
{
}

//_____________________________________________________________________
int DSTTrackInfoReader::InitRun(PHCompositeNode* topNode)
{
  // make a clusterContainer node and a trackmap node
  // assume writing macro removed them
  //  evalNode->addNode(newNode);

  // auto newNode = new PHIODataNode<PHObject>( new DSTContainerv3, "DSTContainer","PHObject");
  // evalNode->addNode(newNode);

  // TrackInfo container
  m_track_info_container = findNode::getClass<TrackInfoContainer>(topNode, "TrackInfoContainer");

  // m_container = findNode::getClass<DSTContainer>(topNode, "DSTContainer");

  // svtxtrackmap constructer is protected
  // auto svtxNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "SVTX"));
  // auto newNode = new PHIODataNode<PHObject>( new SvtxTrackMap, "SvtxTrackMap_2","PHObject");
  // svtxNode->addNode(newNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int DSTTrackInfoReader::process_event(PHCompositeNode* topNode)
{
  // Init(topNode);
  //  load nodes
  auto res = load_nodes(topNode);
  // std::cout <<"This is before I fill cluster map" << std::endl;
  // m_cluster_map->identify();
  if (res != Fun4AllReturnCodes::EVENT_OK)
  {
    return res;
  }
  evaluate_track_info();

  m_track_info_container->Reset();

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int DSTTrackInfoReader::load_nodes(PHCompositeNode* topNode)
{
  // get necessary nodes
  // m_track_map = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");

  m_track_info_container = findNode::getClass<TrackInfoContainer>(topNode, "TrackInfoContainer");

  return Fun4AllReturnCodes::EVENT_OK;
}

void DSTTrackInfoReader::evaluate_track_info()
{
  if (Verbosity() > 1)
  {
    for (int iTrk = 0; iTrk < (int) m_track_info_container->size(); iTrk++)
    {
      std::cout << "size: " << m_track_info_container->size() << std::endl;
      std::cout << "hitbitmap: " << (m_track_info_container->get_trackinfo(iTrk))->get_hitbitmap() << std::endl;
      std::cout << "crossing: " << (m_track_info_container->get_trackinfo(iTrk))->get_crossing() << std::endl;
      std::cout << "chi^2: " << (m_track_info_container->get_trackinfo(iTrk))->get_chisq() << std::endl;
      std::cout << "NDF: " << unsigned((m_track_info_container->get_trackinfo(iTrk))->get_ndf()) << std::endl;
    }
  }
}
