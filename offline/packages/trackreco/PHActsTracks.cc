#include "PHActsTracks.h"
#include "ActsCovarianceRotater.h"

/// Fun4All includes
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <Acts/EventData/ChargePolicy.hpp>
#include <Acts/EventData/SingleCurvilinearTrackParameters.hpp>
#include <Acts/Utilities/Units.hpp>

/// Tracking includes
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>

/// std (and the like) includes
#include <cmath>
#include <iostream>
#include <memory>
#include <utility>

#include <TMatrixDSym.h>


PHActsTracks::PHActsTracks(const std::string &name)
  : SubsysReco(name)
  , m_actsTrackMap(nullptr)
  , m_trackMap(nullptr)
  , m_hitIdClusKey(nullptr)
  , m_sourceLinks(nullptr)
  , m_tGeometry(nullptr)
{
  Verbosity(0);
}

int PHActsTracks::End(PHCompositeNode *topNode)
{

  if (Verbosity() > 10)
  {
    std::cout << "Finished PHActsTracks" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsTracks::Init(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsTracks::InitRun(PHCompositeNode *topNode)
{
  createNodes(topNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsTracks::process_event(PHCompositeNode *topNode)
{
  if (Verbosity() > 1)
  {
    std::cout << "Start process_event in PHActsTracks" << std::endl;
  }

  /// Check to get the nodes needed
  if (getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
    return Fun4AllReturnCodes::ABORTEVENT;

  /// Vector to hold source links for a particular track
  std::vector<SourceLink> trackSourceLinks;
  std::vector<FW::TrackParameters> trackSeeds;

  ActsCovarianceRotater *rotater = new ActsCovarianceRotater();
  
  for (SvtxTrackMap::Iter trackIter = m_trackMap->begin();
       trackIter != m_trackMap->end(); ++trackIter)
  {
    const SvtxTrack *track = trackIter->second;
    const unsigned int trackKey = trackIter->first;
    if (!track)
      continue;

    if (Verbosity() > 1)
    {
      std::cout << "found SvtxTrack " << trackIter->first << std::endl;
      track->identify();
    }

    /// Get the necessary parameters and values for the TrackParameters
    const Acts::BoundSymMatrix seedCov = 
          rotater->rotateSvtxTrackCovToActs(track);
    const Acts::Vector3D seedPos(track->get_x()  * Acts::UnitConstants::cm,
                                 track->get_y()  * Acts::UnitConstants::cm,
                                 track->get_z()  * Acts::UnitConstants::cm);
    const Acts::Vector3D seedMom(track->get_px() * Acts::UnitConstants::GeV,
                                 track->get_py() * Acts::UnitConstants::GeV,
                                 track->get_pz() * Acts::UnitConstants::GeV);

    // just set to 10 ns for now?
    const double trackTime = 10 * Acts::UnitConstants::ns;
    const int trackQ = track->get_charge();

    const FW::TrackParameters trackSeed(seedCov, seedPos, seedMom, 
					trackQ * Acts::UnitConstants::e, 
					trackTime);

    /// Start fresh for this track
    trackSourceLinks.clear();
    for (SvtxTrack::ConstClusterKeyIter clusIter = track->begin_cluster_keys();
         clusIter != track->end_cluster_keys();
         ++clusIter)
    {
      const TrkrDefs::cluskey key = *clusIter;

      const unsigned int hitId = m_hitIdClusKey->find(key)->second;
 
      trackSourceLinks.push_back(m_sourceLinks->find(hitId)->second);
      
      if (Verbosity() > 100)
	{

	  
	  std::cout << std::endl << "cluskey " << key
		    << " has hitid " << hitId
		    << std::endl;
	  std::cout << "Adding the following surface for this SL" << std::endl;
	  m_sourceLinks->find(hitId)->second.referenceSurface().toStream(m_tGeometry->geoContext, std::cout);
	}
    }

    if (Verbosity() > 0)
    {
      for (unsigned int i = 0; i < trackSourceLinks.size(); ++i)
      {
        std::cout << "proto_track readback: hitid " << trackSourceLinks.at(i).hitID() << std::endl;
      }
    }

    ActsTrack actsTrack(trackSeed, trackSourceLinks);
    m_actsTrackMap->insert(std::pair<unsigned int, ActsTrack>(trackKey, actsTrack));
  }

  if (Verbosity() > 20)
    std::cout << "Finished PHActsTrack::process_event" << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}


int PHActsTracks::ResetEvent(PHCompositeNode *topNode)
{
  /// Reset the proto track vector after each event
  m_actsTrackMap->clear();
  return Fun4AllReturnCodes::EVENT_OK;
}


void PHActsTracks::createNodes(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));

  if (!dstNode)
  {
    std::cerr << "DST node is missing, quitting" << std::endl;
    throw std::runtime_error("Failed to find DST node in PHActsTracks::createNodes");
  }

  PHCompositeNode *svtxNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "SVTX"));

  if (!svtxNode)
  {
    svtxNode = new PHCompositeNode("SVTX");
    dstNode->addNode(svtxNode);
  }

  m_actsTrackMap = findNode::getClass<std::map<unsigned int, ActsTrack>>(topNode, "ActsTrackMap");
  
  if(!m_actsTrackMap)
    {
      m_actsTrackMap = new std::map<unsigned int, ActsTrack>;
      PHDataNode<std::map<unsigned int, ActsTrack>> *actsTrackMapNode = 
	new PHDataNode<std::map<unsigned int, ActsTrack>>(m_actsTrackMap, "ActsTrackMap");
      svtxNode->addNode(actsTrackMapNode);
    }


  return;
}

int PHActsTracks::getNodes(PHCompositeNode *topNode)
{
  m_trackMap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");

  if (!m_trackMap)
  {
    std::cout << PHWHERE << "SvtxTrackMap not found on node tree. Exiting."
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_sourceLinks = findNode::getClass<std::map<unsigned int, SourceLink>>(topNode, "TrkrClusterSourceLinks");

  if (!m_sourceLinks)
  {
    std::cout << PHWHERE << "TrkrClusterSourceLinks node not found on node tree. Exiting."
              << std::endl;

    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_tGeometry = findNode::getClass<ActsTrackingGeometry>(topNode, "ActsTrackingGeometry");

  if(!m_tGeometry)
    {
      std::cout << PHWHERE << "ActsTrackingGeometry not on node tree, exiting."
		<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  m_hitIdClusKey = findNode::getClass<std::map<TrkrDefs::cluskey, unsigned int>>(topNode, "HitIDClusIDActsMap");

  if (!m_hitIdClusKey)
  {
    std::cout << PHWHERE << "HitID cluster key map not found on node tree. Exiting. "
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
