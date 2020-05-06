#include "PHActsTracks.h"

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
  , m_actsProtoTracks(nullptr)
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

  for (SvtxTrackMap::Iter trackIter = m_trackMap->begin();
       trackIter != m_trackMap->end(); ++trackIter)
  {
    const SvtxTrack *track = trackIter->second;

    if (!track)
      continue;

    if (Verbosity() > 1)
    {
      std::cout << "found SvtxTrack " << trackIter->first << std::endl;
      track->identify();
    }

    /// Get the necessary parameters and values for the TrackParameters
    const Acts::BoundSymMatrix seedCov = getActsCovMatrix(track);
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
    m_actsProtoTracks->push_back(actsTrack);
  }

  if (Verbosity() > 20)
    std::cout << "Finished PHActsTrack::process_event" << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

Acts::BoundSymMatrix PHActsTracks::getActsCovMatrix(const SvtxTrack *track)
{
  Acts::BoundSymMatrix matrix = Acts::BoundSymMatrix::Zero();
  
  const double px = track->get_px();
  const double py = track->get_py();
  const double pz = track->get_pz();
  const double p = sqrt(px * px + py * py + pz * pz);
  const double phiPos = atan2(track->get_x(), track->get_y());
  const int charge = track->get_charge();

  // Get the track seed covariance matrix
  // These are the variances, so the std devs are sqrt(seedCov[i][j])
  Acts::BoundSymMatrix seedCov = Acts::BoundSymMatrix::Zero();

  for (int i = 0; i < 6; i++)
  {
    for (int j = 0; j < 6; j++)
    {
      /// Track covariance matrix is in basis (x,y,z,px,py,pz). Need to put
      /// it in form of (x,y,px,py,pz,time) for acts
   int row = -1;
      int col = -1;
      if( i < 2)
	row = i;
      else if ( i < 5)
	row = i+1;
      else if (i == 5)
	row = 2;

      if( j < 2 )
	col = j;
      else if ( j < 5 )
	col = j+1;
      else if ( j == 5 )
	col = 2;

      seedCov(i,j) = track->get_error(row, col);
      
      /// get the units right, since acts works in mm
      /// Don't need to multiply by units of GeV since that is Acts default
      if(i < 2 && j < 2)
	seedCov(i,j) *= Acts::UnitConstants::cm2;
      else if (i < 2 && j < 5 )
	seedCov(i,j) *= Acts::UnitConstants::cm;
      else if (i < 5 && j < 2 )
	seedCov(i,j) *= Acts::UnitConstants::cm;

    }
  }
  
  /// Set the time column/row to 0
  for(int i = 5; i< 6; i++)
    {
      for(int j = 0; j <6; j++)
	{
	  seedCov(i,j) = 0;
	  seedCov(j,i) = 0;
	}
    }
  /// Just give time a 10 ns covariance for now
  seedCov(5,5) = 10 * Acts::UnitConstants::ns;

  if(Verbosity() > 10)
    {
      std::cout << "Initial covariance: " << std::endl;
      for(int i =0; i<6; i++)
	{
	  std::cout << std::endl;
	  for(int j= 0; j<6; j++)
	    {
	      std::cout << track->get_error(i,j) << ", ";
	    }

	}

      std::cout << std::endl << "Shifted covariance " << std::endl;
      for(int i = 0; i < 6; i++)
	{
	  std::cout<<std::endl;
	  for(int j=0; j<6; j++)
	    {
	      std::cout << seedCov(i,j) << ", ";

	    }
	}
    }


  /// Need to transform from global to local coordinate frame. 
  /// Amounts to the local transformation as in PHActsSourceLinks as well as
  /// a rotation from cartesian to spherical coordinates for the momentum
  /// Rotating from (x_G, y_G, px, py, pz, time) to (x_L, y_L, phi, theta, q/p,time)

  /// Make a unit p vector for the rotation
  const double uPx = px / p;
  const double uPy = py / p;
  const double uPz = pz / p;
  const double uP = sqrt(uPx * uPx + uPy * uPy + uPz * uPz);
  
  /// This needs to rotate to (x_L, y_l, phi, theta, q/p, t)
  Acts::BoundSymMatrix rotation = Acts::BoundSymMatrix::Zero();

  /// Local position rotations
  rotation(0,0) = cos(phiPos);
  rotation(0,1) = sin(phiPos);
  rotation(1,0) = -1 * sin(phiPos);
  rotation(1,1) = cos(phiPos);

  /// Momentum vector rotations
  /// phi rotation
  rotation(2,2) = -1 * uPy / (uPx * uPx + uPy * uPy);
  rotation(2,3) = -1 * uPx / (uPx * uPx + uPy * uPy);

  /// theta rotation
  /// Leave uP in for clarity, even though it is trivially unity
  rotation(3,2) = (uPx * uPz) / (uP * uP * sqrt( uPx * uPx + uPy * uPy) );
  rotation(3,3) = (uPy * uPz) / (uP * uP * sqrt( uPx * uPx + uPy * uPy) );
  rotation(3,4) = (-1 * sqrt(uPx * uPx + uPy * uPy)) / (uP * uP);
  
  /// p rotation
  rotation(4,2) = uPx / uP;
  rotation(4,3) = uPy / uP;
  rotation(4,4) = uPz / uP;

  /// time rotation
  rotation(5,5) = 1;

  if(Verbosity() > 10 )
    {
      std::cout << std::endl << "Rotation matrix is : " << std::endl;
      for(int i = 0; i < 6; i++)
	{
	  std::cout << std::endl;
	  for(int j = 0 ; j < 6 ; j++)
	    {
	      std::cout<< rotation(i,j) << ", ";
	    }
	}

    }

  /// Rotate the covariance matrix by the jacobian rotation matrix
  matrix = rotation * seedCov * rotation.transpose();

  if(Verbosity() > 10)
    {
      std::cout << std::endl << "Acts rotated covariance " << std::endl;
      for(int i = 0; i < 6; i++)
	{
	  for(int j=0; j<6; j++)
	    {
	      std::cout << matrix(i,j) << ", ";
	      
	    }
	  
	  std::cout << std::endl;
	}
    }

  /// Rotate again to get q/p instead of p
  Acts::BoundSymMatrix qprotation = Acts::BoundSymMatrix::Zero();
  qprotation(0,0) = 1;
  qprotation(1,1) = 1;
  qprotation(2,2) = 1;
  qprotation(3,3) = 1;
  qprotation(4,4) = charge * charge / (p * p * p * p);
  qprotation(5,5) = 1;
  
  Acts::BoundSymMatrix finalMatrix = Acts::BoundSymMatrix::Zero();
  finalMatrix = qprotation * matrix * qprotation.transpose();
  
    if(Verbosity() > 10)
    {
      std::cout << std::endl << "Acts rotated rotated covariance " << std::endl;
      for(int i = 0; i < 6; i++)
	{
	  std::cout << std::endl;
	  for(int j=0; j<6; j++)
	    {
	      std::cout << finalMatrix(i,j) << ", ";
	      
	    }
	}
    }
  
  return finalMatrix;
}

void PHActsTracks::createNodes(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));

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

  m_actsProtoTracks = findNode::getClass<std::vector<ActsTrack>>(topNode, "ActsProtoTracks");

  if (!m_actsProtoTracks)
  {
    m_actsProtoTracks = new std::vector<ActsTrack>;

    PHDataNode<std::vector<ActsTrack>> *protoTrackNode =
        new PHDataNode<std::vector<ActsTrack>>(m_actsProtoTracks, "ActsProtoTracks");

    svtxNode->addNode(protoTrackNode);
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
