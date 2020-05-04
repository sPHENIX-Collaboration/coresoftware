
#include "PHActsTrkProp.h"
#include "MakeActsGeometry.h"
#include "ActsTrack.h"

/// Fun4All includes
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <phool/phool.h>

/// Tracking includes
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>

#include <Acts/EventData/ChargePolicy.hpp>
#include <Acts/EventData/SingleCurvilinearTrackParameters.hpp>
#include <Acts/EventData/TrackParameters.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>

#include <Acts/Propagator/Propagator.hpp>
#include <Acts/Propagator/Navigator.hpp>
#include <Acts/Propagator/AbortList.hpp>
#include <Acts/Propagator/ActionList.hpp>
#include <Acts/Utilities/Helpers.hpp>
#include <Acts/Utilities/Units.hpp>
#include <Acts/TrackFinder/CombinatorialKalmanFilter.hpp>

#include <ACTFW/Fitting/TrkrClusterFindingAlgorithm.hpp>
#include <ACTFW/Plugins/BField/BFieldOptions.hpp>
#include <ACTFW/Plugins/BField/ScalableBField.hpp>
#include <ACTFW/Framework/ProcessCode.hpp>
#include <ACTFW/Framework/WhiteBoard.hpp>
#include <ACTFW/EventData/Track.hpp>
#include <ACTFW/Framework/AlgorithmContext.hpp>


#include <TMatrixDSym.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <utility>

PHActsTrkProp::PHActsTrkProp(const std::string& name)
  : PHTrackPropagating(name)
  , m_event(0)
  , m_trackMap(nullptr)
  , m_actsProtoTracks(nullptr)
  , m_sourceLinks(nullptr)
{
  Verbosity(0);
}

 PHActsTrkProp::~PHActsTrkProp()
{
}

int PHActsTrkProp::Setup(PHCompositeNode* topNode)
{
  createNodes(topNode);
  
  if (getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
    return Fun4AllReturnCodes::ABORTEVENT;

  // Implement different chi2 criteria for different pixel (volumeID: 2)
  // layers:
  m_sourceLinkSelectorConfig.layerMaxChi2 = {{2, {{2, 8}, {4, 7}}}};
  // Implement different chi2 criteria for pixel (volumeID: 2) and strip
  // (volumeID: 3):
  m_sourceLinkSelectorConfig.volumeMaxChi2 = {{2, 7}, {3, 8}};
  m_sourceLinkSelectorConfig.maxChi2 = 8;
  // Set the allowed maximum number of source links to be large enough
  m_sourceLinkSelectorConfig.maxNumSourcelinksOnSurface = 100;
 
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsTrkProp::Process()
{

  m_event++;

  if (Verbosity() > 1)
  {
    std::cout << PHWHERE << "Events processed: " << m_event << std::endl;
    std::cout << "Start PHActsTrkProp::process_event" << std::endl;
  }

  FW::TrkrClusterFindingAlgorithm::Config findCfg;
  findCfg.finder = FW::TrkrClusterFindingAlgorithm::makeFinderFunction(
                   m_tGeometry->tGeometry,
		   m_tGeometry->magField,
		   Acts::Logging::VERBOSE);

  PerigeeSurface pSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>
    (Acts::Vector3D(0., 0., 0.));

  std::vector<SourceLink> sourceLinks;

  std::map<unsigned int, SourceLink>::iterator slIter = m_sourceLinks->begin();
  while(slIter != m_sourceLinks->end())
    {
      sourceLinks.push_back(slIter->second);
      ++slIter;
    }

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
      const Acts::Vector3D seedPos(track->get_x(),
				   track->get_y(),
				   track->get_z());
      const Acts::Vector3D seedMom(track->get_px(),
				   track->get_py(),
				   track->get_pz());
      
      // just set to 0 for now?
      const double trackTime = 0;
      const int trackQ = track->get_charge();
      
      const FW::TrackParameters trackSeed(seedCov, seedPos,
					  seedMom, trackQ, trackTime);
      
      
      Acts::CombinatorialKalmanFilterOptions<SourceLinkSelector> ckfOptions(
	        m_tGeometry->geoContext, 
		m_tGeometry->magFieldContext, 
		m_tGeometry->calibContext, 
		m_sourceLinkSelectorConfig, 
		&(*pSurface));
      
      auto result = findCfg.finder(sourceLinks, trackSeed, ckfOptions);
      
      if(result.ok())
	{
	  const auto& fitOutput = result.value();
	  auto parameterMap = fitOutput.fittedParameters;
	  /*
	  for(auto element : parameterMap)
	    {
	      if(element.second)
		{
		  const auto& params = element.second.value();
		}

	    }
	  */
	}
    }


  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsTrkProp::End()
{
  if (Verbosity() > 10)
  {
    std::cout << "Finished PHActsTrkProp" << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

Acts::BoundSymMatrix PHActsTrkProp::getActsCovMatrix(const SvtxTrack *track)
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
  for (int i = 0; i < 5; i++)
  {
    for (int j = 0; j < 6; j++)
    {
      /// Track covariance matrix is in basis (x,y,z,px,py,pz). Need to put
      /// it in form of (x,y,px,py,pz,time) for acts
      if(i < 2) /// get x,y components
	seedCov(i, j) = track->get_error(i, j);
      else if(i < 5) /// get px,py,pz components 1 row up
	seedCov(i,j) = track->get_error(i+1, j);
      else if (i == 5) /// Get z components, scale by drift velocity
	seedCov(i,j) = track->get_error(2, j) * 8.; //cm per millisec drift vel
    }
  }
  std::cout<<track->get_x()<<std::endl;
  /// convert the global z position covariances to timing covariances
  /// TPC z position resolution is 0.05 cm, drift velocity is 8cm/ms
  /// --> therefore timing resolution is ~6 microseconds
  seedCov(5,5) = 6 * Acts::UnitConstants::us;

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
  rotation(2,3) = -1 * uPy / (uPx * uPx + uPy * uPy);
  rotation(2,4) = -1 * uPx / (uPx * uPx + uPy * uPy);

  /// theta rotation
  /// Leave uP in for clarity, even though it is trivially unity
  rotation(3,3) = (uPx * uPz) / (uP * uP * sqrt( uPx * uPx + uPy * uPy) );
  rotation(3,4) = (uPy * uPz) / (uP * uP * sqrt( uPx * uPx + uPy * uPy) );
  rotation(3,5) = (-1 * sqrt(uPx * uPx + uPy * uPy)) / (uP * uP);
  
  /// q/p rotation
  rotation(4,3) = charge / uPx;
  rotation(4,4) = charge / uPy;
  rotation(4,5) = charge / uPz;

  /// time rotation
  rotation(5,5) = 1;

  matrix = rotation * seedCov * rotation.transpose();

  return matrix;
}

void PHActsTrkProp::createNodes(PHCompositeNode* topNode)
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

  if(!m_actsProtoTracks)
    {
      m_actsProtoTracks = new std::vector<ActsTrack>;

      PHDataNode<std::vector<ActsTrack>> *protoTrackNode =
        new PHDataNode<std::vector<ActsTrack>>(m_actsProtoTracks, "ActsProtoTracks");
      
      svtxNode->addNode(protoTrackNode);

    }


  return;
}

int PHActsTrkProp::getNodes(PHCompositeNode* topNode)
{

  m_sourceLinks = findNode::getClass<std::map<unsigned int, SourceLink>>(topNode, "TrkrClusterSourceLinks");

  if (!m_sourceLinks)
    {
      std::cout << PHWHERE << "TrkrClusterSourceLinks node not found on node tree. Exiting."
		<< std::endl;
      
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  
  m_tGeometry = findNode::getClass<ActsTrackingGeometry>(topNode, "ActsTrackingGeometry");

  if (!m_tGeometry)
  {
    std::cout << "ActsGeometry not on node tree. Exiting."
              << std::endl;

    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_trackMap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if (!m_trackMap)
    {
      std::cout << PHWHERE << " ERROR: Can't find SvtxTrackMap. Exiting " 
		<< std::endl;

      return Fun4AllReturnCodes::ABORTEVENT;
    }


  return Fun4AllReturnCodes::EVENT_OK;
}

