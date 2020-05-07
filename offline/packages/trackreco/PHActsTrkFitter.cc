/*!
 *  \file		PHActsTrkFitter.C
 *  \brief		Refit SvtxTracks with PHActs.
 *  \details	Refit SvtxTracks with PHActs.
 *  \author	        Tony Frawley <afrawley@fsu.edu>
 */

#include "PHActsTrkFitter.h"
#include "MakeActsGeometry.h"
#include "ActsTrack.h"

/// Tracking includes
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <Acts/EventData/TrackParameters.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>
#include <Acts/Surfaces/PlaneSurface.hpp>
#include <Acts/Surfaces/Surface.hpp>

#include <ACTFW/EventData/Track.hpp>
#include <ACTFW/Framework/AlgorithmContext.hpp>

#include <cmath>
#include <iostream>
#include <vector>

PHActsTrkFitter::PHActsTrkFitter(const std::string& name)
  : PHTrackFitting(name)
  , m_event(0)
  , m_actsProtoTracks(nullptr)
  , m_tGeometry(nullptr)
{
  Verbosity(0);
}

PHActsTrkFitter::~PHActsTrkFitter()
{
}

int PHActsTrkFitter::Setup(PHCompositeNode* topNode)
{
  if (getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
    return Fun4AllReturnCodes::ABORTEVENT;
  
  fitCfg.fit = FW::TrkrClusterFittingAlgorithm::makeFitterFunction(
               m_tGeometry->tGeometry,
	       m_tGeometry->magField,
	       Acts::Logging::INFO);

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsTrkFitter::Process()
{

  m_event++;

  if (Verbosity() > 1)
  {
    std::cout << PHWHERE << "Events processed: " << m_event << std::endl;
    std::cout << "Start PHActsTrkFitter::process_event" << std::endl;
  }


  /// Construct a perigee surface as the target surface (?)
  auto pSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(
	          Acts::Vector3D{0., 0., 0.});

  std::map<unsigned int, ActsTrack>::iterator trackIter;

  for (trackIter = m_actsProtoTracks->begin();
       trackIter != m_actsProtoTracks->end();
       ++trackIter)
  {
    ActsTrack track = trackIter->second;
    /// Can correlate with the SvtxTrackMap with the key
    const unsigned int trackKey = trackIter->first;

    std::vector<SourceLink> sourceLinks = track.getSourceLinks();
    FW::TrackParameters trackSeed = track.getTrackParams();
      
    if(Verbosity() > 10)
      {
	std::cout << " Processing proto track with position:" 
		  << trackSeed.position() << std::endl 
		  << "momentum: " << trackSeed.momentum() << std::endl;

      }
    /// Call KF now. Have a vector of sourceLinks corresponding to clusters
    /// associated to this track and the corresponding track seed which
    /// corresponds to the PHGenFitTrkProp track seeds
    Acts::KalmanFitterOptions<Acts::VoidOutlierFinder> kfOptions(
      m_tGeometry->geoContext,
      m_tGeometry->magFieldContext,
      m_tGeometry->calibContext,
      Acts::VoidOutlierFinder(),
      &(*pSurface));
  
    auto result = fitCfg.fit(sourceLinks, trackSeed, kfOptions);

    /// Check that the result is okay
    if (result.ok())
    {  
      const Acts::KalmanFitterResult<SourceLink>& fitOutput = result.value();
      if (fitOutput.fittedParameters)
      {
        updateSvtxTrack(fitOutput, trackKey);

        /// Get position, momentum from params
        if (Verbosity() > 10)
        {
	  const auto& params = fitOutput.fittedParameters.value();
          std::cout << "Fitted parameters for track" << std::endl;
          std::cout << " position : " << params.position().transpose()
                    << std::endl;
          std::cout << " momentum : " << params.momentum().transpose()
                    << std::endl;
        }
      }
    }
    

  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsTrkFitter::End(PHCompositeNode* topNode)
{
  if (Verbosity() > 10)
  {
    std::cout << "Finished PHActsTrkFitter" << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

void PHActsTrkFitter::updateSvtxTrack(const Acts::KalmanFitterResult<SourceLink>& fitOutput, const unsigned int trackKey)
{
  const auto& params = fitOutput.fittedParameters.value();
 
  SvtxTrackMap::Iter trackIter = m_trackMap->find(trackKey);
  SvtxTrack *track = trackIter->second;

  /// Acts default unit is mm. So convert to cm
  track->set_x(params.position()(0) / Acts::UnitConstants::cm);
  track->set_y(params.position()(1) / Acts::UnitConstants::cm);
  track->set_z(params.position()(2) / Acts::UnitConstants::cm);
  track->set_px(params.momentum()(0));
  track->set_py(params.momentum()(1));
  track->set_pz(params.momentum()(2));
  

  if(params.covariance())
    {
   
      Acts::BoundSymMatrix rotatedCov = rotateCovarianceLocalToGlobal(fitOutput);
      for(int i = 0; i < 6; i++)
	for(int j = 0; j < 6; j++)
	  track->set_error(i,j, rotatedCov(i,j));
    
    }
 
  return;

}

Acts::BoundSymMatrix PHActsTrkFitter::rotateCovarianceLocalToGlobal(
                     const Acts::KalmanFitterResult<SourceLink>& fitOutput)
{
  Acts::BoundSymMatrix matrix = Acts::BoundSymMatrix::Zero();
  
  const auto& params = fitOutput.fittedParameters.value();
  auto covarianceMatrix = *params.covariance();

  const double px = params.momentum()(0);
  const double py = params.momentum()(1);
  const double pz = params.momentum()(2);
  const double p = sqrt(px * px + py * py + pz * pz);
  
  const double x = params.position()(0);
  const double y = params.position()(1);

  const int charge = params.charge();
  const double phiPos = atan2(x, y);

  /// We need to rotate the opposite of what was done in PHActsTracks.
  /// So first rotate from (x_l, y_l, phi, theta, q/p, time) to 
  /// (x_l, y_l, phi, theta, p, time)

  Acts::BoundSymMatrix qprotation = Acts::BoundSymMatrix::Zero();
  qprotation(0,0) = 1;
  qprotation(1,1) = 1;
  qprotation(2,2) = 1;
  qprotation(3,3) = 1;
  qprotation(4,4) = charge * charge / (p * p * p * p);
  qprotation(5,5) = 1;
  
  /// Want the inverse of the rotation matrix from PHActsTracks 
  /// because we are rotating back from local to global. So we do R^TCR 
  /// rather than RCR^T
  matrix = qprotation.transpose() * covarianceMatrix * qprotation;

  /// Now rotate to (x_g, y_g, px, py, pz, t)
  /// Make a unit p vector for the rotation
  const double uPx = px / p;
  const double uPy = py / p;
  const double uPz = pz / p;
  const double uP = sqrt(uPx * uPx + uPy * uPy + uPz * uPz);
  
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
  /// Undoing the rotation, so R^TCR instead of RCR^T
  matrix = rotation.transpose() * matrix * rotation;
  
  if(Verbosity() > 20)
    {
      std::cout<<"Rotated Matrix"<<std::endl;
      for(int i =0;i<6; i++)
	{
	  for(int j =0; j<6; j++)
	    std::cout<<matrix(i,j)<<", ";
	  std::cout<<std::endl;
	}
    }

  /// Now matrix is in basis (x_g, y_g, px, py, pz, t)
  /// Shift rows and columns to get like (x_g, y_g, z, px, py, pz)
  Acts::BoundSymMatrix svtxCovariance = Acts::BoundSymMatrix::Zero();
  for(int i = 0; i < 6; i++ )
    {
      for(int j = 0; j < 6; j++ )
	{
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
	  
	  svtxCovariance(row,col) = matrix(i, j);
	  
	  if(row < 2 && col < 2)
	    svtxCovariance(row,col) /= Acts::UnitConstants::cm2;
	  else if(row < 2 && col < 5)
	    svtxCovariance(row,col) /= Acts::UnitConstants::cm;
	  else if (row < 5 && col < 2)
	    svtxCovariance(row,col) /= Acts::UnitConstants::cm;
	}
    }
  if(Verbosity() > 20)
    {
      std::cout << "shiftd matrix"<<std::endl;
      for(int i =0;i<6;i++)
	{
	  for(int j =0;j<6; j++)
	    std::cout<<svtxCovariance(i,j)<<", ";
	  std::cout<<std::endl;
	}
    }

  return svtxCovariance;
}

int PHActsTrkFitter::createNodes(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsTrkFitter::getNodes(PHCompositeNode* topNode)
{
  m_actsProtoTracks = findNode::getClass<std::map<unsigned int, ActsTrack>>(topNode, "ActsTrackMap");

  if (!m_actsProtoTracks)
  {
    std::cout << "Acts proto tracks not on node tree. Exiting."
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_tGeometry = findNode::getClass<ActsTrackingGeometry>(topNode, "ActsTrackingGeometry");
  if(!m_tGeometry)
    {
      std::cout << "ActsTrackingGeometry not on node tree. Exiting."
		<< std::endl;
      
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  m_trackMap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  
  if(!m_trackMap)
    {
      std::cout << PHWHERE << "SvtxTrackMap not found on node tree. Exiting."
		<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  return Fun4AllReturnCodes::EVENT_OK;
}
