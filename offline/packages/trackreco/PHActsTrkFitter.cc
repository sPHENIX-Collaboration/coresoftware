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
  if(m_event % 10 == 0)
    std::cout << "Processed " << m_event << " events " << std::endl;
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
	SvtxTrack *track = convertActsToSvtx(fitOutput, trackKey);
	if(!track)
	  std::cout<<std::endl;
        const auto& params = fitOutput.fittedParameters.value();
        /// Get position, momentum from params
        if (Verbosity() > 10)
        {
          std::cout << "Fitted parameters for track" << std::endl;
          std::cout << " position : " << params.position().transpose()
                    << std::endl;
          std::cout << " momentum : " << params.momentum().transpose()
                    << std::endl;
        }
      }
    }

    /// Update the acts track node on the node tree
    


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

SvtxTrack* PHActsTrkFitter::convertActsToSvtx(const Acts::KalmanFitterResult<SourceLink>& fitOutput, const unsigned int trackKey)
{

  return nullptr;
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

  return Fun4AllReturnCodes::EVENT_OK;
}
