#include "PHActsSiliconMicromegasFitter.h"
#include "MakeActsGeometry.h"


/// Tracking includes
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackState_v1.h>
#include <trackbase_historic/SvtxTrackMap.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/PHDataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>

#include <Acts/EventData/TrackParameters.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>
#include <Acts/Surfaces/PlaneSurface.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/EventData/MultiTrajectory.hpp>
#include <Acts/EventData/MultiTrajectoryHelpers.hpp>

#include <ActsExamples/EventData/Track.hpp>
#include <ActsExamples/Framework/AlgorithmContext.hpp>

PHActsSiliconMicromegasFitter::PHActsSiliconMicromegasFitter(const std::string &name) :
  PHTrackFitting(name)
  , m_event(0)
  , m_trackMap(nullptr)
  , m_hitIdClusKey(nullptr)
  , m_actsProtoTracks(nullptr)
  , m_tGeometry(nullptr)
 
{

}


int PHActsSiliconMicromegasFitter::Setup(PHCompositeNode *topNode)
{
  if(Verbosity() > 0)
    std::cout << "Setup PHActsSiliconMicromegasFitter" << std::endl;
  
  if(getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
    return Fun4AllReturnCodes::ABORTEVENT;

  m_fitCfg.fit = ActsExamples::TrkrClusterFittingAlgorithm::makeFitterFunction(
        m_tGeometry->tGeometry,  
	m_tGeometry->magField);

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsSiliconMicromegasFitter::Process()
{
  if(Verbosity() > 0)
    {
      std::cout << "Starting PHActsSiliconMicromegasFitter at event " 
		<< m_event << std::endl;
    }
  
  /// Set the Acts fitting logging level
  auto logLevel = Acts::Logging::INFO;
  if(Verbosity() > 3)
    logLevel = Acts::Logging::VERBOSE;
  
  auto logger = Acts::getDefaultLogger("PHActsSiliconMicromegasFitter",
				       logLevel);
  
  for(auto trackIter : *m_actsProtoTracks)
    {
      auto track = trackIter.second;
      auto sourceLinks = track.getSourceLinks();
      auto trackSeed = track.getTrackParams();

      /// Get the surfaces and SLs for the silicon+MMs
      SurfaceVec surfaces;
      auto siliconMMsSls = getSiliconMMsSls(sourceLinks, surfaces);

      /// Reset track covariance matrix to something Acts can work with
      Acts::BoundSymMatrix cov;
      cov << 1000 * Acts::UnitConstants::um, 0., 0., 0., 0., 0.,
           0., 1000 * Acts::UnitConstants::um, 0., 0., 0., 0.,
           0., 0., 0.05, 0., 0., 0.,
           0., 0., 0., 0.05, 0., 0.,
           0., 0., 0., 0., 0.00005 , 0.,
           0., 0., 0., 0., 0., 1.;
    
      ActsExamples::TrackParameters newTrackSeed(
		   trackSeed.fourPosition(m_tGeometry->geoContext),
		   trackSeed.momentum(),
		   trackSeed.absoluteMomentum(),
		   trackSeed.charge(),
		   cov);

      auto pSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(
        	      track.getVertex());

      /// Call direct navigation Acts KF
      Acts::KalmanFitterOptions<Acts::VoidOutlierFinder> kfOptions(
            m_tGeometry->geoContext,
	    m_tGeometry->magFieldContext,
	    m_tGeometry->calibContext,
	    Acts::VoidOutlierFinder(),
	    Acts::LoggerWrapper(*logger),
	    Acts::PropagatorPlainOptions(),
	    &(*pSurface));
      
      auto result = m_fitCfg.fit(siliconMMsSls, newTrackSeed,
      			 kfOptions, surfaces);

      if(result.ok())
	{
	  auto fitOutput = result.value();
	}
      
    }

  if(Verbosity() > 0)
    {
      std::cout << "Finish PHActsSiliconMicromegasFitter event "
		<< m_event << std::endl;
    }

  return Fun4AllReturnCodes::EVENT_OK;

}

int PHActsSiliconMicromegasFitter::ResetEvent(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsSiliconMicromegasFitter::End(PHCompositeNode *topNode)
{

  return Fun4AllReturnCodes::EVENT_OK;

}

SourceLinkVec PHActsSiliconMicromegasFitter::getSiliconMMsSls(SourceLinkVec trackSls, 
							      SurfaceVec &surfaces)
{

  SourceLinkVec siliconMMSls;

  for(auto sl : trackSls)
    {
      auto volume = sl.referenceSurface().geometryId().volume();

      /// If volume is not the TPC add it to the list
      if(volume != 14)
	{
	  siliconMMSls.push_back(sl);
	  surfaces.push_back(&sl.referenceSurface());
	}
    }

  return siliconMMSls;
}

int PHActsSiliconMicromegasFitter::getNodes(PHCompositeNode *topNode)
{
  m_tGeometry = findNode::getClass<ActsTrackingGeometry>(topNode, 
							 "ActsTrackingGeometry");
  
  if(!m_tGeometry)
    {
      std::cout << PHWHERE << " ActsTrackingGeometry not on node tree, exiting" 
		<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  m_actsProtoTracks = findNode::getClass<std::map<unsigned int, ActsTrack>>
      (topNode, "ActsTrackMap");

  if(!m_actsProtoTracks)
    {
      std::cout << PHWHERE << " Acts proto tracks not on node tree, exiting." 
		<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;

    }

  m_trackMap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if(!m_trackMap)
    {
      std::cout << PHWHERE << " SvtxTrackMap not on node tree, exiting." 
		<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  m_hitIdClusKey = findNode::getClass<std::map<TrkrDefs::cluskey, unsigned int>>(topNode, "HitIDClusIDActsMap");
  if(!m_hitIdClusKey)
    {
      std::cout << PHWHERE << " hit id cluskey Acts map not on node tree, exiting." 
		<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  return Fun4AllReturnCodes::EVENT_OK;
}
