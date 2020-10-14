#include "PHActsSiliconMicromegasFitter.h"
#include "MakeActsGeometry.h"
#include "ActsTransformations.h"

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

#include <TH2.h>
#include <TFile.h>

PHActsSiliconMicromegasFitter::PHActsSiliconMicromegasFitter(
    const std::string &name) : PHTrackFitting(name)
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

  if(createNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
    return Fun4AllReturnCodes::ABORTEVENT;

  /// Call the directed navigator fitter function
  m_fitCfg.dFit = ActsExamples::TrkrClusterFittingAlgorithm::makeFitterFunction(
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
      auto trackKey = trackIter.first;
      auto track = trackIter.second;
      auto sourceLinks = track.getSourceLinks();
      auto trackSeed = track.getTrackParams();

      /// Get the surfaces and SLs for the silicon+MMs
      SurfaceVec surfaces;
      auto siliconMMsSls = getSiliconMMsSls(sourceLinks, surfaces);

      /// If no silicon+MM surfaces in this track, continue to next track
      if(surfaces.size() == 0)
	continue;

      /// We only care about fitting tracks with a MM surface
      bool MMsurface = false;
      for(auto surf : surfaces)
	{
	  if(surf->geometryId().volume() == 16)
	    MMsurface = true;
	}
      if(!MMsurface)
	continue;

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
      
      auto result = m_fitCfg.dFit(siliconMMsSls, newTrackSeed,
				  kfOptions, surfaces);

      if(result.ok())
	{
	  auto fitOutput = result.value();
	  /// Make a trajectory state for storage, which conforms to Acts track fit
	  /// analysis tool
	  std::vector<size_t> trackTips;
	  trackTips.push_back(fitOutput.trackTip);
	  ActsExamples::IndexedParams indexedParams;
	  if (fitOutput.fittedParameters)
	    {
	      indexedParams.emplace(fitOutput.trackTip, 
				    fitOutput.fittedParameters.value());
	    }
	  
	  Trajectory trajectory(fitOutput.fittedStates, trackTips, indexedParams);
	  
	  if(fitOutput.fittedParameters)
	    updateSvtxTrack(trajectory, trackKey, track.getVertex());
	}
      
    }

  if(Verbosity() > 0)
    {
      std::cout << "Finish PHActsSiliconMicromegasFitter event "
		<< m_event << std::endl;
    }

  return Fun4AllReturnCodes::EVENT_OK;

}

void PHActsSiliconMicromegasFitter::updateSvtxTrack(Trajectory traj, 
						    const unsigned int trackKey,
						    Acts::Vector3D vertex)
{
  const auto &[trackTips, mj] = traj.trajectory();
  /// only one track tip in the track fit Trajectory
  auto &trackTip = trackTips.front();

  SvtxTrackMap::Iter trackIter = m_trackMap->find(trackKey);
  SvtxTrack *track = trackIter->second;
  
  if(Verbosity() > 2)
    {
      std::cout << "Identify (proto) track before updating with acts results " << std::endl;
      track->identify();
      std::cout << " cluster keys size " << track->size_cluster_keys() << std::endl;  
    }

  // The number of associated clusters may have changed - start over
  track->clear_states();
  track->clear_cluster_keys();

  // create a state at pathlength = 0.0
  // This state holds the track parameters, which will be updated below
  float pathlength = 0.0;
  SvtxTrackState_v1 out(pathlength);
  out.set_x(0.0);
  out.set_y(0.0);
  out.set_z(0.0);
  track->insert_state(&out);   

  auto trajState =
    Acts::MultiTrajectoryHelpers::trajectoryState(mj, trackTip);
 
  const auto& params = traj.trackParameters(trackTip);

  /// Acts default unit is mm. So convert to cm
  track->set_x(params.position(m_tGeometry->geoContext)(0)
	       / Acts::UnitConstants::cm);
  track->set_y(params.position(m_tGeometry->geoContext)(1)
	       / Acts::UnitConstants::cm);
  track->set_z(params.position(m_tGeometry->geoContext)(2)
	       / Acts::UnitConstants::cm);

  track->set_px(params.momentum()(0));
  track->set_py(params.momentum()(1));
  track->set_pz(params.momentum()(2));
  
  track->set_charge(params.charge());
  track->set_chisq(trajState.chi2Sum);
  track->set_ndf(trajState.NDF);

  ActsTransformations *rotater = new ActsTransformations();
  rotater->setVerbosity(Verbosity());
  
  if(params.covariance())
    {
   
      Acts::BoundSymMatrix rotatedCov = 
	rotater->rotateActsCovToSvtxTrack(params,
					  m_tGeometry->geoContext);
      
      for(int i = 0; i < 6; i++)
	{
	  for(int j = 0; j < 6; j++)
	    {
	      track->set_error(i,j, rotatedCov(i,j));
	    }
	}
    }
 
  float dca3Dxy = -9999.;
  float dca3Dz = -9999.;
  float dca3DxyCov = -9999.;
  float dca3DzCov = -9999.;

  rotater->calculateDCA(params, vertex, m_tGeometry->geoContext, 
			dca3Dxy, dca3Dz, dca3DxyCov, dca3DzCov);

  // convert from mm to cm
  track->set_dca3d_xy(dca3Dxy / Acts::UnitConstants::cm);
  track->set_dca3d_z(dca3Dz / Acts::UnitConstants::cm);
  track->set_dca3d_xy_error(dca3DxyCov / Acts::UnitConstants::cm);
  track->set_dca3d_z_error(dca3DzCov / Acts::UnitConstants::cm);
  
  // Also need to update the state list and cluster ID list for all measurements associated with the acts track  
  // loop over acts track states, copy over to SvtxTrackStates, and add to SvtxTrack

  rotater->fillSvtxTrackStates(traj, trackTip, track,
			       m_tGeometry->geoContext,
			       m_hitIdClusKey);  

  if(Verbosity() > 2)
    {  
      std::cout << " Identify fitted track after updating track states:" << std::endl;
      track->identify();
      std::cout << " cluster keys size " << track->size_cluster_keys() << std::endl;  
    }
 
 return;
  
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

  if(Verbosity() > 0)
    std::cout << "Sorting " << trackSls.size() << " SLs" << std::endl;
  
  for(auto sl : trackSls)
    {
      auto volume = sl.referenceSurface().geometryId().volume();
      
      /// If volume is not the TPC add it to the list
      if(volume != 14)
	{
	  siliconMMSls.push_back(sl);
	  surfaces.push_back(&sl.referenceSurface());
	  if(Verbosity() > 0)
	    std::cout << "Adding surface to sequence with geoID : "
		      << sl.referenceSurface().geometryId() << std::endl;
	}
    }

  /// Surfaces need to be sorted in order, i.e. from smallest to
  /// largest radius extending from target surface
  /// Add a check to ensure this
  if(surfaces.size() > 0)
    checkSurfaceVec(surfaces);

  return siliconMMSls;
}

void PHActsSiliconMicromegasFitter::checkSurfaceVec(SurfaceVec &surfaces)
{
  for(int i = 0; i < surfaces.size() - 1; i++)
    {
      auto surface = surfaces.at(i);
      auto thisVolume = surface->geometryId().volume();
      auto thisLayer  = surface->geometryId().layer();
      
      auto nextSurface = surfaces.at(i+1);
      auto nextVolume = nextSurface->geometryId().volume();
      auto nextLayer = nextSurface->geometryId().layer();
      
      /// Implement a check to ensure surfaces are sorted
      if(nextVolume == thisVolume) 
	{
	  if(nextLayer < thisLayer)
	    {
	      std::cout << PHWHERE 
			<< "Surface not in order... removing surface" 
			<< std::endl;
	      surfaces.erase(surfaces.begin() + i);
	      /// Subtract one so we don't skip a surface
	      i--;
	      continue;
	    }
	}
      else 
	{
	  if(nextVolume < thisVolume)
	    {
	      std::cout << PHWHERE 
			<< "Volume not in order... removing surface" 
			<< std::endl;
	      surfaces.erase(surfaces.begin() + i);
	      /// Subtract one so we don't skip a surface
	      i--;
	      continue;
	    }
	}
    } 
}

int PHActsSiliconMicromegasFitter::createNodes(PHCompositeNode *topNode)
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

  m_actsFitResults = findNode::getClass<std::map<const unsigned int, Trajectory>>(topNode, "ActsFitResults");
  
  if(!m_actsFitResults)
    {
      m_actsFitResults = new std::map<const unsigned int, Trajectory>;

      PHDataNode<std::map<const unsigned int, 
			  Trajectory>> *fitNode = 
		 new PHDataNode<std::map<const unsigned int, 
				    Trajectory>>
		 (m_actsFitResults, "ActsFitResults");

      svtxNode->addNode(fitNode);
      
    }

  return Fun4AllReturnCodes::EVENT_OK;

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
