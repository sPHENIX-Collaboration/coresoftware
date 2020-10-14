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

  outfile = new TFile("siliconMMoutfile.root","recreate");
  h_nClus = new TH2I("h_nClus",";N_{clus}^{in}; N_{clus}^{out}",
		     10,0,10,10,0,10);

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

      int nInSurf = surfaces.size();
      
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
	      //const auto& params = fitOutput.fittedParameters.value();
	      
	    }
	  
	  Trajectory trajectory(fitOutput.fittedStates, trackTips, indexedParams);
	  int nOutSurf = checkClusterKeys(trajectory,fitOutput.trackTip);
	 
	  h_nClus->Fill(nInSurf, nOutSurf);
	}
      
    }

  if(Verbosity() > 0)
    {
      std::cout << "Finish PHActsSiliconMicromegasFitter event "
		<< m_event << std::endl;
    }

  return Fun4AllReturnCodes::EVENT_OK;

}

int PHActsSiliconMicromegasFitter::checkClusterKeys(Trajectory traj, 
						     const size_t trackTip)
{
  int nOutSurf = 0;
  const auto &[trackTips, mj] = traj.trajectory();
 
  mj.visitBackwards(trackTip, [&](const auto &state) {
      /// Only fill the track states with non-outlier measurement
      auto typeFlags = state.typeFlags();
      if (not typeFlags.test(Acts::TrackStateFlag::MeasurementFlag))
	{
	  return true;
	}
      
      auto meas = std::get<Measurement>(*state.uncalibrated());

      nOutSurf++;

      if(Verbosity() > 0)
	std::cout << "obtained surface with geoID : "
		  << meas.referenceObject().geometryId() << std::endl;
      /// Get local position
      Acts::Vector2D local(meas.parameters()[Acts::eBoundLoc0],
			   meas.parameters()[Acts::eBoundLoc1]);

      /// This is an arbitrary vector. Doesn't matter in coordinate transformation
      /// in Acts code
      Acts::Vector3D mom(1., 1., 1.);
      Acts::Vector3D global = meas.referenceObject().localToGlobal(
					    m_tGeometry->geoContext,
					    local, mom);
   
      return true;
    });

  return nOutSurf;
}


int PHActsSiliconMicromegasFitter::ResetEvent(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsSiliconMicromegasFitter::End(PHCompositeNode *topNode)
{
  outfile->cd();
  h_nClus->Write();
  outfile->Close();
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
