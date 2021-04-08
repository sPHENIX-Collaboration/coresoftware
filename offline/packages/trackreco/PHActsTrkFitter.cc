/*!
 *  \file		PHActsTrkFitter.C
 *  \brief		Refit SvtxTracks with PHActs.
 *  \details	Refit SvtxTracks with PHActs.
 *  \author	        Tony Frawley <afrawley@fsu.edu>
 */

#include "PHActsTrkFitter.h"
#include "ActsTrack.h"
#include "ActsTransformations.h"

/// Tracking includes
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackState_v1.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxVertexMap.h>
#include <trackbase_historic/SvtxVertex.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/PHDataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/PHTimer.h>

#include <Acts/EventData/TrackParameters.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>
#include <Acts/Surfaces/PlaneSurface.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/EventData/MultiTrajectory.hpp>
#include <Acts/EventData/MultiTrajectoryHelpers.hpp>

#include <ActsExamples/EventData/Track.hpp>
#include <ActsExamples/Framework/AlgorithmContext.hpp>

#include <cmath>
#include <iostream>
#include <vector>

PHActsTrkFitter::PHActsTrkFitter(const std::string& name)
  : PHTrackFitting(name)
  , m_event(0)
  , m_tGeometry(nullptr)
  , m_trackMap(nullptr)
  , m_vertexMap(nullptr)
  , m_clusterContainer(nullptr)
  , m_surfMaps(nullptr)
  , m_nBadFits(0)
  , m_fitSiliconMMs(false)
  , m_fillSvtxTrackStates(true)
  , m_timeAnalysis(false)
  , m_timeFile(nullptr)
  , h_eventTime(nullptr)
  , h_fitTime(nullptr)
  , h_updateTime(nullptr)
  , h_stateTime(nullptr)
  , h_rotTime(nullptr)
{
  Verbosity(0);
}

PHActsTrkFitter::~PHActsTrkFitter()
{
}

int PHActsTrkFitter::Setup(PHCompositeNode* topNode)
{
  if(Verbosity() > 1)
    std::cout << "Setup PHActsTrkFitter" << std::endl;
  
  if(createNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
    return Fun4AllReturnCodes::ABORTEVENT;

  if (getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
    return Fun4AllReturnCodes::ABORTEVENT;

  m_fitCfg.fit = ActsExamples::TrkrClusterFittingAlgorithm::makeFitterFunction(
               m_tGeometry->tGeometry,
	       m_tGeometry->magField);

  m_fitCfg.dFit = ActsExamples::TrkrClusterFittingAlgorithm::makeFitterFunction(
	       m_tGeometry->magField);

  if(m_timeAnalysis)
    {
      m_timeFile = new TFile(std::string(Name() + ".root").c_str(), 
			     "RECREATE");
      h_eventTime = new TH1F("h_eventTime", ";time [ms]",
			     100000, 0, 10000);
      h_fitTime = new TH2F("h_fitTime",";p_{T} [GeV];time [ms]",
			   80, 0, 40, 100000, 0, 1000);
      h_updateTime = new TH1F("h_updateTime",";time [ms]",
			      100000, 0, 1000);
      
      h_rotTime = new TH1F("h_rotTime", ";time [ms]",
			   100000, 0, 1000);
      h_stateTime = new TH1F("h_stateTime", ";time [ms]",
			     100000, 0, 1000);			     
    }		 
  
  if(Verbosity() > 1)
    std::cout << "Finish PHActsTrkFitter Setup" << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsTrkFitter::Process()
{
  auto eventTimer = std::make_unique<PHTimer>("eventTimer");
  eventTimer->stop();
  eventTimer->restart();
  
  m_event++;

  auto logLevel = Acts::Logging::FATAL;

  if (Verbosity() > 1)
  {
    std::cout << PHWHERE << "Events processed: " << m_event << std::endl;
    std::cout << "Start PHActsTrkFitter::process_event" << std::endl;
    if(Verbosity() > 4)
      logLevel = Acts::Logging::VERBOSE;
  }

  loopTracks(logLevel);

  eventTimer->stop();
  auto eventTime = eventTimer->get_accumulated_time();

  if(Verbosity() > 0)
    std::cout << "PHActsTrkFitter total event time " 
	      << eventTime << std::endl;

  if(m_timeAnalysis)     
    h_eventTime->Fill(eventTime);
    

  if(Verbosity() > 1)
    std::cout << "PHActsTrkFitter::process_event finished" 
	      << std::endl;

  // put this in the output file
  if(Verbosity() > 0)
    std::cout << " SvtxTrackMap size is now " << m_trackMap->size() 
	      << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsTrkFitter::ResetEvent(PHCompositeNode *topNode)
{

  if(Verbosity() > 1)
    {
      std::cout << "Reset PHActsTrkFitter" << std::endl;

    }
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsTrkFitter::End(PHCompositeNode *topNode)
{
  if(m_timeAnalysis)
    {
      m_timeFile->cd();
      h_fitTime->Write();
      h_eventTime->Write();
      h_rotTime->Write();
      h_stateTime->Write();
      h_updateTime->Write();
      m_timeFile->Write();
      m_timeFile->Close();
    } 

  if (Verbosity() > 0)
  {
    std::cout << "The Acts track fitter had " << m_nBadFits 
	      << " fits return an error" << std::endl;

    std::cout << "Finished PHActsTrkFitter" << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

void PHActsTrkFitter::loopTracks(Acts::Logging::Level logLevel)
{
  auto logger = Acts::getDefaultLogger("PHActsTrkFitter", logLevel);

  /// Store a vector of track fits that fail to erase, so that the
  /// track map iterator doesn't crash
  std::vector<unsigned int> badTracks;

  ActsTransformations transformer;

  for(auto& [trackKey, track] : *m_trackMap)
    {
  
      if(!track)
	{
	  continue;
	}

      auto trackTimer = std::make_unique<PHTimer>("TrackTimer");
      trackTimer->stop();
      trackTimer->restart();

      auto sourceLinks = getSourceLinks(track);
      /// If using directed navigation, collect surface list to navigate
   
      SurfacePtrVec surfaces;
      if(m_fitSiliconMMs)
	{	
	  sourceLinks = getSurfaceVector(sourceLinks, surfaces);
	  /// Check to see if there is a track to fit, if not skip it
	  if(surfaces.size() == 0)
	    continue;
	  bool MMsurface = false;
	  for(auto surf : surfaces)
	    {
	      if(surf->geometryId().volume() == 16)
		{
		  MMsurface = true;
		  break;
		}
	    }
	  /// If there's not a MM surface, we don't want to fit only
	  /// the silicon
	  if(!MMsurface)
	    continue;
	}
  
      Acts::Vector3D momentum(track->get_px(), 
			      track->get_py(), 
			      track->get_pz());
   
      auto actsVertex = getVertex(track);
      auto pSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(
					  actsVertex);
      auto actsFourPos = Acts::Vector4D(actsVertex(0), actsVertex(1),
					actsVertex(2),
					10 * Acts::UnitConstants::ns);
      Acts::BoundSymMatrix cov = setDefaultCovariance();
 
      /// Reset the track seed with the dummy covariance and the 
      /// primary vertex as the track position
      ActsExamples::TrackParameters seed(actsFourPos,
					 momentum,
					 track->get_p(),
					 track->get_charge(),
					 cov);

      if(Verbosity() > 2)
	printTrackSeed(seed);
     
      /// Call KF now. Have a vector of sourceLinks corresponding to clusters
      /// associated to this track and the corresponding track seed which
      /// corresponds to the PHGenFitTrkProp track seeds
      Acts::KalmanFitterOptions<Acts::VoidOutlierFinder> kfOptions(
			        m_tGeometry->geoContext,
				m_tGeometry->magFieldContext,
				m_tGeometry->calibContext,
				Acts::VoidOutlierFinder(),
				Acts::LoggerWrapper(*logger),
				Acts::PropagatorPlainOptions(),
				&(*pSurface));
 
      auto fitTimer = std::make_unique<PHTimer>("FitTimer");
      fitTimer->stop();
      fitTimer->restart();
      auto result = fitTrack(sourceLinks, seed, kfOptions,
			     surfaces);
      fitTimer->stop();
      auto fitTime = fitTimer->get_accumulated_time();
      
      if(Verbosity() > 0)
	std::cout << "PHActsTrkFitter Acts fit time "
		  << fitTime << std::endl;

      /// Check that the track fit result did not return an error
      if (result.ok())
	{  
	  const FitResult& fitOutput = result.value();
	  
	  if(m_timeAnalysis)
	    {
	      h_fitTime->Fill(fitOutput.fittedParameters.value()
			      .transverseMomentum(), 
			      fitTime);
	    }

	  getTrackFitResult(fitOutput, track);

	}
      else
	{
	  /// Track fit failed, get rid of the track from the map
	  badTracks.push_back(trackKey);
	  m_nBadFits++;

	}

      trackTimer->stop();
      auto trackTime = trackTimer->get_accumulated_time();
      
      if(Verbosity() > 0)
	std::cout << "PHActsTrkFitter total single track time "
		  << trackTime << std::endl;
    
    }

  /// Now erase bad tracks from the track map
  for(auto key : badTracks)
    {
      m_trackMap->erase(key);
    }

  return;

}
void PHActsTrkFitter::printTrackSeed(ActsExamples::TrackParameters seed)
{
  std::cout << PHWHERE << " Processing proto track with position:" 
	    << seed.position(m_tGeometry->geoContext) 
	    << std::endl 
	    << "momentum: " << seed.momentum() 
	    << std::endl
	    << "charge : " << seed.charge() 
	    << std::endl;
  std::cout << "proto track covariance " << std::endl
	    << seed.covariance().value() << std::endl;
  
}
Acts::Vector3D PHActsTrkFitter::getVertex(SvtxTrack *track)
{
  auto vertexId = track->get_vertex_id();
  const SvtxVertex* svtxVertex = m_vertexMap->get(vertexId);
  Acts::Vector3D vertex(svtxVertex->get_x() * Acts::UnitConstants::cm, 
			svtxVertex->get_y() * Acts::UnitConstants::cm, 
			svtxVertex->get_z() * Acts::UnitConstants::cm);
  return vertex;
}
Surface PHActsTrkFitter::getSurface(TrkrDefs::cluskey cluskey, 
				    TrkrDefs::subsurfkey surfkey)
{

  auto trkrid = TrkrDefs::getTrkrId(cluskey);
  auto hitsetkey = TrkrDefs::getHitSetKeyFromClusKey(cluskey);

  if(trkrid == TrkrDefs::TrkrId::mvtxId or
     trkrid == TrkrDefs::TrkrId::inttId)
    {
      return getSiliconSurface(hitsetkey);
    }
  else
    {
      return getTpcMMSurface(hitsetkey, surfkey);
    }

}
Surface PHActsTrkFitter::getSiliconSurface(TrkrDefs::hitsetkey hitsetkey)
{
  auto surfMap = m_surfMaps->siliconSurfaceMap;
  auto iter = surfMap.find(hitsetkey);
  if(iter != surfMap.end())
    {
      return iter->second;
    }
  
  /// If it can't be found, return nullptr
  return nullptr;

}
Surface PHActsTrkFitter::getTpcMMSurface(TrkrDefs::hitsetkey hitsetkey,
					 TrkrDefs::subsurfkey surfkey)
{
  auto surfMap = m_surfMaps->tpcSurfaceMap;
  auto iter = surfMap.find(hitsetkey);
  if(iter != surfMap.end())
    {
      auto surfvec = iter->second;
      return surfvec.at(surfkey);
    }
  
  /// If it can't be found, return nullptr to skip this cluster
  return nullptr;

}
SourceLinkVec PHActsTrkFitter::getSourceLinks(SvtxTrack* track)
{

  SourceLinkVec sourcelinks;

  for (SvtxTrack::ConstClusterKeyIter clusIter = track->begin_cluster_keys();
       clusIter != track->end_cluster_keys();
       ++clusIter)
    {
      auto key = *clusIter;
      auto cluster = m_clusterContainer->findCluster(key);

      auto subsurfkey = cluster->getSubSurfKey();
      
      /// Make a safety check for clusters that couldn't be attached
      /// to a surface
      auto surf = getSurface(key, subsurfkey);
      if(!surf)
	continue;
  
      Acts::BoundVector loc = Acts::BoundVector::Zero();
      loc[Acts::eBoundLoc0] = cluster->getLocalX() * Acts::UnitConstants::cm;
      loc[Acts::eBoundLoc1] = cluster->getLocalY() * Acts::UnitConstants::cm;
      
      Acts::BoundMatrix cov = Acts::BoundMatrix::Zero();
      cov(Acts::eBoundLoc0, Acts::eBoundLoc0) = 
	cluster->getActsLocalError(0,0) * Acts::UnitConstants::cm2;
      cov(Acts::eBoundLoc0, Acts::eBoundLoc1) =
	cluster->getActsLocalError(0,1) * Acts::UnitConstants::cm2;
      cov(Acts::eBoundLoc1, Acts::eBoundLoc0) = 
	cluster->getActsLocalError(1,0) * Acts::UnitConstants::cm2;
      cov(Acts::eBoundLoc1, Acts::eBoundLoc1) = 
	cluster->getActsLocalError(1,1) * Acts::UnitConstants::cm2;
    
      SourceLink sl(key, surf, loc, cov);
      
     
      if(Verbosity() > 4)
	{
	  std::cout << "Adding SL to track: "<< std::endl;
	  std::cout << "SL : " << sl.location().transpose() << std::endl
		    << sl.covariance().transpose() << std::endl 
		    << sl.geoId() << std::endl;
	}
      sourcelinks.push_back(sl);      
    }
 
  return sourcelinks;

}

void PHActsTrkFitter::getTrackFitResult(const FitResult &fitOutput,
				        SvtxTrack* track)
{
  /// Make a trajectory state for storage, which conforms to Acts track fit
  /// analysis tool
  std::vector<size_t> trackTips;
  trackTips.push_back(fitOutput.trackTip);
  ActsExamples::IndexedParams indexedParams;
  if (fitOutput.fittedParameters)
    {
       indexedParams.emplace(fitOutput.trackTip, 
			    fitOutput.fittedParameters.value());

      if (Verbosity() > 2)
        {
	  const auto& params = fitOutput.fittedParameters.value();
      
          std::cout << "Fitted parameters for track" << std::endl;
          std::cout << " position : " << params.position(m_tGeometry->geoContext).transpose()
	    
                    << std::endl;
	  std::cout << "charge: "<<params.charge()<<std::endl;
          std::cout << " momentum : " << params.momentum().transpose()
                    << std::endl;
	  std::cout << "For trackTip == " << fitOutput.trackTip << std::endl;
        }
    }
  auto trajectory = std::make_unique<Trajectory>(fitOutput.fittedStates,
						 trackTips, indexedParams, 
						 track->get_vertex_id());

  /// Get position, momentum from the Acts output. Update the values of
  /// the proto track
  auto updateTrackTimer = std::make_unique<PHTimer>("UpdateTrackTimer");
  updateTrackTimer->stop();
  updateTrackTimer->restart();
  if(fitOutput.fittedParameters)
    updateSvtxTrack(*trajectory, track);
  
  updateTrackTimer->stop();
  auto updateTime = updateTrackTimer->get_accumulated_time();
  
  if(Verbosity() > 0)
    std::cout << "PHActsTrkFitter update SvtxTrack time "
	      << updateTime << std::endl;

  if(m_timeAnalysis)
    h_updateTime->Fill(updateTime);
  
  return;
}

ActsExamples::TrkrClusterFittingAlgorithm::FitterResult PHActsTrkFitter::fitTrack(
          const SourceLinkVec& sourceLinks, 
	  const ActsExamples::TrackParameters& seed,
	  const Acts::KalmanFitterOptions<Acts::VoidOutlierFinder>& 
	         kfOptions, 
	  const SurfacePtrVec& surfSequence)
{
  if(m_fitSiliconMMs) 
    return m_fitCfg.dFit(sourceLinks, seed, kfOptions, surfSequence);  
  else
    return m_fitCfg.fit(sourceLinks, seed, kfOptions);
}

SourceLinkVec PHActsTrkFitter::getSurfaceVector(SourceLinkVec sourceLinks,
						SurfacePtrVec& surfaces)
{
   SourceLinkVec siliconMMSls;

  if(Verbosity() > 1)
    std::cout << "Sorting " << sourceLinks.size() << " SLs" << std::endl;
  
  for(auto sl : sourceLinks)
    {
      auto volume = sl.referenceSurface().geometryId().volume();

      if(Verbosity() > 1)
	std::cout<<"SL available on : " << sl.referenceSurface().geometryId()<<std::endl;
    
      /// If volume is not the TPC add the SL to the list
      if(volume != 14)
	{
	  siliconMMSls.push_back(sl);	
	  surfaces.push_back(&sl.referenceSurface());
	}
    }

  /// Surfaces need to be sorted in order, i.e. from smallest to
  /// largest radius extending from target surface
  /// Add a check to ensure this
  if(surfaces.size() > 0)
    checkSurfaceVec(surfaces);

  if(Verbosity() > 1)
    {
      for(auto surf : surfaces)
	{
	  std::cout << "Surface vector : " << surf->geometryId() << std::endl;
	}
    }

  return siliconMMSls;

}

void PHActsTrkFitter::checkSurfaceVec(SurfacePtrVec &surfaces)
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
	      if(Verbosity() > 2)
		std::cout << PHWHERE 
			  << "Surface not in order... removing surface" 
			  << surface->geometryId() << std::endl;
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
	      if(Verbosity() > 2)
		std::cout << PHWHERE 
			  << "Volume not in order... removing surface" 
			  << surface->geometryId() << std::endl;
	      surfaces.erase(surfaces.begin() + i);
	      /// Subtract one so we don't skip a surface
	      i--;
	      continue;
	    }
	}
    } 

}

void PHActsTrkFitter::updateSvtxTrack(Trajectory traj, 
				      SvtxTrack* track)
{
  const auto &[trackTips, mj] = traj.trajectory();
  /// only one track tip in the track fit Trajectory
  auto &trackTip = trackTips.front();

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
  SvtxTrackState_v1 out( pathlength);
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

  auto rotater = std::make_unique<ActsTransformations>();
  rotater->setVerbosity(Verbosity());
  
  float dca3Dxy = NAN;
  float dca3Dz = NAN;
  float dca3DxyCov = NAN;
  float dca3DzCov = NAN;
      
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
	      track->set_acts_covariance(i,j, params.covariance().value()(i,j));
	    }
	}
    
      unsigned int vertexId = track->get_vertex_id();
      const SvtxVertex *svtxVertex = m_vertexMap->get(vertexId);
      Acts::Vector3D vertex(
		  svtxVertex->get_x() * Acts::UnitConstants::cm, 
		  svtxVertex->get_y() * Acts::UnitConstants::cm, 
		  svtxVertex->get_z() * Acts::UnitConstants::cm);
      rotater->calculateDCA(params, vertex, rotatedCov,
			    m_tGeometry->geoContext, 
			    dca3Dxy, dca3Dz, 
			    dca3DxyCov, dca3DzCov);
    }
 
  /// Set the DCA here. The DCA will be updated after the final
  /// vertex fitting in PHActsVertexFinder
  track->set_dca3d_xy(dca3Dxy / Acts::UnitConstants::cm);
  track->set_dca3d_z(dca3Dz / Acts::UnitConstants::cm);

  /// The covariance that goes into the rotater is already in sphenix
  /// units, so we don't need to convert back
  track->set_dca3d_xy_error(sqrt(dca3DxyCov));
  track->set_dca3d_z_error(sqrt(dca3DzCov));
  
  // Also need to update the state list and cluster ID list for all measurements associated with the acts track  
  // loop over acts track states, copy over to SvtxTrackStates, and add to SvtxTrack

  auto trackStateTimer = std::make_unique<PHTimer>("TrackStateTimer");
  trackStateTimer->stop();
  trackStateTimer->restart();
  
  if(m_fillSvtxTrackStates)
    rotater->fillSvtxTrackStates(traj, trackTip, track,
				 m_tGeometry->geoContext);  

  trackStateTimer->stop();
  auto stateTime = trackStateTimer->get_accumulated_time();
  
  if(Verbosity() > 0)
    std::cout << "PHActsTrkFitter update SvtxTrackStates time "
	      << stateTime << std::endl;

  if(m_timeAnalysis)
    h_stateTime->Fill(stateTime);

  if(Verbosity() > 2)
    {  
      std::cout << " Identify fitted track after updating track states:" 
		<< std::endl;
      track->identify();
      std::cout << " cluster keys size " << track->size_cluster_keys() 
		<< std::endl;  
    }
 
 return;
  
}

Acts::BoundSymMatrix PHActsTrkFitter::setDefaultCovariance()
{
  Acts::BoundSymMatrix cov;
   
  /// Acts cares about the track covariance as it helps the KF
  /// know whether or not to trust the initial track seed or not.
  /// We reset it here to some loose values as it helps Acts improve
  /// the fitting. 
  /// If the covariance is too loose, it won't be able to propagate,
  /// but if it is too tight, it will just "believe" the track seed over
  /// the hit data
 
  /// If we are using distortions, then we need to blow up the covariance
  /// a bit since the seed was created with distorted TPC clusters
  if(m_fitSiliconMMs)
    cov << 1000 * Acts::UnitConstants::um, 0., 0., 0., 0., 0.,
           0., 1000 * Acts::UnitConstants::um, 0., 0., 0., 0.,
           0., 0., 0.1, 0., 0., 0.,
           0., 0., 0., 0.1, 0., 0.,
           0., 0., 0., 0., 0.005 , 0.,
           0., 0., 0., 0., 0., 1.;
  else
    cov << 1000 * Acts::UnitConstants::um, 0., 0., 0., 0., 0.,
           0., 1000 * Acts::UnitConstants::um, 0., 0., 0., 0.,
           0., 0., 0.05, 0., 0., 0.,
           0., 0., 0., 0.05, 0., 0.,
           0., 0., 0., 0., 0.00005 , 0.,
           0., 0., 0., 0., 0., 1.;

  return cov;
}
    

int PHActsTrkFitter::createNodes(PHCompositeNode* topNode)
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

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsTrkFitter::getNodes(PHCompositeNode* topNode)
{
  m_surfMaps = findNode::getClass<ActsSurfaceMaps>(topNode, "ActsSurfaceMaps");
  if(!m_surfMaps)
    {
      std::cout << PHWHERE << "ActsSurfaceMaps not on node tree, bailing."
		<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  m_clusterContainer = findNode::getClass<TrkrClusterContainer>(topNode,"TRKR_CLUSTER");
  if(!m_clusterContainer)
    {
      std::cout << PHWHERE 
		<< "No trkr cluster container, exiting." << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  m_tGeometry = findNode::getClass<ActsTrackingGeometry>(topNode, "ActsTrackingGeometry");
  if(!m_tGeometry)
    {
      std::cout << "ActsTrackingGeometry not on node tree. Exiting."
		<< std::endl;
      
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  
  m_vertexMap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
  if(!m_vertexMap)
    {
      std::cout << PHWHERE << "SvtxVertexMap not on node tree, bailing"
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

