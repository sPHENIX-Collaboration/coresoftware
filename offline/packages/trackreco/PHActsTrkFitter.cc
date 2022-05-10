/*!
 *  \file		PHActsTrkFitter.C
 *  \brief		Refit SvtxTracks with PHActs.
 *  \details	Refit SvtxTracks with PHActs.
 *  \author	        Tony Frawley <afrawley@fsu.edu>
 */

#include "PHActsTrkFitter.h"

/// Tracking includes
#include <tpc/TpcDefs.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/MvtxDefs.h>
#include <trackbase/InttDefs.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrack_v3.h>
#include <trackbase_historic/SvtxTrackState_v1.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackMap_v1.h>
#include <trackbase_historic/ActsTransformations.h>

#include <micromegas/MicromegasDefs.h>

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
#include <Acts/TrackFitting/GainMatrixSmoother.hpp>
#include <Acts/TrackFitting/GainMatrixUpdater.hpp>

#include <ActsExamples/EventData/Index.hpp>

#include <TDatabasePDG.h>

#include <cmath>
#include <iostream>
#include <vector>



PHActsTrkFitter::PHActsTrkFitter(const std::string& name)
  : SubsysReco(name)
  , m_trajectories(nullptr)
{}

int PHActsTrkFitter::InitRun(PHCompositeNode* topNode)
{
  if(Verbosity() > 1)
    std::cout << "Setup PHActsTrkFitter" << std::endl;

   if(createNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
    return Fun4AllReturnCodes::ABORTEVENT;

  if (getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
    return Fun4AllReturnCodes::ABORTEVENT;
  
  m_fitCfg.fit = ActsExamples::TrackFittingAlgorithm::makeTrackFitterFunction(
               m_tGeometry->tGeometry,
	       m_tGeometry->magField);

  m_fitCfg.dFit = ActsExamples::TrackFittingAlgorithm::makeTrackFitterFunction(
	       m_tGeometry->magField);

  m_outlierFinder.verbosity = Verbosity();
  std::map<long unsigned int, float> chi2Cuts;
  chi2Cuts.insert(std::make_pair(10,4));
  chi2Cuts.insert(std::make_pair(12,4));
  chi2Cuts.insert(std::make_pair(14,9));
  chi2Cuts.insert(std::make_pair(16,4));
  m_outlierFinder.chi2Cuts = chi2Cuts;

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

int PHActsTrkFitter::process_event(PHCompositeNode */*topNode*/)
{
  PHTimer eventTimer("eventTimer");
  eventTimer.stop();
  eventTimer.restart();
  
  m_event++;

  auto logLevel = Acts::Logging::FATAL;

  if (Verbosity() > 1)
  {
    std::cout << PHWHERE << "Events processed: " << m_event << std::endl;
    std::cout << "Start PHActsTrkFitter::process_event" << std::endl;
    if(Verbosity() > 4)
      logLevel = Acts::Logging::VERBOSE;
  }

  /// Fill an additional track map if using the acts evaluator
  /// for proto track comparison to fitted track
  if(m_actsEvaluator)
    {
      /// wipe at the beginning of every new fit pass, so that the seeds 
      /// are whatever is currently in SvtxTrackMap
      m_seedTracks->clear();
      for(const auto& [key, track] : *m_trackMap)
	{
	  m_seedTracks->insert(track);
	}
    }

  loopTracks(logLevel);

  eventTimer.stop();
  auto eventTime = eventTimer.get_accumulated_time();

  if(Verbosity() > 1)
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

int PHActsTrkFitter::ResetEvent(PHCompositeNode */*topNode*/)
{
  
  if(Verbosity() > 1)
    {
      std::cout << "Reset PHActsTrkFitter" << std::endl;

    }
  
  m_trajectories->clear();
    
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsTrkFitter::End(PHCompositeNode */*topNode*/)
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

  for(const auto& [trackKey, track] : *m_trackMap)
    {
      if(!track)
	{ continue; }

      PHTimer trackTimer("TrackTimer");
      trackTimer.stop();
      trackTimer.restart();
      ActsExamples::MeasurementContainer measurements;
      auto sourceLinks = getSourceLinks(track, measurements);
  
      if(sourceLinks.size() == 0) { continue; }

      /// If using directed navigation, collect surface list to navigate
      SurfacePtrVec surfaces;
      if(m_fitSiliconMMs)
      {
        sourceLinks = getSurfaceVector(sourceLinks, surfaces);
        
        // skip if there is no surfaces
        if( surfaces.empty() ) continue;
        
        // make sure micromegas are in the tracks, if required
        if( m_useMicromegas &&
          std::none_of( surfaces.begin(), surfaces.end(), [this]( const auto& surface )
          { return m_surfMaps->isMicromegasSurface( surface ); } ) )
        { continue; }
      }

      Acts::Vector3 momentum(track->get_px(), 
			     track->get_py(), 
			     track->get_pz());

      // z value of track is from silicon stub, so it is independent of bunch crossing
      Acts::Vector3 position(track->get_x() * Acts::UnitConstants::cm,
			     track->get_y() * Acts::UnitConstants::cm,
			     track->get_z() * Acts::UnitConstants::cm);

      auto pSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(
					  position);
      auto actsFourPos = Acts::Vector4(position(0), position(1),
				       position(2),
				       10 * Acts::UnitConstants::ns);
      Acts::BoundSymMatrix cov = setDefaultCovariance();
 
      int charge = track->get_charge();
      if(m_fieldMap.find("3d") != std::string::npos)
	{ charge *= -1; }

      /// Acts requires a wrapped vector, so we need to replace the
      /// std::vector contents with a wrapper vector to get the memory
      /// access correct
      std::vector<std::reference_wrapper<const SourceLink>>  wrappedSls;
      for(const auto& sl : sourceLinks)
	{ wrappedSls.push_back(std::cref(sl)); }
           
      if(Verbosity() > 2)
	{ printTrackSeed(track); }

      /// Reset the track seed with the dummy covariance
      auto seed = ActsExamples::TrackParameters::create(pSurface,
							m_tGeometry->geoContext,
							actsFourPos,
							momentum,
							charge / track->get_p(),
							cov).value();
      
      /// Set host of propagator options for Acts to do e.g. material integration
      Acts::PropagatorPlainOptions ppPlainOptions;
      ppPlainOptions.absPdgCode = m_pHypothesis;
      ppPlainOptions.mass = 
	TDatabasePDG::Instance()->GetParticle(m_pHypothesis)->Mass() * Acts::UnitConstants::GeV;
       
      Acts::KalmanFitterExtensions extensions;
      ActsExamples::MeasurementCalibrator calibrator{measurements};
      extensions.calibrator.connect<&ActsExamples::MeasurementCalibrator::calibrate>(&calibrator);
     
      if(m_useOutlierFinder)
	{ 
	  extensions.outlierFinder.connect<&ResidualOutlierFinder::operator()>(&m_outlierFinder);
	}

      Acts::GainMatrixUpdater kfUpdater;
      Acts::GainMatrixSmoother kfSmoother;
      extensions.updater.connect<&Acts::GainMatrixUpdater::operator()>(&kfUpdater);
      extensions.smoother.connect<&Acts::GainMatrixSmoother::operator()>(&kfSmoother);

      Acts::KalmanFitterOptions kfOptions(m_tGeometry->geoContext,
					  m_tGeometry->magFieldContext,
					  m_tGeometry->calibContext,
					  extensions,
					  Acts::LoggerWrapper(*logger),
					  ppPlainOptions,
					  &(*pSurface));
 
      kfOptions.multipleScattering = true;
      kfOptions.energyLoss = true;

      PHTimer fitTimer("FitTimer");
      fitTimer.stop();
      fitTimer.restart();
      auto result = fitTrack(wrappedSls, seed, kfOptions,
			     surfaces);
      fitTimer.stop();
      auto fitTime = fitTimer.get_accumulated_time();
      
      if(Verbosity() > 1)
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
	  
	  if(m_fitSiliconMMs)
	    {
        std::unique_ptr<SvtxTrack_v3> newTrack( static_cast<SvtxTrack_v3*>(track->CloneMe()) );
	      if( getTrackFitResult(fitOutput, newTrack.get()) )
        { m_directedTrackMap->insert(newTrack.release()); } 
	    }
	  else
	    {
	      getTrackFitResult(fitOutput, track);
	    }
	}
      else if (!m_fitSiliconMMs)
	{
	  /// Track fit failed, get rid of the track from the map
	  badTracks.push_back(trackKey);
	  m_nBadFits++;
	  if(Verbosity() > 1)
	    { 
	      std::cout << "Track fit failed for track " << trackKey 
			<< " with Acts error message " 
			<< result.error() << ", " << result.error().message()
			<< std::endl;
	    }
	}

      trackTimer.stop();
      auto trackTime = trackTimer.get_accumulated_time();
      
      if(Verbosity() > 1)
	std::cout << "PHActsTrkFitter total single track time "
		  << trackTime << std::endl;
    
    }

  /// Now erase bad tracks from the track map
  for(const auto& key : badTracks)
    {
      if(Verbosity() > 2)
	{ std::cout << "Erasing bad track " << key << std::endl; }
      m_trackMap->erase(key);
    }

  return;

}


//___________________________________________________________________________________
Surface PHActsTrkFitter::getSurface(TrkrDefs::cluskey cluskey, TrkrDefs::subsurfkey surfkey) const
{
  const auto trkrid = TrkrDefs::getTrkrId(cluskey);
  const auto hitsetkey = TrkrDefs::getHitSetKeyFromClusKey(cluskey);

  switch( trkrid )
  {
    case TrkrDefs::TrkrId::micromegasId: return getMMSurface( hitsetkey );
    case TrkrDefs::TrkrId::tpcId: return getTpcSurface(hitsetkey, surfkey);
    case TrkrDefs::TrkrId::mvtxId:
    case TrkrDefs::TrkrId::inttId:
    {
      return getSiliconSurface(hitsetkey);
    }
  }
  
  // unreachable
  return nullptr;
  
}

//___________________________________________________________________________________
Surface PHActsTrkFitter::getSiliconSurface(TrkrDefs::hitsetkey hitsetkey) const
{
  unsigned int trkrid =  TrkrDefs::getTrkrId(hitsetkey);
  TrkrDefs::hitsetkey tmpkey = hitsetkey;

  if(trkrid == TrkrDefs::inttId)
    {
      // Set the hitsetkey crossing to zero
      tmpkey = InttDefs::resetCrossingHitSetKey(hitsetkey);
    }

  if(trkrid == TrkrDefs::mvtxId)
    {
      // Set the hitsetkey crossing to zero
      tmpkey = MvtxDefs::resetStrobeHitSetKey(hitsetkey);
    }

  auto surfMap = m_surfMaps->siliconSurfaceMap;
  auto iter = surfMap.find(tmpkey);
  if(iter != surfMap.end())
    {
      return iter->second;
    }
  
  /// If it can't be found, return nullptr
  {
    std::cout << PHWHERE << "Failed to find silicon surface for hitsetkey " << hitsetkey << " tmpkey " << tmpkey << std::endl;
  }
  
  return nullptr;

}

//___________________________________________________________________________________
Surface PHActsTrkFitter::getTpcSurface(TrkrDefs::hitsetkey hitsetkey, TrkrDefs::subsurfkey surfkey) const
{
  unsigned int layer = TrkrDefs::getLayer(hitsetkey);
  const auto iter = m_surfMaps->tpcSurfaceMap.find(layer);
  if(iter != m_surfMaps->tpcSurfaceMap.end())
  {
    auto surfvec = iter->second;
    return surfvec.at(surfkey);
  }
  
  /// If it can't be found, return nullptr to skip this cluster
  return nullptr;
}

//___________________________________________________________________________________
Surface PHActsTrkFitter::getMMSurface(TrkrDefs::hitsetkey hitsetkey) const
{
  const auto iter = m_surfMaps->mmSurfaceMap.find( hitsetkey );
  return (iter == m_surfMaps->mmSurfaceMap.end()) ? nullptr:iter->second;
}



//___________________________________________________________________________________
SourceLinkVec PHActsTrkFitter::getSourceLinks(SvtxTrack* track,
				   ActsExamples::MeasurementContainer& measurements)
{

  SourceLinkVec sourcelinks;

  short int crossing = track->get_crossing();
  if(Verbosity() > 0) 
    {
      std::cout << "tpcid " << track->get_id() << " crossing " << crossing << std::endl; 
      //track->identify();
    }

  if(crossing == SHRT_MAX) return sourcelinks;

  // loop over all clusters
  std::vector<std::pair<TrkrDefs::cluskey, Acts::Vector3>> global_raw;

  //int iter = 0;
  for (SvtxTrack::ConstClusterKeyIter clusIter = track->begin_cluster_keys();
       clusIter != track->end_cluster_keys();
       ++clusIter)
    {
      auto key = *clusIter;
      auto cluster = m_clusterContainer->findCluster(key);
      if(!cluster)
	{
	  if(Verbosity() > 0) std::cout << "Failed to get cluster with key " << key << " for track " << track->get_id() << std::endl;
	  continue;
	}

      auto subsurfkey = cluster->getSubSurfKey();
      
      /// Make a safety check for clusters that couldn't be attached
      /// to a surface
      auto surf = getSurface(key, subsurfkey);
      if(!surf)
	continue;

      unsigned int trkrid = TrkrDefs::getTrkrId(key);
      unsigned int side = TpcDefs::getSide(key);

      // For the TPC, cluster z has to be corrected for the crossing z offset, distortion, and TOF z offset 
      // we do this locally here and do not modify the cluster, since the cluster may be associated with multiple silicon tracks  
      
      // transform to global coordinates for z correction 
      ActsTransformations transformer;
      auto global = transformer.getGlobalPosition(key, cluster,
						  m_surfMaps,
						  m_tGeometry);
      
      if(Verbosity() > 0)
	{
	  std::cout << " zinit " << global[2] << " xinit " << global[0] << " yinit " << global[1] << " side " << side << " crossing " << crossing 
		    << " cluskey " << key << " subsurfkey " << subsurfkey << std::endl;
	}
      
      if(trkrid ==  TrkrDefs::tpcId)
	{	  
	  // make all corrections to global position of TPC cluster

	  float z = m_clusterCrossingCorrection.correctZ(global[2], side, crossing);
	  global[2] = z;
	  
	  // apply distortion corrections
	  if(_dcc_static) { global = _distortionCorrection.get_corrected_position( global, _dcc_static ); }
	  if(_dcc_average) { global = _distortionCorrection.get_corrected_position( global, _dcc_average ); }
	  if(_dcc_fluctuation) { global = _distortionCorrection.get_corrected_position( global, _dcc_fluctuation ); }
	 
	  // add the global positions to a vector to give to the cluster mover
	  global_raw.push_back(std::make_pair(key, global));
	}
      else
	{
	  // silicon cluster, no corrections needed
	  global_raw.push_back(std::make_pair(key, global));
	}
      
    }	  // end loop over clusters here
  
  // move the cluster positions back to the original readout surface
  std::vector<std::pair<TrkrDefs::cluskey,Acts::Vector3>> global_moved = _clusterMover.processTrack(global_raw);
  
  // loop over global positions returned by cluster mover
  for(int i=0; i<global_moved.size(); ++i)
    {
      TrkrDefs::cluskey cluskey = global_moved[i].first;
      Acts::Vector3 global = global_moved[i].second;

      Surface surf;
      TrkrDefs::subsurfkey subsurfkey;

      unsigned int trkrid = TrkrDefs::getTrkrId(cluskey);
      if(trkrid ==  TrkrDefs::tpcId)
	{      
	  // get the new surface corresponding to this global position
	  TrkrDefs::hitsetkey tpcHitSetKey = TrkrDefs::getHitSetKeyFromClusKey(cluskey);
	  surf = get_tpc_surface_from_coords(tpcHitSetKey,
					     global,
					     m_surfMaps,
					     m_tGeometry,
					     subsurfkey);
	  
	  if(!surf)
	    {
	      /// If the surface can't be found, we can't track with it. So 
	      /// just continue and don't add the cluster to the source links
	      if(Verbosity() > 0) std::cout << PHWHERE << "Failed to find surface for cluster " << cluskey << std::endl;
	      continue;
	    }
	  
	  if(Verbosity() > 0)
	    {
	      unsigned int side = TpcDefs::getSide(cluskey);       
	      std::cout << "      global z corrected and moved " << global[2] << " xcorr " << global[0] << " ycorr " << global[1] << " side " << side << " crossing " << crossing 
			<< " cluskey " << cluskey << " subsurfkey " << subsurfkey << std::endl;
	    }
	}
      else
	{
	  // silicon cluster, no changes possible
	  auto cluster = m_clusterContainer->findCluster(cluskey);
	  subsurfkey = cluster->getSubSurfKey();
	  surf = getSurface(cluskey, subsurfkey);	  
	}

      // get local coordinates
      Acts::Vector2 localPos;
      Acts::Vector3 normal = surf->normal(m_tGeometry->geoContext);
      auto local = surf->globalToLocal(m_tGeometry->geoContext,
				       global * Acts::UnitConstants::cm,
				       normal);
      
      if(local.ok())
	{
	  localPos = local.value() / Acts::UnitConstants::cm;
	}
      else
	{
	  /// otherwise take the manual calculation
	  Acts::Vector3 center = surf->center(m_tGeometry->geoContext)/Acts::UnitConstants::cm;
	  double clusRadius = sqrt(global[0]*global[0] + global[1]*global[1]);
	  double clusphi = atan2(global[1], global[0]);
	  double rClusPhi = clusRadius * clusphi;
	  double surfRadius = sqrt(center(0)*center(0) + center(1)*center(1));
	  double surfPhiCenter = atan2(center[1], center[0]);
	  double surfRphiCenter = surfPhiCenter * surfRadius;
	  double surfZCenter = center[2];
	  
	  localPos(0) = rClusPhi - surfRphiCenter;
	  localPos(1) = global[2] - surfZCenter; 
	}
      
      auto cluster = m_clusterContainer->findCluster(cluskey);
      if(Verbosity() > 0)
	{
	  std::cout << " cluster local X " << cluster->getLocalX() << " cluster local Y " << cluster->getLocalY() << std::endl;
	  std::cout << " new      local X " << localPos(0) << " new       local Y " << localPos(1) << std::endl;
	}
      
      
      Acts::ActsVector<2> loc;
      loc[Acts::eBoundLoc0] = localPos(0) * Acts::UnitConstants::cm;
      loc[Acts::eBoundLoc1] = localPos(1) * Acts::UnitConstants::cm;
      std::array<Acts::BoundIndices,2> indices;
      indices[0] = Acts::BoundIndices::eBoundLoc0;
      indices[1] = Acts::BoundIndices::eBoundLoc1;
      Acts::ActsSymMatrix<2> cov = Acts::ActsSymMatrix<2>::Zero();
      cov(Acts::eBoundLoc0, Acts::eBoundLoc0) = 
	cluster->getActsLocalError(0,0) * Acts::UnitConstants::cm2;
      cov(Acts::eBoundLoc0, Acts::eBoundLoc1) =
	cluster->getActsLocalError(0,1) * Acts::UnitConstants::cm2;
      cov(Acts::eBoundLoc1, Acts::eBoundLoc0) = 
	cluster->getActsLocalError(1,0) * Acts::UnitConstants::cm2;
      cov(Acts::eBoundLoc1, Acts::eBoundLoc1) = 
	cluster->getActsLocalError(1,1) * Acts::UnitConstants::cm2;
      ActsExamples::Index index = measurements.size();
      
      SourceLink sl(surf->geometryId(), index, cluskey);
      
      Acts::Measurement<Acts::BoundIndices,2> meas(sl, indices, loc, cov);
      if(Verbosity() > 3)
	{
	  std::cout << "source link " << sl.index() << ", loc : " 
		    << loc.transpose() << std::endl 
		    << ", cov : " << cov.transpose() << std::endl
		    << " geo id " << sl.geometryId() << std::endl;
	  std::cout << "Surface : " << std::endl;
	  surf.get()->toStream(m_tGeometry->geoContext, std::cout);
	  std::cout << std::endl;
	  std::cout << "Cluster error " << cluster->getRPhiError() << " , " << cluster->getZError() << std::endl;
	  std::cout << "For key " << cluskey << " with local pos " << std::endl
		    << localPos(0) << ", " << localPos(1)
		    << std::endl;
	}
      
      sourcelinks.push_back(sl);
      measurements.push_back(meas);
      //iter++;
    }
  
    
  return sourcelinks;
}

bool PHActsTrkFitter::getTrackFitResult(const FitResult &fitOutput,
				        SvtxTrack* track)
{
  /// Make a trajectory state for storage, which conforms to Acts track fit
  /// analysis tool
  std::vector<size_t> trackTips;
  trackTips.reserve(1);
  trackTips.emplace_back(fitOutput.lastMeasurementIndex);
  ActsExamples::Trajectories::IndexedParameters indexedParams;
  if (fitOutput.fittedParameters)
    {
       indexedParams.emplace(fitOutput.lastMeasurementIndex, 
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
	  std::cout << "For trackTip == " << fitOutput.lastMeasurementIndex << std::endl;
        }
    }
  else 
    {
      /// Track fit failed in some way if there are no fit parameters. Remove
      m_trackMap->erase(track->get_id());
      std::cout << " track fit failed for track " << track->get_id() << std::endl;
      return false;
    }

  Trajectory trajectory(fitOutput.fittedStates,
			trackTips, indexedParams);
 
  m_trajectories->insert(std::make_pair(track->get_id(), trajectory));
    

  /// Get position, momentum from the Acts output. Update the values of
  /// the proto track
  PHTimer updateTrackTimer("UpdateTrackTimer");
  updateTrackTimer.stop();
  updateTrackTimer.restart();
  if(fitOutput.fittedParameters)
    { updateSvtxTrack(trajectory, track); }
  
  updateTrackTimer.stop();
  auto updateTime = updateTrackTimer.get_accumulated_time();
  
  if(Verbosity() > 1)
    std::cout << "PHActsTrkFitter update SvtxTrack time "
	      << updateTime << std::endl;

  if(m_timeAnalysis)
    h_updateTime->Fill(updateTime);
  
  return true;
}

ActsExamples::TrackFittingAlgorithm::TrackFitterResult PHActsTrkFitter::fitTrack(
    const std::vector<std::reference_wrapper<const SourceLink>>& sourceLinks, 
    const ActsExamples::TrackParameters& seed,
    const ActsExamples::TrackFittingAlgorithm::TrackFitterOptions& kfOptions, 
    const SurfacePtrVec& surfSequence)
{

  if(m_fitSiliconMMs) 
    { return (*m_fitCfg.dFit)(sourceLinks, seed, kfOptions, surfSequence); }

  return (*m_fitCfg.fit)(sourceLinks, seed, kfOptions); 
}

SourceLinkVec PHActsTrkFitter::getSurfaceVector(const SourceLinkVec& sourceLinks,
						SurfacePtrVec& surfaces) const
{
  SourceLinkVec siliconMMSls;

  if(Verbosity() > 1)
    std::cout << "Sorting " << sourceLinks.size() << " SLs" << std::endl;
  
  for(const auto& sl : sourceLinks)
    {
      if(Verbosity() > 1)
	{ std::cout << "SL available on : " << sl.geometryId() << std::endl; } 
      const auto surf = m_tGeometry->tGeometry->findSurface(sl.geometryId());
      // skip TPC surfaces
      if( m_surfMaps->isTpcSurface( surf ) ) continue;
      
      // also skip micromegas surfaces if not used
      if( m_surfMaps->isMicromegasSurface( surf ) && !m_useMicromegas ) continue;

      // update vectors
      siliconMMSls.push_back(sl);
      surfaces.push_back(surf);
    }

  /// Surfaces need to be sorted in order, i.e. from smallest to
  /// largest radius extending from target surface
  /// Add a check to ensure this
  if(!surfaces.empty())
  { checkSurfaceVec(surfaces); }

  if(Verbosity() > 1)
    {
      for(const auto& surf : surfaces)
	{
	  std::cout << "Surface vector : " << surf->geometryId() << std::endl;
	}
    }

  return siliconMMSls;
}

void PHActsTrkFitter::checkSurfaceVec(SurfacePtrVec &surfaces) const
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
  const auto& mj = traj.multiTrajectory();
  const auto& tips = traj.tips();

  /// only one track tip in the track fit Trajectory
  auto &trackTip = tips.front();

  if(Verbosity() > 2)
    {
      std::cout << "Identify (proto) track before updating with acts results " << std::endl;
      track->identify();
      std::cout << " cluster keys size " << track->size_cluster_keys() << std::endl;  
    }

  
  if(!m_fitSiliconMMs)
    { track->clear_states(); }

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

  ActsTransformations rotater;
  rotater.setVerbosity(Verbosity());
 
  if(params.covariance())
    {     
      auto rotatedCov = rotater.rotateActsCovToSvtxTrack(params);
      
      for(int i = 0; i < 6; i++)
	{
	  for(int j = 0; j < 6; j++)
	    {
	      track->set_error(i,j, rotatedCov(i,j));
	      track->set_acts_covariance(i,j, 
					 params.covariance().value()(i,j));
	    }
	} 
    }

  // Also need to update the state list and cluster ID list for all measurements associated with the acts track  
  // loop over acts track states, copy over to SvtxTrackStates, and add to SvtxTrack

  PHTimer trackStateTimer("TrackStateTimer");
  trackStateTimer.stop();
  trackStateTimer.restart();

  if(m_fillSvtxTrackStates)
    { 
      rotater.fillSvtxTrackStates(mj, trackTip, track,
				  m_tGeometry->geoContext);  
    }
  
  trackStateTimer.stop();
  auto stateTime = trackStateTimer.get_accumulated_time();
  
  if(Verbosity() > 1)
    std::cout << "PHActsTrkFitter update SvtxTrackStates time "
	      << stateTime << std::endl;

  if(m_timeAnalysis)
    { h_stateTime->Fill(stateTime); }

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

Acts::BoundSymMatrix PHActsTrkFitter::setDefaultCovariance() const
{
  Acts::BoundSymMatrix cov = Acts::BoundSymMatrix::Zero();
   
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
    {
      cov << 1000 * Acts::UnitConstants::um, 0., 0., 0., 0., 0.,
           0., 1000 * Acts::UnitConstants::um, 0., 0., 0., 0.,
           0., 0., 0.1, 0., 0., 0.,
           0., 0., 0., 0.1, 0., 0.,
           0., 0., 0., 0., 0.005 , 0.,
           0., 0., 0., 0., 0., 1.;
    }
  else
    {
      double sigmaD0 = 50 * Acts::UnitConstants::um;
      double sigmaZ0 = 50 * Acts::UnitConstants::um;
      double sigmaPhi = 1 * Acts::UnitConstants::degree;
      double sigmaTheta = 1 * Acts::UnitConstants::degree;
      double sigmaT = 1. * Acts::UnitConstants::ns;
     
      cov(Acts::eBoundLoc0, Acts::eBoundLoc0) = sigmaD0 * sigmaD0;
      cov(Acts::eBoundLoc1, Acts::eBoundLoc1) = sigmaZ0 * sigmaZ0;
      cov(Acts::eBoundTime, Acts::eBoundTime) = sigmaT * sigmaT;
      cov(Acts::eBoundPhi, Acts::eBoundPhi) = sigmaPhi * sigmaPhi;
      cov(Acts::eBoundTheta, Acts::eBoundTheta) = sigmaTheta * sigmaTheta;
      /// Acts takes this value very seriously - tuned to be in a "sweet spot"
      cov(Acts::eBoundQOverP, Acts::eBoundQOverP) = 0.0001;

    }

  return cov;
}

void PHActsTrkFitter::printTrackSeed(const SvtxTrack* seed) const
{
  std::cout << PHWHERE << " Processing proto track with id: " 
	    << seed->get_id() << " and  position:" 
	    << seed->get_x() << ", " << seed->get_y() << ", " 
	    << seed->get_z() << std::endl 
	    << "momentum: " << seed->get_px() << ", " << seed->get_py()
	    << ", " << seed->get_pz() << std::endl
	    << "charge : " << seed->get_charge()
	    << std::endl;  
}
    
int PHActsTrkFitter::createNodes(PHCompositeNode* topNode)
{

  PHNodeIterator iter(topNode);
  
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));

  if (!dstNode)
  {
    std::cerr << "DST node is missing, quitting" << std::endl;
    throw std::runtime_error("Failed to find DST node in PHActsTrkFitter::createNodes");
  }
  
  PHCompositeNode *svtxNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "SVTX"));

  if (!svtxNode)
  {
    svtxNode = new PHCompositeNode("SVTX");
    dstNode->addNode(svtxNode);
  }

  if(m_fitSiliconMMs)
    {
      m_directedTrackMap = findNode::getClass<SvtxTrackMap>(topNode,
							    "SvtxSiliconMMTrackMap");
      if(!m_directedTrackMap)
	{
	  /// Copy this trackmap, then use it for the rest of processing
	  m_directedTrackMap = new SvtxTrackMap_v1;

	  PHIODataNode<PHObject> *trackNode = 
	    new PHIODataNode<PHObject>(m_directedTrackMap,"SvtxSiliconMMTrackMap","PHObject");
	  svtxNode->addNode(trackNode);
	} 
    }

  m_trajectories = findNode::getClass<std::map<const unsigned int, Trajectory>>(topNode, "ActsTrajectories");
  if(!m_trajectories)
    {
      m_trajectories = new std::map<const unsigned int, Trajectory>;
      PHDataNode<std::map<const unsigned int, Trajectory>> *node = 
	new PHDataNode<std::map<const unsigned int, Trajectory>>(m_trajectories, "ActsTrajectories");
      svtxNode->addNode(node);
      
    }
  
  if(m_actsEvaluator)
    {
      m_seedTracks = findNode::getClass<SvtxTrackMap>(topNode,_seed_track_map_name);
      
      if(!m_seedTracks)
	{
	  m_seedTracks = new SvtxTrackMap_v1;
	  
	  PHIODataNode<PHObject> *seedNode = 
	    new PHIODataNode<PHObject>(m_seedTracks,_seed_track_map_name,"PHObject");
	  svtxNode->addNode(seedNode);
	}
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
 
  m_trackMap = findNode::getClass<SvtxTrackMap>(topNode, _track_map_name);
  
  if(!m_trackMap)
    {
      std::cout << PHWHERE << "SvtxTrackMap not found on node tree. Exiting."
		<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

 // tpc distortion corrections
  _dcc_static = findNode::getClass<TpcDistortionCorrectionContainer>(topNode,"TpcDistortionCorrectionContainerStatic");
  if( _dcc_static )
    { 
      std::cout << PHWHERE << "  found static TPC distortion correction container" << std::endl; 
    }
  _dcc_average = findNode::getClass<TpcDistortionCorrectionContainer>(topNode,"TpcDistortionCorrectionContainerAverage");
  if( _dcc_average )
    { 
      std::cout << PHWHERE << "  found average TPC distortion correction container" << std::endl; 
    }
  _dcc_fluctuation = findNode::getClass<TpcDistortionCorrectionContainer>(topNode,"TpcDistortionCorrectionContainerFluctuation");
  if( _dcc_fluctuation )
    { 
      std::cout << PHWHERE << "  found fluctuation TPC distortion correction container" << std::endl; 
    }

  return Fun4AllReturnCodes::EVENT_OK;
}

Surface PHActsTrkFitter::get_tpc_surface_from_coords(TrkrDefs::hitsetkey hitsetkey,
						       Acts::Vector3 world,
						       ActsSurfaceMaps *surfMaps,
						       ActsTrackingGeometry *tGeometry,
						       TrkrDefs::subsurfkey& subsurfkey)
{
  unsigned int layer = TrkrDefs::getLayer(hitsetkey);
  std::map<unsigned int, std::vector<Surface>>::iterator mapIter;
  mapIter = surfMaps->tpcSurfaceMap.find(layer);
  
  if(mapIter == surfMaps->tpcSurfaceMap.end())
    {
      std::cout << PHWHERE 
		<< "Error: hitsetkey not found in clusterSurfaceMap, hitsetkey = "
		<< hitsetkey << std::endl;
      return nullptr;
    }
  
  double world_phi = atan2(world[1], world[0]);
  double world_z = world[2];
  
  std::vector<Surface> surf_vec = mapIter->second;
  unsigned int surf_index = 999;
    
  // Predict which surface index this phi and z will correspond to
  // assumes that the vector elements are ordered positive z, -pi to pi, then negative z, -pi to pi
  double fraction =  (world_phi + M_PI) / (2.0 * M_PI);
  double rounded_nsurf = round( (double) (surf_vec.size()/2) * fraction  - 0.5);
  unsigned int nsurf = (unsigned int) rounded_nsurf; 
  if(world_z < 0)
    nsurf += surf_vec.size()/2;

  Surface this_surf = surf_vec[nsurf];
      
  auto vec3d = this_surf->center(tGeometry->geoContext);
  std::vector<double> surf_center = {vec3d(0) / 10.0, vec3d(1) / 10.0, vec3d(2) / 10.0};  // convert from mm to cm
  double surf_z = surf_center[2];
  double surf_phi = atan2(surf_center[1], surf_center[0]);
  double surfStepPhi = tGeometry->tpcSurfStepPhi;
  double surfStepZ = tGeometry->tpcSurfStepZ;

  if( (world_phi > surf_phi - surfStepPhi / 2.0 && world_phi < surf_phi + surfStepPhi / 2.0 ) &&
      (world_z > surf_z - surfStepZ / 2.0 && world_z < surf_z + surfStepZ / 2.0) )	
    {
      if(Verbosity() > 2)
	std::cout <<  "     got it:  surf_phi " << surf_phi << " surf_z " << surf_z 
		  << " surfStepPhi/2 " << surfStepPhi/2.0 << " surfStepZ/2 " << surfStepZ/2.0  
		  << " world_phi " << world_phi << " world_z " << world_z 
		  << " rounded_nsurf "<< rounded_nsurf << " surf_index " << nsurf
		  << std::endl;        
      
      surf_index = nsurf;
      subsurfkey = nsurf;
    }    
  else
    {
      if(Verbosity() > 1)
	{
	  std::cout << PHWHERE 
		    << "Error: TPC surface index not defined, skipping cluster!" 
		    << std::endl;
	  std::cout << "     coordinates: " << world[0] << "  " << world[1] << "  " << world[2] 
		    << " radius " << sqrt(world[0]*world[0]+world[1]*world[1]) << std::endl;
	  std::cout << "     world_phi " << world_phi << " world_z " << world_z << std::endl;
	  std::cout << "     surf coords: " << surf_center[0] << "  " << surf_center[1] << "  " << surf_center[2] << std::endl;
	  std::cout << "     surf_phi " << surf_phi << " surf_z " << surf_z << std::endl; 
	  std::cout << " number of surfaces " << surf_vec.size() << std::endl;
	}
      return nullptr;
    }
  
  return surf_vec[surf_index];
  
}
