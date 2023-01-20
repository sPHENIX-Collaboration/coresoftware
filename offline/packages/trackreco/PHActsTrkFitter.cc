/*!
 *  \file		PHActsTrkFitter.C
 *  \brief		Refit SvtxTracks with PHActs.
 *  \details	Refit SvtxTracks with PHActs.
 *  \author	        Tony Frawley <afrawley@fsu.edu>
 */

#include "PHActsTrkFitter.h"

/// Tracking includes
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/MvtxDefs.h>
#include <trackbase/InttDefs.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/ClusterErrorPara.h>
#include <trackbase/Calibrator.h>

#include <trackbase_historic/ActsTransformations.h>
#include <trackbase_historic/SvtxTrack_v4.h>
#include <trackbase_historic/SvtxTrackState_v1.h>
#include <trackbase_historic/SvtxTrackMap_v2.h>
#include <trackbase_historic/TrackSeed.h>
#include <trackbase_historic/TrackSeedContainer.h>
#include <trackbase_historic/SvtxAlignmentStateMap_v1.h>

#include <g4detectors/PHG4TpcCylinderGeomContainer.h>

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

#include <tpc/TpcDistortionCorrectionContainer.h>

#include <Acts/EventData/TrackParameters.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>
#include <Acts/Surfaces/PlaneSurface.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/EventData/MultiTrajectory.hpp>
#include <Acts/EventData/MultiTrajectoryHelpers.hpp>
#include <Acts/TrackFitting/GainMatrixSmoother.hpp>
#include <Acts/TrackFitting/GainMatrixUpdater.hpp>


#include <TDatabasePDG.h>

#include <cmath>
#include <iostream>
#include <vector>

namespace
{
  // check vector validity
  inline bool is_valid( const Acts::Vector3 vec )
  {  return !( std::isnan( vec.x() ) || std::isnan( vec.y() ) || std::isnan( vec.z() ) ); }  
}

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <trackbase/alignmentTransformationContainer.h>


PHActsTrkFitter::PHActsTrkFitter(const std::string& name)
  : SubsysReco(name)
  , m_trajectories(nullptr)
{}

int PHActsTrkFitter::InitRun(PHCompositeNode* topNode)
{
  if(Verbosity() > 1)
    { std::cout << "Setup PHActsTrkFitter" << std::endl; }

  if(createNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
    { return Fun4AllReturnCodes::ABORTEVENT; }

  if (getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
    { return Fun4AllReturnCodes::ABORTEVENT; }
  
  m_alignStates.distortionContainers(_dcc_static, _dcc_average, _dcc_fluctuation);
  m_alignStates.actsGeometry(m_tGeometry);
  m_alignStates.clusters(m_clusterContainer);
  m_alignStates.stateMap(m_alignmentStateMap);
  m_alignStates.verbosity(Verbosity());

  m_fitCfg.fit = ActsTrackFittingAlgorithm::makeKalmanFitterFunction(
    m_tGeometry->geometry().tGeometry,
    m_tGeometry->geometry().magField);

  m_fitCfg.dFit = ActsTrackFittingAlgorithm::makeKalmanFitterFunction(m_tGeometry->geometry().magField);

  m_outlierFinder.verbosity = Verbosity();
  std::map<long unsigned int, float> chi2Cuts;
  chi2Cuts.insert(std::make_pair(10,4));
  chi2Cuts.insert(std::make_pair(12,4));
  chi2Cuts.insert(std::make_pair(14,9));
  chi2Cuts.insert(std::make_pair(16,4));
  m_outlierFinder.chi2Cuts = chi2Cuts;
  if(m_useOutlierFinder)
    {
      m_fitCfg.fit->outlierFinder(m_outlierFinder);
    }

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

   auto cellgeo =
      findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");

  if (cellgeo)
  {
    _clusterMover.initialize_geometry(cellgeo);
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
    {
      std::cout << " SvtxTrackMap size is now " << m_trackMap->size() 
	      << std::endl;

    }

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

  if(Verbosity()>0)
    {
      std::cout << " seed map size " << m_seedMap->size() << std::endl;
    }

  for(auto trackiter = m_seedMap->begin(); trackiter != m_seedMap->end();
      ++trackiter)
    {
      TrackSeed *track = *trackiter;
      if(!track)
	{ continue; }
      
      unsigned int tpcid = track->get_tpc_seed_index();
      unsigned int siid = track->get_silicon_seed_index();
      
      if(Verbosity() >1)
	{ std::cout << "tpc and si id " << tpcid << ", " << siid << std::endl; }

      /// A track seed is made for every tpc seed. Not every tpc seed
      /// has a silicon match, we skip those cases completely in pp running
      if(m_pp_mode && siid == std::numeric_limits<unsigned int>::max()) 
	{
	  if(Verbosity() > 1) std::cout << "Running in pp mode and SvtxSeedTrack has no silicon match, skip it" << std::endl;
	  continue;
	}

      // get the crossing number
      auto siseed = m_siliconSeeds->get(siid);
      short crossing = SHRT_MAX;
      if(siseed) crossing = siseed->get_crossing();
      else if(!m_pp_mode) crossing = 0;

      // if the crossing was not determined in pp running, skip this case completely
      if(m_pp_mode && crossing == SHRT_MAX) 
	{
	  // Skip this in the pp case.
	  if(Verbosity() > 1) std::cout << "Crossing not determined, skipping track" << std::endl;
	  continue;
	}

      auto tpcseed = m_tpcSeeds->get(tpcid);

      /// Need to also check that the tpc seed wasn't removed by the ghost finder
      if(!tpcseed)
	{ std::cout << "no tpc seed"<<std::endl; continue; }

      if(Verbosity() > 0) 
	{
	  if(siseed) std::cout << " silicon seed position is (x,y,z) = " << siseed->get_x() << "  " << siseed->get_y() << "  " << siseed->get_z() << std::endl;
	  std::cout << " tpc seed position is (x,y,z) = " << tpcseed->get_x() << "  " << tpcseed->get_y() << "  " << tpcseed->get_z() << std::endl;
	}

      PHTimer trackTimer("TrackTimer");
      trackTimer.stop();
      trackTimer.restart();
      ActsTrackFittingAlgorithm::MeasurementContainer measurements;
  
      SourceLinkVec sourceLinks;
      if(siseed) sourceLinks = getSourceLinks(siseed, measurements, crossing);
      const auto tpcSourceLinks = getSourceLinks(tpcseed, measurements, crossing);
      sourceLinks.insert( sourceLinks.end(), tpcSourceLinks.begin(), tpcSourceLinks.end() );

      // position comes from the silicon seed, unless there is no silicon seed
      Acts::Vector3 position(0,0,0);
      if(siseed)
        {
          position(0) = siseed->get_x() * Acts::UnitConstants::cm;
          position(1) = siseed->get_y() * Acts::UnitConstants::cm;
          position(2) = siseed->get_z() * Acts::UnitConstants::cm;
        }
      else
        {
          position(0) = tpcseed->get_x() * Acts::UnitConstants::cm;
          position(1) = tpcseed->get_y() * Acts::UnitConstants::cm;
          position(2) = tpcseed->get_z() * Acts::UnitConstants::cm;
        }
      if( !is_valid( position ) ) continue;

      if(sourceLinks.empty()) { continue; }

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
          { return m_tGeometry->maps().isMicromegasSurface( surface ); } ) )
        { continue; }
      }

      Acts::Vector3 momentum(
	       tpcseed->get_px(m_clusterContainer, m_tGeometry), 
	       tpcseed->get_py(m_clusterContainer, m_tGeometry),
	       tpcseed->get_pz());
      if( !is_valid( momentum ) ) continue;
 
      auto pSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(
					  position);
      auto actsFourPos = Acts::Vector4(position(0), position(1),
				       position(2),
				       10 * Acts::UnitConstants::ns);
      Acts::BoundSymMatrix cov = setDefaultCovariance();
 
      int charge = tpcseed->get_charge();
      if(m_fieldMap.find("3d") != std::string::npos)
	{ charge *= -1; }

      /// Acts requires a wrapped vector, so we need to replace the
      /// std::vector contents with a wrapper vector to get the memory
      /// access correct
      std::vector<std::reference_wrapper<const SourceLink>>  wrappedSls;
      for(const auto& sl : sourceLinks)
      { wrappedSls.push_back(std::cref(sl)); }
           
      /// Reset the track seed with the dummy covariance
      auto seed = ActsTrackFittingAlgorithm::TrackParameters::create(
        pSurface,
        m_tGeometry->geometry().getGeoContext(),
        actsFourPos,
        momentum,
        charge / momentum.norm(),
        cov).value();
      
      if(Verbosity() > 2)
      { printTrackSeed(seed); }

      /// Set host of propagator options for Acts to do e.g. material integration
      Acts::PropagatorPlainOptions ppPlainOptions;
      ppPlainOptions.absPdgCode = m_pHypothesis;
      ppPlainOptions.mass = TDatabasePDG::Instance()->GetParticle(
        m_pHypothesis)->Mass() * Acts::UnitConstants::GeV;
       
      Calibrator calibrator{measurements};

      auto magcontext = m_tGeometry->geometry().magFieldContext;
      auto calibcontext = m_tGeometry->geometry().calibContext;

      ActsTrackFittingAlgorithm::GeneralFitterOptions 
        kfOptions{
	m_tGeometry->geometry().getGeoContext(),
        magcontext,
        calibcontext,
        calibrator,
        &(*pSurface),
        Acts::LoggerWrapper(*logger),ppPlainOptions};
      
      PHTimer fitTimer("FitTimer");
      fitTimer.stop();
      fitTimer.restart();
      auto mtj = std::make_shared<Acts::VectorMultiTrajectory>();
      auto result = fitTrack(wrappedSls, seed, kfOptions, surfaces,mtj);
      fitTimer.stop();
      auto fitTime = fitTimer.get_accumulated_time();
   
      if(Verbosity() > 1)
	{ std::cout << "PHActsTrkFitter Acts fit time " << fitTime << std::endl; }

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
	  
        SvtxTrack_v4 newTrack;
        newTrack.set_tpc_seed(tpcseed);
        newTrack.set_crossing(crossing);
        newTrack.set_silicon_seed(siseed);

        if(m_fitSiliconMMs)
        {
          
          unsigned int trid = m_directedTrackMap->size();
          newTrack.set_id(trid);

          if( getTrackFitResult(fitOutput, &newTrack) )
          { m_directedTrackMap->insertWithKey(&newTrack, trid); }
          
        } else {
          
          unsigned int trid = m_trackMap->size();
          newTrack.set_id(trid);

          if( getTrackFitResult(fitOutput, &newTrack))
          { m_trackMap->insertWithKey(&newTrack, trid); }
        
        }
        
      } else if (!m_fitSiliconMMs) {

        /// Track fit failed, get rid of the track from the map
        m_nBadFits++;
        if(Verbosity() > 1)
        { 
          std::cout << "Track fit failed for track " << m_seedMap->find(track) 
            << " with Acts error message " 
            << result.error() << ", " << result.error().message()
            << std::endl;
        }
      }

      trackTimer.stop();
      auto trackTime = trackTimer.get_accumulated_time();
      
      if(Verbosity() > 1)
	{ std::cout << "PHActsTrkFitter total single track time " << trackTime << std::endl; }
    
    }

  return;

}

//___________________________________________________________________________________
SourceLinkVec PHActsTrkFitter::getSourceLinks(TrackSeed* track,
					      ActsTrackFittingAlgorithm::MeasurementContainer& measurements,
				   short int crossing )
{

  SourceLinkVec sourcelinks;

  if(m_pp_mode && crossing == SHRT_MAX) 
    {
      // Need to skip this in the pp case, for AuAu it should not happen
      return sourcelinks; 
    }

  PHTimer SLTrackTimer("SLTrackTimer");
  SLTrackTimer.stop();
  SLTrackTimer.restart();
    
  // loop over all clusters
  std::vector<std::pair<TrkrDefs::cluskey, Acts::Vector3>> global_raw;

  for (auto clusIter = track->begin_cluster_keys();
       clusIter != track->end_cluster_keys();
       ++clusIter)
    {
      auto key = *clusIter;
      auto cluster = m_clusterContainer->findCluster(key);
      if(!cluster)
	{
	  if(Verbosity() > 0) std::cout << "Failed to get cluster with key " << key << " for track " << m_seedMap->find(track) << std::endl;
    else std::cout<< "PHActsTrkFitter :: Key: "<< key << " for track " << m_seedMap->find(track) <<std::endl;
	  continue;
	}

      auto subsurfkey = cluster->getSubSurfKey();
      
      /// Make a safety check for clusters that couldn't be attached
      /// to a surface
      auto surf = m_tGeometry->maps().getSurface(key, cluster);
      if(!surf)
	{ continue; }

      const unsigned int trkrid = TrkrDefs::getTrkrId(key);
      const unsigned int side = TpcDefs::getSide(key);

      // For the TPC, cluster z has to be corrected for the crossing z offset, distortion, and TOF z offset 
      // we do this locally here and do not modify the cluster, since the cluster may be associated with multiple silicon tracks  
      Acts::Vector3 global  = m_tGeometry->getGlobalPosition(key, cluster);

      // temporary for testing transforms 
      //=========================
      bool test_transforms = false;

      if(test_transforms)  
	{

	  // Alignment transformation testing purposes
	  auto hitsetkey = TrkrDefs::getHitSetKeyFromClusKey(key);

	  float globphi = atan2(global(1),global(0))*180.0/M_PI;
	  std::cout << "Check in TrkFitter: global phi " << globphi << " hitsetkey: " << hitsetkey <<" global: " << std::endl << global << std::endl;

	  auto x = cluster->getLocalX() * 10.0;   // mm
	  auto y = cluster->getLocalY() * 10.0;
	  
	  if(trkrid == TrkrDefs::tpcId)
	    {
	      // must convert local Y from cluster average time of arival to local cluster z position
	      double drift_velocity = m_tGeometry->get_drift_velocity();
	      double zdriftlength = cluster->getLocalY() * drift_velocity;
	      double surfCenterZ = 52.89; // 52.89 is where G4 thinks the surface center is
	      double zloc = surfCenterZ - zdriftlength;   // converts z drift length to local z position in the TPC in north
	      unsigned int side = TpcDefs::getSide(key);
	      if(side == 0) zloc = -zloc;
	      y = zloc * 10.0;
	    }

	  Eigen::Vector3d clusterLocalPosition (x,y,0);  // follows the convention for the acts transform of local = (x,z,y)
	  std::cout << "local: "<< std::endl <<clusterLocalPosition << std::endl;

	  if (trkrid == TrkrDefs::inttId)
	    {
	      unsigned int layer     = TrkrDefs::getLayer(hitsetkey);
	      unsigned int ladderz   = InttDefs::getLadderZId(hitsetkey);
	      unsigned int ladderphi = InttDefs::getLadderPhiId(hitsetkey);
	      std::cout << "layer: "<<layer<< " ladderZ: "<<ladderz<< " ladderPhi: " << ladderphi<<std::endl;
	    }
	  else if (trkrid == TrkrDefs::mvtxId)
	    {
	      unsigned int layer                          = TrkrDefs::getLayer(hitsetkey);
	      unsigned int stave                          = MvtxDefs::getStaveId(hitsetkey);
	      unsigned int chip                           = MvtxDefs::getChipId(hitsetkey);
	      std::cout << "layer: " << layer << " stave: " << stave << "chip: " << chip << std::endl;
	    }
	  else if(trkrid == TrkrDefs::tpcId)
	    {
	      unsigned int layer                          = TrkrDefs::getLayer(hitsetkey);
	      unsigned int sector                         = TpcDefs::getSectorId(hitsetkey);
	      unsigned int side                           = TpcDefs::getSide(hitsetkey);
	      std::cout<< "subsurfkey: "<< subsurfkey << " layer: " << layer << " sector: " << sector 
		       << " side: " << side << std::endl;
	    }
	  else if(trkrid == TrkrDefs::micromegasId)
	    {
	      unsigned int layer                            = TrkrDefs::getLayer(hitsetkey);
	      unsigned short segmentation = (unsigned short) MicromegasDefs::getSegmentationType(hitsetkey);
	      unsigned int tile                             = MicromegasDefs::getTileId(hitsetkey);
	      std::cout<< " layer: " << layer << " segmentation: "<< segmentation  << " tile: " << tile << std::endl;
	    }

          Acts::GeometryIdentifier id = surf->geometryId();
	  std::cout << " Geometry Id: " << id << std::endl;

	  auto alignmentTransformation = m_alignmentTransformationMap->getTransform(id);      

	  std::cout << " Transform: " << std::endl << alignmentTransformation.matrix() << std::endl;

	  Eigen::Vector3d finalCoords = alignmentTransformation*clusterLocalPosition;
	  float phi = atan2(finalCoords(1),finalCoords(0))*180.0/M_PI;

	  finalCoords /= 10.0;
	  float deltaX = finalCoords(0)-global(0);
	  float deltaY = finalCoords(1)-global(1);

	  std::cout<< "deltax: "<<deltaX << " deltaY: " << deltaY << std::endl;
	  std::cout << " phi: "<< phi <<" Final Alignment Transform Coordinates: " << finalCoords << std::endl << std::endl;

	}  // end testing transforms
      //=========================

      if(trkrid ==  TrkrDefs::tpcId)
	{	  
	  // make all corrections to global position of TPC cluster
	  float z = m_clusterCrossingCorrection.correctZ(global[2], side, crossing);
	  global[2] = z;
	  
	  // apply distortion corrections
	  if(_dcc_static) { global = _distortionCorrection.get_corrected_position( global, _dcc_static ); }
	  if(_dcc_average) { global = _distortionCorrection.get_corrected_position( global, _dcc_average ); }
	  if(_dcc_fluctuation) { global = _distortionCorrection.get_corrected_position( global, _dcc_fluctuation ); }
	}
   
      if(Verbosity() > 0)
	{
	  std::cout << " zinit " << global[2] << " xinit " << global[0] << " yinit " << global[1] << " side " << side << " crossing " << crossing 
		    << " cluskey " << key << " subsurfkey " << subsurfkey << std::endl;
	}

      // add the global positions to a vector to give to the cluster mover
      global_raw.push_back(std::make_pair(key, global));
      
    }	  // end loop over clusters here
  
  // move the cluster positions back to the original readout surface
  auto global_moved = _clusterMover.processTrack(global_raw);
  
  // loop over global positions returned by cluster mover
  for(int i=0; i<global_moved.size(); ++i)
    {
      TrkrDefs::cluskey cluskey = global_moved[i].first;
      Acts::Vector3 global = global_moved[i].second;
   
      auto cluster = m_clusterContainer->findCluster(cluskey);
      Surface surf = m_tGeometry->maps().getSurface(cluskey, cluster);

      // if this is a TPC cluster, the crossing correction may have moved it across the central membrane, check the surface
      auto trkrid = TrkrDefs::getTrkrId(cluskey);
      if(trkrid == TrkrDefs::tpcId)
	{
	  TrkrDefs::hitsetkey hitsetkey = TrkrDefs::getHitSetKeyFromClusKey(cluskey);
	  TrkrDefs::subsurfkey new_subsurfkey = 0;    
	  surf = m_tGeometry->get_tpc_surface_from_coords(hitsetkey,  global, new_subsurfkey);
	}
    
      if(!surf)
	{ continue; }

      // get local coordinates
      Acts::Vector2 localPos;
      global *= Acts::UnitConstants::cm;

      Acts::Vector3 normal = surf->normal(m_tGeometry->geometry().getGeoContext());
      auto local = surf->globalToLocal(m_tGeometry->geometry().getGeoContext(),
				       global, normal);
     
      if(local.ok())
	{ localPos = local.value() / Acts::UnitConstants::cm; }
      else
	{
	  /// otherwise take the manual calculation for the TPC
	  Acts::Vector3 loct = surf->transform(m_tGeometry->geometry().getGeoContext()).inverse() * global;
	  loct /= Acts::UnitConstants::cm;

	  localPos(0) = loct(0);
	  localPos(1) = loct(1);
	}
      
      if(Verbosity() > 0)
	{
	  std::cout << " cluster global after mover: " << global << std::endl; 
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

      if(m_cluster_version==3){
	cov(Acts::eBoundLoc0, Acts::eBoundLoc0) = 
	  cluster->getActsLocalError(0,0) * Acts::UnitConstants::cm2;
	cov(Acts::eBoundLoc0, Acts::eBoundLoc1) =
	  cluster->getActsLocalError(0,1) * Acts::UnitConstants::cm2;
	cov(Acts::eBoundLoc1, Acts::eBoundLoc0) = 
	  cluster->getActsLocalError(1,0) * Acts::UnitConstants::cm2;
	cov(Acts::eBoundLoc1, Acts::eBoundLoc1) = 
	  cluster->getActsLocalError(1,1) * Acts::UnitConstants::cm2;
      }else if(m_cluster_version==4){
	double clusRadius = sqrt(global[0]*global[0] + global[1]*global[1]);
	auto para_errors = _ClusErrPara.get_cluster_error(track,cluster,clusRadius,cluskey);
	cov(Acts::eBoundLoc0, Acts::eBoundLoc0) = para_errors.first * Acts::UnitConstants::cm2;
	cov(Acts::eBoundLoc0, Acts::eBoundLoc1) = 0;
	cov(Acts::eBoundLoc1, Acts::eBoundLoc0) = 0;
	cov(Acts::eBoundLoc1, Acts::eBoundLoc1) = para_errors.second * Acts::UnitConstants::cm2;
      }

      ActsSourceLink::Index index = measurements.size();
      
      SourceLink sl(surf->geometryId(), index, cluskey);
      
      Acts::Measurement<Acts::BoundIndices,2> meas(sl, indices, loc, cov);
      if(Verbosity() > 3)
	{
	  std::cout << "source link " << sl.index() << ", loc : " 
		    << loc.transpose() << std::endl 
		    << ", cov : " << cov.transpose() << std::endl
		    << " geo id " << sl.geometryId() << std::endl;
	  std::cout << "Surface : " << std::endl;
	  surf.get()->toStream(m_tGeometry->geometry().getGeoContext(), std::cout);
	  std::cout << std::endl;
	  std::cout << "Cluster error " << cluster->getRPhiError() << " , " << cluster->getZError() << std::endl;
	  std::cout << "For key " << cluskey << " with local pos " << std::endl
		    << localPos(0) << ", " << localPos(1)
		    << std::endl;
	}
      
      sourcelinks.push_back(sl);
      measurements.push_back(meas);
 
    }
  
  SLTrackTimer.stop();
  auto SLTime = SLTrackTimer.get_accumulated_time();
  
  if(Verbosity() > 1)
    std::cout << "PHActsTrkFitter Source Links generation time:  "
	      << SLTime << std::endl;


  return sourcelinks;
}

bool PHActsTrkFitter::getTrackFitResult(const FitResult &fitOutput, SvtxTrack* track)
{
  /// Make a trajectory state for storage, which conforms to Acts track fit
  /// analysis tool
  std::vector<Acts::MultiTrajectoryTraits::IndexType> trackTips;
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
          std::cout << " position : " << params.position(m_tGeometry->geometry().getGeoContext()).transpose()
	    
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
      if(Verbosity() > 2)
	{ std::cout << " track fit failed for track " << track->get_id() << std::endl; }
	
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
  
  if(m_commissioning)
    {
      if(track->get_silicon_seed() && track->get_tpc_seed())
	{
	  m_alignStates.fillAlignmentStateMap(trajectory, track);
	}
    }
    
  updateTrackTimer.stop();
  auto updateTime = updateTrackTimer.get_accumulated_time();
  
  if(Verbosity() > 1)
    std::cout << "PHActsTrkFitter update SvtxTrack time "
	      << updateTime << std::endl;

  if(m_timeAnalysis)
    h_updateTime->Fill(updateTime);
  
  return true;
}

ActsTrackFittingAlgorithm::TrackFitterResult PHActsTrkFitter::fitTrack(
    const std::vector<std::reference_wrapper<const SourceLink>>& sourceLinks, 
    const ActsTrackFittingAlgorithm::TrackParameters& seed,
    const ActsTrackFittingAlgorithm::GeneralFitterOptions& kfOptions, 
    const SurfacePtrVec& surfSequence,
    std::shared_ptr<Acts::VectorMultiTrajectory>& mtj)
{
  if(m_fitSiliconMMs) 
  { 
    return (*m_fitCfg.dFit)(sourceLinks, seed, kfOptions, surfSequence, mtj); 
  } else {
    return (*m_fitCfg.fit)(sourceLinks, seed, kfOptions, mtj); 
  }
}

SourceLinkVec PHActsTrkFitter::getSurfaceVector(const SourceLinkVec& sourceLinks,
						SurfacePtrVec& surfaces) const
{
  SourceLinkVec siliconMMSls;

//   if(Verbosity() > 1)
//     std::cout << "Sorting " << sourceLinks.size() << " SLs" << std::endl;
  
  for(const auto& sl : sourceLinks)
  {
    if(Verbosity() > 1)
    { std::cout << "SL available on : " << sl.geometryId() << std::endl; }
      
    const auto surf = m_tGeometry->geometry().tGeometry->findSurface(sl.geometryId());
    // skip TPC surfaces
    if( m_tGeometry->maps().isTpcSurface( surf ) ) continue;
    
    // also skip micromegas surfaces if not used
    if( m_tGeometry->maps().isMicromegasSurface( surf ) && !m_useMicromegas ) continue;
    
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
    const auto& surface = surfaces.at(i);
    const auto thisVolume = surface->geometryId().volume();
    const auto thisLayer  = surface->geometryId().layer();
      
    const auto nextSurface = surfaces.at(i+1);
    const auto nextVolume = nextSurface->geometryId().volume();
    const auto nextLayer = nextSurface->geometryId().layer();
    
    /// Implement a check to ensure surfaces are sorted
    if(nextVolume == thisVolume) 
    {
      if(nextLayer < thisLayer)
      {
        std::cout 
          << "PHActsTrkFitter::checkSurfaceVec - " 
          << "Surface not in order... removing surface" 
          << surface->geometryId() << std::endl;
        
        surfaces.erase(surfaces.begin() + i);
	      
        /// Subtract one so we don't skip a surface
	      --i;
	      continue;
	    }
      
    } else {

      if(nextVolume < thisVolume)
      {
        std::cout 
          << "PHActsTrkFitter::checkSurfaceVec - " 
          << "Volume not in order... removing surface" 
          << surface->geometryId() << std::endl;
        
        surfaces.erase(surfaces.begin() + i);

        /// Subtract one so we don't skip a surface
	      --i;
	      continue;
      
      }
    }
  } 

}

void PHActsTrkFitter::updateSvtxTrack(Trajectory traj, SvtxTrack* track)
{
  const auto& mj = traj.multiTrajectory();
  const auto& tips = traj.tips();

  /// only one track tip in the track fit Trajectory
  auto &trackTip = tips.front();

  if(Verbosity() > 2)
    {
      std::cout << "Identify (proto) track before updating with acts results " << std::endl;
      track->identify();      
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
  track->set_x(params.position(m_tGeometry->geometry().getGeoContext())(0)
	       / Acts::UnitConstants::cm);
  track->set_y(params.position(m_tGeometry->geometry().getGeoContext())(1)
	       / Acts::UnitConstants::cm);
  track->set_z(params.position(m_tGeometry->geometry().getGeoContext())(2)
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
	    { track->set_error(i,j, rotatedCov(i,j)); }
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
				  m_tGeometry->geometry().getGeoContext());  
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

void PHActsTrkFitter::printTrackSeed(const ActsTrackFittingAlgorithm::TrackParameters& seed) const
{
  std::cout 
    << PHWHERE 
    << " Processing proto track:"
    << std::endl;  

  std::cout 
    << "position: " << seed.position(m_tGeometry->geometry().getGeoContext()).transpose() 
    << std::endl
    << "momentum: " << seed.momentum().transpose()
    << std::endl;

  std::cout << "charge : " << seed.charge() << std::endl;
  std::cout << "absolutemom : " << seed.absoluteMomentum() << std::endl;
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
  
  PHNodeIterator dstIter(topNode);
  PHCompositeNode *svtxNode = dynamic_cast<PHCompositeNode *>(dstIter.findFirst("PHCompositeNode", "SVTX"));

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
	  m_directedTrackMap = new SvtxTrackMap_v2;

	  PHIODataNode<PHObject> *trackNode = 
	    new PHIODataNode<PHObject>(m_directedTrackMap,"SvtxSiliconMMTrackMap","PHObject");
	  svtxNode->addNode(trackNode);
	} 
    }

  m_trajectories = findNode::getClass<std::map<const unsigned int, Trajectory>>(topNode, "ActsTrajectories");
  if(!m_trajectories)
    {
      m_trajectories = new std::map<const unsigned int, Trajectory>;
      auto node = 
	new PHDataNode<std::map<const unsigned int, Trajectory>>(m_trajectories, "ActsTrajectories");
      svtxNode->addNode(node);
      
    }
  
   m_trackMap = findNode::getClass<SvtxTrackMap>(topNode, _track_map_name);
  
  if(!m_trackMap)
    {
      m_trackMap = new SvtxTrackMap_v2;
      PHIODataNode<PHObject>* node = new PHIODataNode<PHObject>(m_trackMap,_track_map_name,"PHObject");
      svtxNode->addNode(node);
    }

  m_alignmentStateMap = findNode::getClass<SvtxAlignmentStateMap>(topNode, "SvtxAlignmentStateMap");
  if(!m_alignmentStateMap)
    {
      m_alignmentStateMap = new SvtxAlignmentStateMap_v1;
      auto node = new PHDataNode<SvtxAlignmentStateMap>(m_alignmentStateMap,"SvtxAlignmentStateMap","PHObject");
      svtxNode->addNode(node);
    }

  if(m_actsEvaluator)
    {
      m_seedTracks = findNode::getClass<SvtxTrackMap>(topNode,_seed_track_map_name);
      
      if(!m_seedTracks)
	{
	  m_seedTracks = new SvtxTrackMap_v2;
	  
	  PHIODataNode<PHObject> *seedNode = 
	    new PHIODataNode<PHObject>(m_seedTracks,_seed_track_map_name,"PHObject");
	  svtxNode->addNode(seedNode);
	}
    }
  
  return Fun4AllReturnCodes::EVENT_OK;
}




int PHActsTrkFitter::getNodes(PHCompositeNode* topNode)
{
  /*
  m_alignmentTransformationMap = findNode::getClass<alignmentTransformationContainer>(topNode, "alignmentTransformationContainer");
  if(!m_alignmentTransformationMap)
    {
      std::cout << PHWHERE << "alignmentTransformationContainer not on node tree. Bailing"
		<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  */



  m_tpcSeeds = findNode::getClass<TrackSeedContainer>(topNode, "TpcTrackSeedContainer");
  if(!m_tpcSeeds)
    {
      std::cout << PHWHERE << "TpcTrackSeedContainer not on node tree. Bailing"
		<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

    m_siliconSeeds = findNode::getClass<TrackSeedContainer>(topNode, "SiliconTrackSeedContainer");
  if(!m_siliconSeeds)
    {
      std::cout << PHWHERE << "SiliconTrackSeedContainer not on node tree. Bailing"
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

  m_tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if(!m_tGeometry)
    {
      std::cout << "ActsGeometry not on node tree. Exiting."
		<< std::endl;
      
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  m_seedMap = findNode::getClass<TrackSeedContainer>(topNode,"SvtxTrackSeedContainer");
  if(!m_seedMap)
    {
      std::cout << "No Svtx seed map on node tree. Exiting."
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
