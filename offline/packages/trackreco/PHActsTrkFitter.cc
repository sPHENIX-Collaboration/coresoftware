/*!
 *  \file		PHActsTrkFitter.C
 *  \brief		Refit SvtxTracks with PHActs.
 *  \details	Refit SvtxTracks with PHActs.
 *  \author	        Tony Frawley <afrawley@fsu.edu>
 */

#include "PHActsTrkFitter.h"
#include "MakeActsGeometry.h"
#include "ActsTrack.h"
#include "ActsCovarianceRotater.h"

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
#include <Acts/Utilities/Definitions.hpp>

#include <ACTFW/EventData/Track.hpp>
#include <ACTFW/Framework/AlgorithmContext.hpp>

#include <cmath>
#include <iostream>
#include <vector>
#include <chrono>

using namespace std::chrono;

PHActsTrkFitter::PHActsTrkFitter(const std::string& name)
  : PHTrackFitting(name)
  , m_event(0)
  , m_actsFitResults(nullptr)
  , m_actsProtoTracks(nullptr)
  , m_tGeometry(nullptr)
  , m_trackMap(nullptr)
  , m_hitIdClusKey(nullptr)
  , m_nBadFits(0)
  , m_timeAnalysis(false)
  , m_timeFile(nullptr)
  , h_eventTime(nullptr)
  
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
  
  auto logger = Acts::Logging::INFO;
  if(Verbosity() > 2)
    logger = Acts::Logging::VERBOSE;

  fitCfg.fit = FW::TrkrClusterFittingAlgorithm::makeFitterFunction(
               m_tGeometry->tGeometry,
	       m_tGeometry->magField,
	       logger);

  if(m_timeAnalysis)
    {
      m_timeFile = new TFile("ActsTimeFile.root","RECREATE");
      h_eventTime = new TH1F("h_eventTime",";time [ms]",100,0,100);
    }		 
  
  if(Verbosity() > 1)
    std::cout << "Finish PHActsTrkFitter Setup" << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsTrkFitter::Process()
{
  auto startTime = high_resolution_clock::now();
  m_event++;

  if (Verbosity() > 1)
  {
    std::cout << PHWHERE << "Events processed: " << m_event << std::endl;
    std::cout << "Start PHActsTrkFitter::process_event" << std::endl;
  }

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
    
    /// Construct a perigee surface as the target surface
    /// This surface is what Acts fits with respect to, so we set it to
    /// the initial vertex estimation
    auto pSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(
		          track.getVertex());
   
    if(Verbosity() > 1)
      {
	std::cout << " Processing proto track with position:" 
		  << trackSeed.position() << std::endl 
		  << "momentum: " << trackSeed.momentum() << std::endl
		  << "charge : "<<trackSeed.charge() << std::endl
		  << "initial vertex : "<<track.getVertex()
		  << " corresponding to SvtxTrack key "<< trackKey
		  << std::endl;
	std::cout << "proto track covariance " << std::endl
		  << trackSeed.covariance().value() << std::endl;
     
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

    /// Check that the track fit result did not return an error
    if (result.ok())
    {  
      const FitResult& fitOutput = result.value();

      /// Make a trajectory state for storage, which conforms to Acts track fit
      /// analysis tool
      std::vector<size_t> trackTips;
      trackTips.push_back(fitOutput.trackTip);
      FW::IndexedParams indexedParams;
      if (fitOutput.fittedParameters)
      {
	indexedParams.emplace(fitOutput.trackTip, fitOutput.fittedParameters.value());

        if (Verbosity() > 2)
        {
	  const auto& params = fitOutput.fittedParameters.value();
          std::cout << "Fitted parameters for track" << std::endl;
          std::cout << " position : " << params.position().transpose()
                    << std::endl;
	  std::cout << "charge: "<<params.charge()<<std::endl;
          std::cout << " momentum : " << params.momentum().transpose()
                    << std::endl;
	  std::cout << "For trackTip == " << fitOutput.trackTip << std::endl;
        }
      }

      Trajectory trajectory(fitOutput.fittedStates, trackTips, indexedParams);

      /// Get position, momentum from the Acts output. Update the values of
      /// the proto track
      if(fitOutput.fittedParameters)
	updateSvtxTrack(trajectory, trackKey, track.getVertex());

      /// Insert a new entry into the map
      m_actsFitResults->insert(
	   std::pair<const unsigned int, Trajectory>(trackKey, trajectory));
  
    }
    else
      {
	if(Verbosity() > 10)
	  std::cout<<"Track fit failed"<<std::endl;
	/// Insert an empty track fit output into the map since the fit failed
       	m_actsFitResults->insert(std::pair<const unsigned int, Trajectory>
				 (trackKey, FW::TrkrClusterMultiTrajectory()));

	/// Mark the SvtxTrack as bad, for better analysis
	/// can remove later
	SvtxTrackMap::Iter trackIter = m_trackMap->find(trackKey);
	SvtxTrack *track = trackIter->second;

	track->set_x(-9999);
	track->set_y(-9999);
	track->set_z(-9999);
	track->set_px(-9999);
	track->set_py(-9999);
	track->set_pz(-9999);

	m_nBadFits++;
      }
    

  }

  auto stopTime = high_resolution_clock::now();
  auto eventTime = duration_cast<milliseconds>(stopTime - startTime);

  if(m_timeAnalysis)
    h_eventTime->Fill(eventTime.count());

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsTrkFitter::ResetEvent(PHCompositeNode *topNode)
{
  m_actsFitResults->clear();
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
      h_eventTime->Write();
      m_timeFile->Write();
      m_timeFile->Close();
    } 

  std::cout<<"The Acts track fitter had " << m_nBadFits <<" fits return an error"<<std::endl;

  if (Verbosity() > 10)
  {
    std::cout << "Finished PHActsTrkFitter" << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

void PHActsTrkFitter::updateSvtxTrack(Trajectory traj, 
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
  SvtxTrackState_v1 out( pathlength);
  out.set_x(0.0);
  out.set_y(0.0);
  out.set_z(0.0);
  track->insert_state(&out);   

  auto trajState =
    Acts::MultiTrajectoryHelpers::trajectoryState(mj, trackTip);
 
  const auto& params = traj.trackParameters(trackTip);

  /// Acts default unit is mm. So convert to cm
  track->set_x(params.position()(0) / Acts::UnitConstants::cm);
  track->set_y(params.position()(1) / Acts::UnitConstants::cm);
  track->set_z(params.position()(2) / Acts::UnitConstants::cm);

  track->set_px(params.momentum()(0));
  track->set_py(params.momentum()(1));
  track->set_pz(params.momentum()(2));
  
  track->set_chisq(trajState.chi2Sum);
  track->set_ndf(trajState.NDF);

  if(params.covariance())
    {
      ActsCovarianceRotater *rotater = new ActsCovarianceRotater();
      Acts::BoundSymMatrix rotatedCov = 
	rotater->rotateActsCovToSvtxTrack(params);
      
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

  calculateDCA(params, vertex, 
	       dca3Dxy, dca3Dz, dca3DxyCov, dca3DzCov);
 
  // convert from mm to cm
  track->set_dca3d_xy(dca3Dxy / Acts::UnitConstants::cm);
  track->set_dca3d_z(dca3Dz / Acts::UnitConstants::cm);
  track->set_dca3d_xy_error(dca3DxyCov / Acts::UnitConstants::cm);
  track->set_dca3d_z_error(dca3DzCov / Acts::UnitConstants::cm);
  
  // Also need to update the state list and cluster ID list for all measurements associated with the acts track  
  // loop over acts track states, copy over to SvtxTrackStates, and add to SvtxTrack
  fillSvtxTrackStates(traj, trackTip, track);  

  if(Verbosity() > 2)
    {  
      std::cout << " Identify fitted track after updating track states:" << std::endl;
      track->identify();
      std::cout << " cluster keys size " << track->size_cluster_keys() << std::endl;  
    }
 
 return;
  
}

void PHActsTrkFitter::fillSvtxTrackStates(const Trajectory traj, const size_t &trackTip, SvtxTrack *svtx_track)
{
  const auto &[trackTips, mj] = traj.trajectory();
  
  mj.visitBackwards(trackTip, [&](const auto &state) {
      /// Only fill the track states with non-outlier measurement
      auto typeFlags = state.typeFlags();
      if (not typeFlags.test(Acts::TrackStateFlag::MeasurementFlag))
	{
	  return true;
	}
      
      auto meas = std::get<Measurement>(*state.uncalibrated());

      /// Get the surface, if we need geometry information
      ///auto stateSurface = meas.referenceSurface();

      /// Get local position
      Acts::Vector2D local(meas.parameters()[Acts::ParDef::eLOC_0],
			   meas.parameters()[Acts::ParDef::eLOC_1]);
      /// Get global position
      Acts::Vector3D global(0, 0, 0);
      /// This is an arbitrary vector. Doesn't matter in coordinate transformation
      /// in Acts code
      Acts::Vector3D mom(1, 1, 1);
      meas.referenceObject().localToGlobal(m_tGeometry->geoContext,
					    local, mom, global);
      
      float pathlength = state.pathLength() / Acts::UnitConstants::cm;  
      SvtxTrackState_v1 out( pathlength );
      out.set_x(global.x() / Acts::UnitConstants::cm);
      out.set_y(global.y() / Acts::UnitConstants::cm);
      out.set_z(global.z() /  Acts::UnitConstants::cm);

      // I assume we want the smoothed for the final track states?      
      if (state.hasSmoothed())
	{
	  Acts::BoundParameters parameter(m_tGeometry->geoContext,
					  state.smoothedCovariance(), state.smoothed(),
					  state.referenceSurface().getSharedPtr());
	  
	  out.set_px(parameter.momentum().x());
	  out.set_py(parameter.momentum().y());
	  out.set_pz(parameter.momentum().z());

	  /// Get measurement covariance    
	  ActsCovarianceRotater *rotater = new ActsCovarianceRotater();
	  rotater->setVerbosity(Verbosity());

	  Acts::BoundSymMatrix globalCov = rotater->rotateActsCovToSvtxTrack(parameter);
	  for (int i = 0; i < 6; i++)
	    {
	      for (int j = 0; j < 6; j++)
		{ 
		  out.set_error(i, j, globalCov(i,j)); 
		}
	    }
	  	  
	  const unsigned int hitId = state.uncalibrated().hitID();
	  TrkrDefs::cluskey cluskey = getClusKey(hitId);
	  svtx_track->insert_cluster_key(cluskey);
	  
	  if(Verbosity() > 2)
	    {
	      std::cout << " inserting state with x,y,z = " << global.x() /  Acts::UnitConstants::cm 
			<< "  " << global.y() /  Acts::UnitConstants::cm << "  " 
			<< global.z() /  Acts::UnitConstants::cm 
			<< " pathlength " << pathlength
			<< " momentum px,py,pz = " <<  parameter.momentum().x() << "  " <<  parameter.momentum().y() << "  " << parameter.momentum().y()  
			<< " cluskey " << cluskey << std::endl
			<< "covariance " << globalCov << std::endl; 
	    }
	  
	  svtx_track->insert_state(&out);      
	}
  
      return true;      
    }
    );

  return;
}

TrkrDefs::cluskey PHActsTrkFitter::getClusKey(const unsigned int hitID)
{
  TrkrDefs::cluskey clusKey = 0;
  /// Unfortunately the map is backwards for looking up cluster key from
  /// hit ID. So we need to iterate over it. There won't be duplicates since
  /// the cluster key and hit id are a one-to-one map
  std::map<TrkrDefs::cluskey, unsigned int>::iterator
      hitIter = m_hitIdClusKey->begin();
  while (hitIter != m_hitIdClusKey->end())
  {
    if (hitIter->second == hitID)
    {
      clusKey = hitIter->first;
      break;
    }
    ++hitIter;
  }

  return clusKey;
}
    
void PHActsTrkFitter::calculateDCA(const Acts::BoundParameters param,
				   Acts::Vector3D vertex,
				   float &dca3Dxy,
				   float &dca3Dz,
				   float &dca3DxyCov,
				   float &dca3DzCov)
{
  Acts::Vector3D pos = param.position();
  Acts::Vector3D mom = param.momentum();

  /// Correct for initial vertex estimation
  pos -= vertex;

  Acts::BoundSymMatrix cov = Acts::BoundSymMatrix::Zero();
  if(param.covariance())
    cov = param.covariance().value();

  Acts::ActsSymMatrixD<3> posCov;
  for(int i = 0; i < 3; ++i)
    {
      for(int j = 0; j < 3; ++j)
	{
	  posCov(i,j) = cov(i,j);
	} 
    }

  Acts::Vector3D r = mom.cross(Acts::Vector3D(0.,0.,1.));
  float phi = atan2(r(1), r(0));

  Acts::RotationMatrix3D rot;
  Acts::RotationMatrix3D rot_T;
  rot(0,0) = cos(phi);
  rot(0,1) = -sin(phi);
  rot(0,2) = 0;
  rot(1,0) = sin(phi);
  rot(1,1) = cos(phi);
  rot(1,2) = 0;
  rot(2,0) = 0;
  rot(2,1) = 0;
  rot(2,2) = 1;
  
  rot_T = rot.transpose();

  Acts::Vector3D pos_R = rot * pos;
  Acts::ActsSymMatrixD<3> rotCov = rot * posCov * rot_T;

  dca3Dxy = pos_R(0);
  dca3Dz = pos_R(2);
  dca3DxyCov = rotCov(0,0);
  dca3DzCov = rotCov(2,2);
  
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

  m_actsFitResults = findNode::getClass<std::map<const unsigned int, Trajectory>>(topNode, "ActsTrajectories");
  
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

  m_hitIdClusKey = findNode::getClass<std::map<TrkrDefs::cluskey, unsigned int>>(topNode, "HitIDClusIDActsMap");
  
  if (!m_hitIdClusKey)
    {
      std::cout << PHWHERE << "No HitID:ClusKey map on node tree. Bailing."
		<< std::endl;
      
      return Fun4AllReturnCodes::EVENT_OK;
    }


  return Fun4AllReturnCodes::EVENT_OK;
}

