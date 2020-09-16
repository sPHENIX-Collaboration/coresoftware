/*!
 *  \file		PHActsTrkFitter.C
 *  \brief		Refit SvtxTracks with PHActs.
 *  \details	Refit SvtxTracks with PHActs.
 *  \author	        Tony Frawley <afrawley@fsu.edu>
 */

#include "PHActsTrkFitter.h"
#include "MakeActsGeometry.h"
#include "ActsTrack.h"
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
#include <phool/PHTimer.h>

#include <Acts/EventData/TrackParameters.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>
#include <Acts/Surfaces/PlaneSurface.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/EventData/MultiTrajectory.hpp>
#include <Acts/EventData/MultiTrajectoryHelpers.hpp>
#include <Acts/Utilities/Definitions.hpp>

#include <ActsExamples/EventData/Track.hpp>
#include <ActsExamples/Framework/AlgorithmContext.hpp>

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
  if(Verbosity() > 0)
    std::cout << "Setup PHActsTrkFitter" << std::endl;
  
  if(createNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
    return Fun4AllReturnCodes::ABORTEVENT;
  
  if (getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
    return Fun4AllReturnCodes::ABORTEVENT;
  
  fitCfg.fit = ActsExamples::TrkrClusterFittingAlgorithm::makeFitterFunction(
               m_tGeometry->tGeometry,
	       m_tGeometry->magField);

  if(m_timeAnalysis)
    {
      m_timeFile = new TFile(Name().c_str(), "RECREATE");
      h_eventTime = new TH1F("h_eventTime",";time [ms]",100000,0,10000);
      h_fitTime = new TH2F("h_fitTime",";p_{T} [GeV];time [ms]",80,0,40,100000,0,1000);
      h_updateTime = new TH1F("h_updateTime",";time [ms]",
			      100000,0,1000);
      
      h_rotTime = new TH1F("h_rotTime",";time [ms]",100000,0,1000);
      h_stateTime = new TH1F("h_stateTime",";time [ms]",
			     100000,0,1000);			     
    }		 
  
  if(Verbosity() > 0)
    std::cout << "Finish PHActsTrkFitter Setup" << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsTrkFitter::Process()
{
  auto startEventTime = high_resolution_clock::now();

  m_event++;

  auto logLevel = Acts::Logging::INFO;

  if (Verbosity() > 0)
  {
    std::cout << PHWHERE << "Events processed: " << m_event << std::endl;
    std::cout << "Start PHActsTrkFitter::process_event" << std::endl;
    logLevel = Acts::Logging::VERBOSE;
  }

  auto logger = Acts::getDefaultLogger("PHActsTrkFitter", logLevel);

  std::map<unsigned int, ActsTrack>::iterator trackIter;

  for (trackIter = m_actsProtoTracks->begin();
       trackIter != m_actsProtoTracks->end();
       ++trackIter)
  {
    ActsTrack track = trackIter->second;
    /// Can correlate with the SvtxTrackMap with the key
    const unsigned int trackKey = trackIter->first;

    std::vector<SourceLink> sourceLinks = track.getSourceLinks();
    ActsExamples::TrackParameters trackSeed = track.getTrackParams();
  
    /// Acts cares about the track covariance as it helps the KF
    /// know whether or not to trust the initial track seed or not.
    /// We reset it here to some loose values as it helps Acts improve
    /// the fitting. 
    /// If the covariance is too loose, it won't be able to propagate,
    /// but if it is too tight, it will just "believe" the track seed over
    /// the hit data
    Acts::BoundSymMatrix cov;
    cov << 1000 * Acts::UnitConstants::um, 0., 0., 0., 0., 0.,
           0., 1000 * Acts::UnitConstants::um, 0., 0., 0., 0.,
           0., 0., 0.05, 0., 0., 0.,
           0., 0., 0., 0.05, 0., 0.,
           0., 0., 0., 0., 0.00005 , 0.,
           0., 0., 0., 0., 0., 1.;

    ActsExamples::TrackParameters newTrackSeed(cov,
					       trackSeed.position(),
					       trackSeed.momentum(),
					       -1 *trackSeed.charge(),
					       trackSeed.time());

    /// Construct a perigee surface as the target surface
    /// This surface is what Acts fits with respect to, so we set it to
    /// the initial vertex estimation
    auto pSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(
		          track.getVertex());
   
    if(Verbosity() > 0)
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
      Acts::LoggerWrapper(*logger),
      &(*pSurface));

    auto startTime = high_resolution_clock::now();
    auto result = fitCfg.fit(sourceLinks, newTrackSeed, kfOptions);
    auto stopTime = high_resolution_clock::now();
    auto fitTime = duration_cast<microseconds>(stopTime - startTime);
 
    /// Check that the track fit result did not return an error
    if (result.ok())
    {  
      const FitResult& fitOutput = result.value();
   
      /// Make a trajectory state for storage, which conforms to Acts track fit
      /// analysis tool
      std::vector<size_t> trackTips;
      trackTips.push_back(fitOutput.trackTip);
      ActsExamples::IndexedParams indexedParams;
      if (fitOutput.fittedParameters)
      {
	indexedParams.emplace(fitOutput.trackTip, fitOutput.fittedParameters.value());
	const auto& params = fitOutput.fittedParameters.value();
        
	if(m_timeAnalysis)
	{
	  float px = params.momentum()(0);
	  float py = params.momentum()(1);
	  float pt = sqrt(px*px+py*py);
	  h_fitTime->Fill(pt,fitTime.count() / 1000.);
	}
        
	if (Verbosity() > 2)
        {
          std::cout << "Fitted parameters for track" << std::endl;
          std::cout << " position : " << params.position(m_tGeometry->geoContext).transpose()

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
      
      auto startUpdateTime = high_resolution_clock::now();
      if(fitOutput.fittedParameters)
	updateSvtxTrack(trajectory, trackKey, track.getVertex());

      auto stopUpdateTime = high_resolution_clock::now();
      auto updateTime = duration_cast<microseconds>
	(stopUpdateTime - startUpdateTime);

      if(m_timeAnalysis)
	h_updateTime->Fill(updateTime.count() / 1000.);

      if(Verbosity() > 4)
	std::cout << "Fit time == " << fitTime.count() / 1000. << " ms" 
		  << std::endl
		  << "Update SvtxTrack time == " << updateTime.count() / 1000. 
		  << " ms " << std::endl; 

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
				 (trackKey, ActsExamples::TrkrClusterMultiTrajectory()));

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
  
  
  auto stopEventTime = high_resolution_clock::now();
  if(m_timeAnalysis)
    {    
      auto eventTime = duration_cast<microseconds>(stopEventTime - startEventTime);
  
      h_eventTime->Fill(eventTime.count()/1000.);
    }

  if(Verbosity() > 0)
    std::cout << "PHActsTrkFitter::process_event finished" 
	      << std::endl;

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
      h_fitTime->Write();
      h_eventTime->Write();
      h_rotTime->Write();
      h_stateTime->Write();
      h_updateTime->Write();
      m_timeFile->Write();
      m_timeFile->Close();
    } 

  std::cout<<"The Acts track fitter had " << m_nBadFits <<" fits return an error"<<std::endl;

  if (Verbosity() > 0)
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
  track->set_x(params.position(m_tGeometry->geoContext)(0)
	       / Acts::UnitConstants::cm);
  track->set_y(params.position(m_tGeometry->geoContext)(1)
	       / Acts::UnitConstants::cm);
  track->set_z(params.position(m_tGeometry->geoContext)(2)
	       / Acts::UnitConstants::cm);

  track->set_px(params.momentum()(0));
  track->set_py(params.momentum()(1));
  track->set_pz(params.momentum()(2));
  
  track->set_chisq(trajState.chi2Sum);
  track->set_ndf(trajState.NDF);

  ActsTransformations *rotater = new ActsTransformations();
  rotater->setVerbosity(Verbosity());
  
  auto startRotTime = high_resolution_clock::now();
  
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
 
  auto stopRotTime = high_resolution_clock::now();
  auto rotTime = duration_cast<microseconds>(stopRotTime - startRotTime);

  if(m_timeAnalysis)
    h_rotTime->Fill(rotTime.count() / 1000.);


  // convert from mm to cm
  track->set_dca3d_xy(dca3Dxy / Acts::UnitConstants::cm);
  track->set_dca3d_z(dca3Dz / Acts::UnitConstants::cm);
  track->set_dca3d_xy_error(dca3DxyCov / Acts::UnitConstants::cm);
  track->set_dca3d_z_error(dca3DzCov / Acts::UnitConstants::cm);
  
  // Also need to update the state list and cluster ID list for all measurements associated with the acts track  
  // loop over acts track states, copy over to SvtxTrackStates, and add to SvtxTrack

  auto stateStartTime = high_resolution_clock::now();
  
  rotater->fillSvtxTrackStates(traj, trackTip, track,
			       m_tGeometry->geoContext,
			       m_hitIdClusKey);  

  auto stateStopTime = high_resolution_clock::now();
  auto stateTime = duration_cast<microseconds>(stateStopTime - stateStartTime);
  
  if(m_timeAnalysis)
    h_stateTime->Fill(stateTime.count() / 1000.);

  if(Verbosity() > 2)
    {  
      std::cout << " Identify fitted track after updating track states:" << std::endl;
      track->identify();
      std::cout << " cluster keys size " << track->size_cluster_keys() << std::endl;  
    }
 
 return;
  
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

