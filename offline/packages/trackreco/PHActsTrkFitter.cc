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

#include <ActsExamples/EventData/Track.hpp>
#include <ActsExamples/Framework/AlgorithmContext.hpp>

#include <cmath>
#include <iostream>
#include <vector>

PHActsTrkFitter::PHActsTrkFitter(const std::string& name)
  : PHTrackFitting(name)
  , m_event(0)
  , m_actsFitResults(nullptr)
  , m_actsProtoTracks(nullptr)
  , m_tGeometry(nullptr)
  , m_trackMap(nullptr)
  , m_hitIdClusKey(nullptr)
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

  std::map<unsigned int, ActsTrack>::iterator trackIter;
 
  for (trackIter = m_actsProtoTracks->begin();
       trackIter != m_actsProtoTracks->end();
       ++trackIter)
  {
    auto trackTimer = std::make_unique<PHTimer>("TrackTimer");
    trackTimer->stop();
    trackTimer->restart();

    ActsTrack track = trackIter->second;
    /// Can correlate with the SvtxTrackMap with the key
    const unsigned int trackKey = trackIter->first;

    std::vector<SourceLink> sourceLinks = track.getSourceLinks();
    ActsExamples::TrackParameters trackSeed = track.getTrackParams();

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

    Acts::BoundSymMatrix cov = setDefaultCovariance();

    ActsExamples::TrackParameters newTrackSeed(
                  trackSeed.fourPosition(m_tGeometry->geoContext),
		  trackSeed.momentum(),
		  trackSeed.absoluteMomentum(),
		  trackSeed.charge(),
		  cov);

    /// Construct a perigee surface as the target surface
    /// This surface is what Acts fits with respect to, so we set it to
    /// the initial vertex estimation
    auto pSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(
		          track.getVertex());
   
    if(Verbosity() > 2)
      {
	std::cout << " Processing proto track with position:" 
		  << newTrackSeed.position(m_tGeometry->geoContext) 
		  << std::endl 
		  << "momentum: " << newTrackSeed.momentum() 
		  << std::endl
		  << "charge : " << newTrackSeed.charge() 
		  << std::endl
		  << "initial vertex : " << track.getVertex()
		  << " corresponding to SvtxTrack key " << trackKey
		  << std::endl;
	std::cout << "proto track covariance " << std::endl
		  << newTrackSeed.covariance().value() << std::endl;
     
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
      Acts::PropagatorPlainOptions(),
      &(*pSurface));

    auto fitTimer = std::make_unique<PHTimer>("FitTimer");
    fitTimer->stop();
    fitTimer->restart();
    auto result = fitTrack(sourceLinks, newTrackSeed, kfOptions,
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
   
      if(m_fitSiliconMMs)
	updateActsProtoTrack(fitOutput, trackIter);
      else
	getTrackFitResult(fitOutput, trackKey, track.getVertex());
	
    }
    else
      {
	if(Verbosity() > 10)
	  std::cout<<"Track fit failed"<<std::endl;
	/// Insert an empty track fit output into the map since the fit failed
       	m_actsFitResults->insert(
			  std::pair<const unsigned int, Trajectory>
			  (trackKey, ActsExamples::TrkrClusterMultiTrajectory()));

	// fit failed, delete the junk SvtxTrack from the node tree 
	// so the evaluator does not waste time on it
	m_trackMap->erase(trackKey);

	m_nBadFits++;
      }

    trackTimer->stop();
    auto trackTime = trackTimer->get_accumulated_time();
    
    if(Verbosity() > 0)
      std::cout << "PHActsTrkFitter total single track time "
		<< trackTime << std::endl;
  }
  return;
}

void PHActsTrkFitter::updateActsProtoTrack(const FitResult& fitOutput,
		  std::map<unsigned int, ActsTrack>::iterator iter)
{

  const auto& params = fitOutput.fittedParameters.value();
  
  /// We give the track params the vertex position because sometimes
  /// the fit returns a good momentum but bad position fit
  Acts::Vector4D position = Acts::Vector4D::Zero();
  position(0) = iter->second.getVertex()(0);
  position(1) = iter->second.getVertex()(1);
  position(2) = iter->second.getVertex()(2);
  position(3) = params.fourPosition(m_tGeometry->geoContext)(3);
  
  auto momVec = params.momentum();
  auto p = params.absoluteMomentum();
  auto q = params.charge();
  Acts::BoundSymMatrix cov;
  cov << 1000 * Acts::UnitConstants::um, 0., 0., 0., 0., 0.,
           0., 1000 * Acts::UnitConstants::um, 0., 0., 0., 0.,
           0., 0., 0.05, 0., 0., 0.,
           0., 0., 0., 0.05, 0., 0.,
           0., 0., 0., 0., 0.00005 , 0.,
           0., 0., 0., 0., 0., 1.;

  if(params.covariance())
    cov = params.covariance().value();

  if(Verbosity() > 1)
    std::cout << "Updating track seed with (" 
	      << momVec.x() << ", " << momVec.y() << ", " 
	      << momVec.z() << ")" << std::endl;

  const ActsExamples::TrackParameters siliconMMFit(position, momVec, 
						   p, q, cov);
  
  /// Update the track seed parameters
  iter->second.setTrackParams(siliconMMFit);
}

void PHActsTrkFitter::getTrackFitResult(const FitResult &fitOutput,
					const unsigned int trackKey,
					const Acts::Vector3D vertex)
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
  
  Trajectory trajectory(fitOutput.fittedStates, 
			trackTips, indexedParams);
  
  /// Get position, momentum from the Acts output. Update the values of
  /// the proto track
  
  auto updateTrackTimer = std::make_unique<PHTimer>("UpdateTrackTimer");
  updateTrackTimer->stop();
  updateTrackTimer->restart();
  if(fitOutput.fittedParameters)
    updateSvtxTrack(trajectory, trackKey, vertex);
  
  updateTrackTimer->stop();
  auto updateTime = updateTrackTimer->get_accumulated_time();
  
  if(Verbosity() > 0)
    std::cout << "PHActsTrkFitter update SvtxTrack time "
	      << updateTime << std::endl;

  if(m_timeAnalysis)
    h_updateTime->Fill(updateTime);
  
  /// Insert a new entry into the map
  m_actsFitResults->insert(
        std::pair<const unsigned int, Trajectory>(trackKey, 
						  trajectory));
  
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
	    }
	}
    
 

      rotater->calculateDCA(params, vertex, rotatedCov,
			    m_tGeometry->geoContext, 
			    dca3Dxy, dca3Dz, 
			    dca3DxyCov, dca3DzCov);
    }
 
  /// Set the DCA here. The DCA will be updated after the final
  /// vertex fitting in PHActsVertexFinder
  track->set_dca3d_xy(dca3Dxy / Acts::UnitConstants::cm);
  track->set_dca3d_z(dca3Dz / Acts::UnitConstants::cm);
  track->set_dca3d_xy_error(sqrt(dca3DxyCov) / Acts::UnitConstants::cm);
  track->set_dca3d_z_error(sqrt(dca3DzCov) / Acts::UnitConstants::cm);
  
  // Also need to update the state list and cluster ID list for all measurements associated with the acts track  
  // loop over acts track states, copy over to SvtxTrackStates, and add to SvtxTrack

  auto trackStateTimer = std::make_unique<PHTimer>("TrackStateTimer");
  trackStateTimer->stop();
  trackStateTimer->restart();
  
  if(m_fillSvtxTrackStates)
    rotater->fillSvtxTrackStates(traj, trackTip, track,
				 m_tGeometry->geoContext,
				 m_hitIdClusKey);  

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

  m_hitIdClusKey = findNode::getClass<CluskeyBimap>(topNode, "HitIDClusIDActsMap");
  
  if (!m_hitIdClusKey)
    {
      std::cout << PHWHERE << "No HitID:ClusKey map on node tree. Bailing."
		<< std::endl;
      
      return Fun4AllReturnCodes::EVENT_OK;
    }


  return Fun4AllReturnCodes::EVENT_OK;
}

