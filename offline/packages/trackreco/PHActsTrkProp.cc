
#include "PHActsTrkProp.h"
#include "MakeActsGeometry.h"
#include "ActsTrack.h"
#include "ActsTransformations.h"

/// Fun4All includes
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <phool/phool.h>

/// Tracking includes
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrack_v1.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackState_v1.h>

#include <Acts/EventData/SingleCurvilinearTrackParameters.hpp>
#include <Acts/EventData/TrackParameters.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>

#include <Acts/Propagator/Propagator.hpp>
#include <Acts/Propagator/Navigator.hpp>
#include <Acts/Propagator/AbortList.hpp>
#include <Acts/Propagator/ActionList.hpp>
#include <Acts/Utilities/Helpers.hpp>
#include <Acts/Utilities/Units.hpp>
#include <Acts/TrackFinding/CombinatorialKalmanFilter.hpp>
#include <Acts/EventData/MultiTrajectoryHelpers.hpp>

#include <ActsExamples/Plugins/BField/BFieldOptions.hpp>
#include <ActsExamples/Plugins/BField/ScalableBField.hpp>
#include <ActsExamples/Framework/ProcessCode.hpp>
#include <ActsExamples/Framework/WhiteBoard.hpp>
#include <ActsExamples/EventData/Track.hpp>
#include <ActsExamples/Framework/AlgorithmContext.hpp>

#include <TFile.h>
#include <TH1.h>
#include <TMatrixDSym.h>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <vector>
#include <utility>
#include <chrono>

using namespace std::chrono;

PHActsTrkProp::PHActsTrkProp(const std::string& name)
  : PHTrackPropagating(name)
  , m_event(0)
  , m_timeAnalysis(false)
  , m_timeFile(nullptr)
  , h_eventTime(nullptr)
  , m_nBadFits(0)
  , m_tGeometry(nullptr)
  , m_resetCovariance(false)
  , m_trackMap(nullptr)
  , m_actsProtoTracks(nullptr)
  , m_actsFitResults(nullptr)
  , m_actsTrackKeyMap(nullptr)
  , m_hitIdClusKey(nullptr)
  , m_sourceLinks(nullptr)
  , m_topNode(nullptr)
{
  Verbosity(0);
  initializeLayerSelector();
}

 PHActsTrkProp::~PHActsTrkProp()
{
}

int PHActsTrkProp::Setup(PHCompositeNode* topNode)
{
  if(Verbosity() > 0)
    std::cout << "Setup PHActsTrkProp" << std::endl;
  createNodes(topNode);
  
  if (getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
    return Fun4AllReturnCodes::ABORTEVENT;

  /// Setup the source link selection criteria for the various layers
  setupSourceLinkSelection();
  
  /// Setup the Acts general/generic propagator options
  /// see struct in acts/Core/include/Acts/Propagator/Propagator.hpp
  /// for additional options that can be changed
  m_actsPropPlainOptions = Acts::PropagatorPlainOptions();
  
  findCfg.finder = ActsExamples::TrkrClusterFindingAlgorithm::makeFinderFunction(
			         m_tGeometry->tGeometry,
				 m_tGeometry->magField);

  if(m_timeAnalysis)
    {
      m_timeFile = new TFile("PHActsTrkPropTime.root","recreate");
      h_eventTime = new TH1F("h_eventTime",";time [ms]",10000,0,10000);
    }


  if(Verbosity() > 0)
    std::cout <<" Finish PHActsTrkProp setup" << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

std::vector<SourceLink> PHActsTrkProp::getEventSourceLinks()
{
  std::vector<SourceLink> sourceLinks;
  std::map<unsigned int, SourceLink>::iterator slIter = m_sourceLinks->begin();
  while(slIter != m_sourceLinks->end())
    {
      sourceLinks.push_back(slIter->second);
      
      if(Verbosity() > 10)
	{
	  std::cout << std::endl 
		    << "Adding source link to list for track finding: " 
		    << slIter->second.hitID() <<" and surface : " 
		    << std::endl;
	  slIter->second.referenceSurface().toStream(m_tGeometry->geoContext, std::cout);
	  
	}
      ++slIter;
    }

  return sourceLinks;

}

int PHActsTrkProp::Process()
{
  auto startTime = high_resolution_clock::now();
  
  m_event++;


  auto logLevel = Acts::Logging::INFO;

  if (Verbosity() > 0)
  {
    std::cout << PHWHERE << "Events processed: " << m_event << std::endl;
    std::cout << "Start PHActsTrkProp::process_event" << std::endl;
    logLevel = Acts::Logging::VERBOSE;
  }
  
  auto logger = Acts::getDefaultLogger("PHActsTrkProp", logLevel);
    
  /// Collect all source links for the CKF
  std::vector<SourceLink> sourceLinks = getEventSourceLinks();

  int nTrackSeeds = 0;
  
  std::map<unsigned int, ActsTrack>::iterator trackIter;
  for(trackIter = m_actsProtoTracks->begin();
      trackIter != m_actsProtoTracks->end();
      ++trackIter)
  {
    ActsTrack track = trackIter->second;
    const unsigned int trackKey = trackIter->first;

    ActsExamples::TrackParameters trackSeed = track.getTrackParams();

    if(m_resetCovariance)
      {
	if(Verbosity() > 0)
	  std::cout << "PHActsTrkProp : resetting covariance"<<std::endl;
	Acts::BoundSymMatrix covariance;

	/// If we are resetting the covariance and the space point
	/// to the vertex, the POCA should have the covariance of the
	/// vertex
	covariance << 50 * Acts::UnitConstants::um, 0., 0., 0., 0., 0.,
	              0., 25 * Acts::UnitConstants::um, 0., 0., 0., 0.,
	              0., 0., 0.01, 0., 0., 0.,
	              0., 0., 0., 0.01, 0., 0.,
	              0., 0., 0., 0., 0.0001, 0.,
	              0., 0., 0., 0., 0., 1.;

	Acts::Vector4D new4Vec(track.getVertex().x(),
			       track.getVertex().y(),
			       track.getVertex().z(),
			       trackSeed.time());

	ActsExamples::TrackParameters trackSeedNewCov(
				      new4Vec,
				      trackSeed.momentum(),
				      trackSeed.absoluteMomentum(),
				      trackSeed.charge(),
				      covariance);

	trackSeed = trackSeedNewCov;
      }

    /// Construct a perigee surface as the target surface
    /// This surface is what Acts fits with respect to, so we set it to
    /// the initial vertex estimation
    auto pSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(
	                           track.getVertex());

    if(Verbosity() > 0)
      {
	std::cout << "Processing track seed with positon: "
		  << trackSeed.position(m_tGeometry->geoContext).transpose() 
		  << std::endl
		  << "momentum: " << trackSeed.momentum().transpose() 
		  << std::endl
		  << "charge: " << trackSeed.charge() << std::endl
		  << "initial vertex : " <<track.getVertex()
		  << " corresponding to SvtxTrack key " << trackKey
		  << std::endl
		  << "with initial covariance " << std::endl
		  << trackSeed.covariance().value() << std::endl;
      }


    /// Construct the options to pass to the CKF.
    /// SourceLinkSelector set in Init()
    Acts::CombinatorialKalmanFilterOptions<SourceLinkSelector> ckfOptions(
	        m_tGeometry->geoContext, 
		m_tGeometry->magFieldContext, 
		m_tGeometry->calibContext, 
		m_sourceLinkSelectorConfig,
		Acts::LoggerWrapper(*logger),
		m_actsPropPlainOptions,
		&(*pSurface));

    /// Run the CKF for all source links and the constructed track seed
    /// CKF runs both track finder and KF track fitter
    auto result = findCfg.finder(sourceLinks, trackSeed, ckfOptions);
    
    if(result.ok())
      {
	const CKFFitResult& fitOutput = result.value();

	Trajectory traj(fitOutput.fittedStates,
			fitOutput.trackTips,
			fitOutput.fittedParameters);

	m_actsFitResults->insert(std::pair<const unsigned int, Trajectory>
				 (trackKey, traj));
	
	updateSvtxTrack(traj, trackKey, track.getVertex());
      }
    else
      {
	m_actsFitResults->insert(std::pair<const unsigned int, Trajectory>
				 (trackKey, ActsExamples::TrkrClusterMultiTrajectory()));
		
	m_nBadFits++;
      }

    nTrackSeeds++;
  }
  

  if(Verbosity() > 0)
    {
      std::cout << "Finished process_event for PHActsTrkProp" 
		<< std::endl;
      std::cout << "Processed " << nTrackSeeds << " track seeds "
		<< std::endl;
    }

  auto stopTime = high_resolution_clock::now();
  auto eventTime = duration_cast<microseconds>(stopTime - startTime);

  if(m_timeAnalysis)
    h_eventTime->Fill(eventTime.count() / 1000.);

  return Fun4AllReturnCodes::EVENT_OK;
}
 
int PHActsTrkProp::ResetEvent(PHCompositeNode *topNode)
{
  m_actsFitResults->clear();
  m_actsTrackKeyMap->clear();

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsTrkProp::End()
{

  if(m_timeAnalysis)
    {
      m_timeFile->cd();
      h_eventTime->Write();
      m_timeFile->Write();
      m_timeFile->Close();
    }

  if (Verbosity() > 0)
  {
    std::cout << "Finished PHActsTrkProp" << std::endl;
    std::cout << "PHActsTrkProp had " << m_nBadFits 
	      << " bad fits" << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

void PHActsTrkProp::updateSvtxTrack(Trajectory traj, 
				    const unsigned int trackKey,
				    Acts::Vector3D vertex)
{
  
  ActsTransformations *rotater = new ActsTransformations();
  rotater->setVerbosity(0);

  int iTrack = 0;

  /// There could be multiple tracks per trajectory found from the CKF
  /// Therefore we need to get the info from the original SvtxTrack that
  /// we need, and then create new tracks
  SvtxTrack *origTrack = m_trackMap->find(trackKey)->second;
  
  /// All these tracks should have the same vertexId if they came from
  /// the same track seed
  unsigned int vertexId = origTrack->get_vertex_id();
  
  /// hack for now to just set the vertex id to 0 if it is needed
  if(vertexId == UINT_MAX)
    vertexId = 0;

  origTrack->Reset();


  /// Create a state at 0.0 to propagate out from
  float pathLength = 0.0;
  SvtxTrackState_v1 out(pathLength);
  out.set_x(0.0);
  out.set_y(0.0);
  out.set_z(0.0);
  origTrack->insert_state(&out);

  /// This gets the track indexer and associated tracks 
  /// (Acts::MultiTrajectories)
  const auto& [trackTips, mj] = traj.trajectory();
  
  if(trackTips.empty()) 
    {
      if(Verbosity() > 5)
	std::cout << "Empty multiTrajectory..." << std::endl;
      return;
    }
  
  std::map<const size_t, const unsigned int> trackTipTrackKeyMap;

  /// Iterate over the found tracks via their trackTip index
  for(const size_t& trackTip : trackTips) 
    {
      auto trajState = Acts::MultiTrajectoryHelpers::trajectoryState(mj, trackTip);
      
      size_t nStates = trajState.nStates;
      size_t nMeasurements = trajState.nMeasurements;
      size_t nOutliers = trajState.nOutliers;
      size_t nHoles = trajState.nHoles;
      double chi2sum = trajState.chi2Sum;
      size_t NDF = trajState.NDF;

      /// No trackParameters for this trackTip, so fit failed for this tip
      if( !traj.hasTrackParameters(trackTip) )
	continue;
      
      const auto& fittedParameters = traj.trackParameters(trackTip);
    
      float x  = fittedParameters.position(m_tGeometry->geoContext)(0) 
	/ Acts::UnitConstants::cm;
      float y  = fittedParameters.position(m_tGeometry->geoContext)(1) 
	/ Acts::UnitConstants::cm;
      float z  = fittedParameters.position(m_tGeometry->geoContext)(2) 
	/ Acts::UnitConstants::cm;
      float px = fittedParameters.momentum()(0);
      float py = fittedParameters.momentum()(1);
      float pz = fittedParameters.momentum()(2);
      
      if(Verbosity() > 0)
	{
	  std::cout << "Track fit returned a track with: " << std::endl
		    << "momentum : " 
		    << fittedParameters.momentum().transpose()<< std::endl
		    << "position : " 
		    << fittedParameters.position(m_tGeometry->geoContext).transpose()
		    << std::endl
		    << "charge : " << fittedParameters.charge() << std::endl;
	  std::cout << "Track has " << nStates << " states and " 
		    << nMeasurements << " measurements and " 
		    << nHoles << " holes and " << nOutliers << " outliers "
		    << std::endl;
	}
      
      float qOp = fittedParameters.parameters()[Acts::eBoundQOverP];
      
      Acts::BoundSymMatrix rotatedCov = Acts::BoundSymMatrix::Zero();
      if(fittedParameters.covariance())
	{
	  rotater->setVerbosity(0);
	  rotatedCov = rotater->rotateActsCovToSvtxTrack(fittedParameters, 
							 m_tGeometry->geoContext);
	}
      
      float DCA3Dxy = -9999;
      float DCA3Dz = -9999;
      float DCA3DxyCov = -9999;
      float DCA3DzCov = -9999;
      
      rotater->calculateDCA(fittedParameters, vertex, m_tGeometry->geoContext,
			    DCA3Dxy, DCA3Dz, DCA3DxyCov, DCA3DzCov);
      
      /// If it is the first track, just update the original track seed
      if(iTrack == 0)
	{
	  if(Verbosity() > 0)
	    std::cout << "Replacing SvtxTrack " << trackKey 
		      << " with parameters from trackTip " << trackTip
		      << std::endl;
	  
	  origTrack->set_id(trackKey);
	  origTrack->set_vertex_id(vertexId);
	  origTrack->set_chisq(chi2sum);
	  origTrack->set_ndf(NDF);
	  origTrack->set_charge(1);

	  if(qOp < 0)
	    origTrack->set_charge(-1);
	  for(int i = 0; i < 6; ++i)
	    {
	      for(int j = 0; j < 6; ++j)
		{
		  origTrack->set_error(i,j, rotatedCov(i,j));
		}
	    }

	  origTrack->set_px(px);
	  origTrack->set_py(py);
	  origTrack->set_pz(pz);
	  origTrack->set_x(x);
	  origTrack->set_y(y);
	  origTrack->set_z(z);
	  origTrack->set_dca3d_xy(DCA3Dxy / Acts::UnitConstants::cm);
	  origTrack->set_dca3d_z(DCA3Dz / Acts::UnitConstants::cm);
	  origTrack->set_dca3d_xy_error(DCA3DxyCov / Acts::UnitConstants::cm);
	  origTrack->set_dca3d_z_error(DCA3DzCov / Acts::UnitConstants::cm);
	  rotater->setVerbosity(0);
	  rotater->fillSvtxTrackStates(traj, trackTip, origTrack,
				       m_tGeometry->geoContext,
				       m_hitIdClusKey);
        
	  trackTipTrackKeyMap.insert(std::pair<const size_t, const unsigned int>
				      (trackTip, trackKey));

	
	}
      else
	{
	  /// Get the last track key so that we can add new tracks if needed
	  
	  const unsigned int lastTrackKey = m_trackMap->end()->first;

	  /// Otherwise make a new track
	  if(Verbosity() > 0)
	    std::cout << "Creating new SvtxTrack with trackKey " << lastTrackKey
		      << " corresponding to trackTip " << trackTip << std::endl;

	
	  SvtxTrack_v1 *newTrack = new SvtxTrack_v1();
	  /// Needs to be a new track id
	  newTrack->set_id(lastTrackKey);
	  newTrack->set_vertex_id(vertexId);
	  newTrack->set_chisq(chi2sum);
	  newTrack->set_ndf(NDF);
	  newTrack->set_charge(1);
	  /// If negative QOP then set charge to negative
	  if(qOp < 0)
	    newTrack->set_charge(-1);
	  
	  newTrack->set_px(px);
	  newTrack->set_py(py);
	  newTrack->set_pz(pz);
	  newTrack->set_x(x);
	  newTrack->set_y(y);
	  newTrack->set_z(z);
	  newTrack->set_dca3d_xy(DCA3Dxy / Acts::UnitConstants::cm);
	  newTrack->set_dca3d_z(DCA3Dz / Acts::UnitConstants::cm);
	  newTrack->set_dca3d_xy_error(DCA3DxyCov / Acts::UnitConstants::cm);
	  newTrack->set_dca3d_z_error(DCA3DzCov / Acts::UnitConstants::cm);
	  
	  rotater->setVerbosity(0);
	  rotater->fillSvtxTrackStates(traj, trackTip, 
				       dynamic_cast<SvtxTrack*>(newTrack),
				       m_tGeometry->geoContext,
				       m_hitIdClusKey);
        
	  for(int i = 0; i < 6; ++i)
	    {
	      for(int j = 0; j < 6; ++j)
		{
		  newTrack->set_error(i,j, rotatedCov(i,j));
		}
	    }
	  

	  trackTipTrackKeyMap.insert(std::pair<const size_t, const unsigned int>
				      (trackTip, lastTrackKey));

	  m_trackMap->insert(newTrack);
	}
      
      ++iTrack;
    }

  /// Now insert the trajectory + trackTip-trackKey pair
  m_actsTrackKeyMap->insert(std::pair<const unsigned int, 
			    std::map<const size_t, const unsigned int>>
			    (trackKey, trackTipTrackKeyMap));


}


void PHActsTrkProp::setupSourceLinkSelection()
{
  
  /// The CKF requires a source link selector which helps it identify possible
  /// SLs to the track. The selections can be added like
  /// {makeId(volId, layerId), {maxChi2, numSourceLinks}}
  /// We'll just put the max chi2 to 10 

  std::vector<std::pair<Acts::GeometryIdentifier,
			Acts::SourceLinkSelectorCuts>> sourceLinkSelectors;
  
  /// Global detector criteria
  sourceLinkSelectors.push_back({makeId(), {100.,53}});

  /// Volume criteria
  /// Volume IDs - MVTX = 7, INTT = 9, TPC = 11
  sourceLinkSelectors.push_back({makeId(7), {m_volMaxChi2.find(7)->second, 3}});
  sourceLinkSelectors.push_back({makeId(9), {m_volMaxChi2.find(9)->second, 2}});
  sourceLinkSelectors.push_back({makeId(11), {m_volMaxChi2.find(11)->second, 48}});
  
  /// Set individual layer criteria (e.g. for first layer of TPC)
  /// Individual layers should only have one SL per layer
  for(int vol = 0; vol < m_volLayerMaxChi2.size(); ++vol)
    for(std::pair<const int, const float> element : m_volLayerMaxChi2.at(vol))
      ///Vol IDs are 7, 9, 11, hence vol*2+7
      sourceLinkSelectors.push_back({makeId(vol * 2. + 7, element.first),
	                            {element.second, 1}});

  m_sourceLinkSelectorConfig = SourceLinkSelectorConfig(sourceLinkSelectors);

  if(Verbosity() > 1)
    {
      std::cout << "The source link selection criteria were set to: " << std::endl;
      for(int i = 0; i < sourceLinkSelectors.size(); i++)
	{
	  std::cout << "GeoID : " << sourceLinkSelectors.at(i).first
		    << " has chi sq selection " 
		    << sourceLinkSelectors.at(i).second.chi2CutOff
		    << std::endl;
	}
    }
}

Acts::GeometryIdentifier PHActsTrkProp::makeId(int volume, 
				       int layer, 
				       int sensitive)
{
  return Acts::GeometryIdentifier().setVolume(volume)
                                   .setLayer(layer)
                                   .setSensitive(sensitive);
}

void PHActsTrkProp::setVolumeLayerMaxChi2(const int vol, 
					  const int layer,
					  const float maxChi2)
{
  int volume;
  if(vol == 7)
    volume = 0;
  else if(vol == 9)
    volume = 1;
  else if(vol == 11)
    volume = 2;
  else
    {
      std::cout << "Invalid volume number supplied. Not including" << std::endl;
      return;
    }

  m_volLayerMaxChi2.at(volume).insert(std::make_pair(layer, maxChi2));

}
void PHActsTrkProp::setVolumeMaxChi2(const int vol, const float maxChi2)
{
  m_volMaxChi2.insert(std::make_pair(vol, maxChi2));
}

void PHActsTrkProp::createNodes(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);

  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));

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

  m_actsTrackKeyMap = findNode::getClass<std::map<const unsigned int, 
				std::map<const size_t, const unsigned int>>>
                                (topNode, "ActsTrackKeys");
  if(!m_actsTrackKeyMap)
    {
      m_actsTrackKeyMap = new std::map<const unsigned int, 
				       std::map<const size_t, const unsigned int>>;
      PHDataNode<std::map<const unsigned int, 
			  std::map<const size_t, const unsigned int>>> *fitNode = 
	new PHDataNode<std::map<const unsigned int, 
				std::map<const size_t, const unsigned int>>>
	(m_actsTrackKeyMap, "ActsTrackKeys");

      svtxNode->addNode(fitNode);
    }

  m_actsFitResults = findNode::getClass<std::map<const unsigned int, Trajectory>>(topNode, "ActsFitResults");
  if(!m_actsFitResults)
    {
      m_actsFitResults = new std::map<const unsigned int, Trajectory>;
      PHDataNode<std::map<const unsigned int, Trajectory>> *fitNode = 
	new PHDataNode<std::map<const unsigned int, Trajectory>>(m_actsFitResults, "ActsFitResults");
      
      svtxNode->addNode(fitNode);

    }


  return;
}

int PHActsTrkProp::getNodes(PHCompositeNode* topNode)
{

  m_topNode = topNode;

  m_actsProtoTracks = findNode::getClass<std::map<unsigned int, ActsTrack>>(topNode, "ActsTrackMap");

  if(!m_actsProtoTracks)
    {
      std::cout << PHWHERE << "Acts proto tracks not on node tree. Bailing." 
		<< std::endl;

      return Fun4AllReturnCodes::ABORTEVENT;
    }

  m_sourceLinks = findNode::getClass<std::map<unsigned int, SourceLink>>(topNode, "TrkrClusterSourceLinks");

  if (!m_sourceLinks)
    {
      std::cout << PHWHERE << "TrkrClusterSourceLinks node not found on node tree. Exiting."
		<< std::endl;
      
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  
  m_tGeometry = findNode::getClass<ActsTrackingGeometry>(topNode, "ActsTrackingGeometry");

  if (!m_tGeometry)
  {
    std::cout << "ActsGeometry not on node tree. Exiting."
              << std::endl;

    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_trackMap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if (!m_trackMap)
    {
      std::cout << PHWHERE << " ERROR: Can't find SvtxTrackMap. Exiting " 
		<< std::endl;

      return Fun4AllReturnCodes::ABORTEVENT;
    }

  m_hitIdClusKey = findNode::getClass<CluskeyBimap>(topNode, "HitIDClusIDActsMap");
  if(!m_hitIdClusKey)
    {
      std::cout << PHWHERE << "ERROR: Can't find HitIdClusIdActsMap. Exiting."
		<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;

    }
  

  return Fun4AllReturnCodes::EVENT_OK;
}

void PHActsTrkProp::initializeLayerSelector()
{
  std::map<const int, const float> dumMap;
  m_volLayerMaxChi2.push_back(dumMap);
  m_volLayerMaxChi2.push_back(dumMap);
  m_volLayerMaxChi2.push_back(dumMap);
}
