
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

#include <Acts/EventData/ChargePolicy.hpp>
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
#include <Acts/TrackFinder/CombinatorialKalmanFilter.hpp>
#include <Acts/EventData/MultiTrajectoryHelpers.hpp>

#include <ACTFW/Plugins/BField/BFieldOptions.hpp>
#include <ACTFW/Plugins/BField/ScalableBField.hpp>
#include <ACTFW/Framework/ProcessCode.hpp>
#include <ACTFW/Framework/WhiteBoard.hpp>
#include <ACTFW/EventData/Track.hpp>
#include <ACTFW/Framework/AlgorithmContext.hpp>

#include <TMatrixDSym.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <utility>

PHActsTrkProp::PHActsTrkProp(const std::string& name)
  : PHTrackPropagating(name)
  , m_event(0)
  , m_nBadFits(0)
  , m_tGeometry(nullptr)
  , m_trackMap(nullptr)
  , m_actsProtoTracks(nullptr)
  , m_actsFitResults(nullptr)
  , m_actsTrackKeyMap(nullptr)
  , m_hitIdClusKey(nullptr)
  , m_sourceLinks(nullptr)
  , m_topNode(nullptr)
{
  Verbosity(0);
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

  /// The CKF requires a source link selector which helps it identify possible
  /// SLs to the track. The selections can be added like
  /// {makeId(volId, layerId), {maxChi2, numSourceLinks}}
  /// We'll just put the max chi2 to 10 
  m_sourceLinkSelectorConfig = {
    /// global default values
    {makeId(), {10.0, 55}},

    /// MVTX volume should have max 3 SLs
    {makeId(7), {10.0, 3}},
    
    /// INTT volume should have max 2 SLs
    {makeId(9), {10.0, 2}},
    
    /// TPC volume should have max 50 (?) SLs
    {makeId(11), {10.0,50}}
  };

  auto logger = Acts::Logging::INFO;
  if(Verbosity() > 5)
    logger = Acts::Logging::VERBOSE;

  findCfg.finder = FW::TrkrClusterFindingAlgorithm::makeFinderFunction(
                   m_tGeometry->tGeometry,
		   m_tGeometry->magField,
		   logger);

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

  m_event++;

  if (Verbosity() > 0)
  {
    std::cout << PHWHERE << "Events processed: " << m_event << std::endl;
    std::cout << "Start PHActsTrkProp::process_event" << std::endl;
  }

  /// Collect all source links for the CKF
  std::vector<SourceLink> sourceLinks = getEventSourceLinks();

  std::map<unsigned int, ActsTrack>::iterator trackIter;
  for(trackIter = m_actsProtoTracks->begin();
      trackIter != m_actsProtoTracks->end();
      ++trackIter)
  {
    ActsTrack track = trackIter->second;
    const unsigned int trackKey = trackIter->first;

    FW::TrackParameters trackSeed = track.getTrackParams();

    Acts::BoundSymMatrix covariance;
    covariance << 1000 * Acts::UnitConstants::um, 0., 0., 0., 0., 0.,
                  0., 1000 * Acts::UnitConstants::um, 0., 0., 0., 0.,
                  0., 0., 0.05, 0., 0., 0.,
                  0., 0., 0., 0.05, 0., 0.,
                  0., 0., 0., 0., 0.001, 0.,
                  0., 0., 0., 0., 0., 1.;
    
    FW::TrackParameters trackSeedNewCov(covariance,
					trackSeed.position(),
					trackSeed.momentum(),
					trackSeed.charge(),
					trackSeed.time());

    /// Construct a perigee surface as the target surface
    /// This surface is what Acts fits with respect to, so we set it to
    /// the initial vertex estimation
    auto pSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(
	                           track.getVertex());

    if(Verbosity() > 0)
      {
	std::cout << "Processing track seed with positon: "
		  << trackSeed.position().transpose() << std::endl
		  << "momentum: " << trackSeed.momentum().transpose() 
		  << std::endl
		  << "charge: " << trackSeed.charge() << std::endl
		  << "initial vertex : " <<track.getVertex()
		  << " corresponding to SvtxTrack key " << trackKey
		  << std::endl
		  << "with initial covariance " << std::endl
		  << covariance << std::endl;
      }


    /// Construct the options to pass to the CKF.
    /// SourceLinkSelector set in Init()
    Acts::CombinatorialKalmanFilterOptions<SourceLinkSelector> ckfOptions(
	        m_tGeometry->geoContext, 
		m_tGeometry->magFieldContext, 
		m_tGeometry->calibContext, 
		m_sourceLinkSelectorConfig, 
		&(*pSurface));

    /// Run the CKF for all source links and the constructed track seed
    /// CKF runs both track finder and KF track fitter
    auto result = findCfg.finder(sourceLinks, trackSeedNewCov, ckfOptions);
    
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
				 (trackKey, FW::TrkrClusterMultiTrajectory()));
		
	m_nBadFits++;
      }
    
  }
  

  if(Verbosity() > 0)
    std::cout << "Finished process_event for PHActsTrkProp" << std::endl;


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
  const unsigned int vertexId = origTrack->get_vertex_id();
  
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
      
      /// No trackParameters for this trackTip, so fit failed for this tip
      if( !traj.hasTrackParameters(trackTip) )
	continue;
      
      const auto& fittedParameters = traj.trackParameters(trackTip);
    
      float x  = fittedParameters.position()(0) / Acts::UnitConstants::cm;
      float y  = fittedParameters.position()(1) / Acts::UnitConstants::cm;
      float z  = fittedParameters.position()(2) / Acts::UnitConstants::cm;
      float px = fittedParameters.momentum()(0);
      float py = fittedParameters.momentum()(1);
      float pz = fittedParameters.momentum()(2);
      
      if(Verbosity() > 0)
	{
	  std::cout << "Track fit returned a track with: " << std::endl
		    << "momentum : " << fittedParameters.momentum().transpose()
		    << std::endl
		    << "position : " << fittedParameters.position().transpose()
		    << std::endl
		    << "charge : " << fittedParameters.charge() << std::endl;
	}
      
      float qOp = fittedParameters.parameters()[Acts::ParDef::eQOP];
      
      Acts::BoundSymMatrix rotatedCov = Acts::BoundSymMatrix::Zero();
      if(fittedParameters.covariance())
	{
	  rotater->setVerbosity(0);
	  rotatedCov = rotater->rotateActsCovToSvtxTrack(fittedParameters);
	}
      
      
      //size_t nStates = trajState.nStates;
      //size_t nMeasurements = trajState.nMeasurements;
      double chi2sum = trajState.chi2Sum;
      size_t NDF = trajState.NDF;
      
      float DCA3Dxy = -9999;
      float DCA3Dz = -9999;
      float DCA3DxyCov = -9999;
      float DCA3DzCov = -9999;
      
      rotater->calculateDCA(fittedParameters, vertex,
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

Acts::GeometryID PHActsTrkProp::makeId(int volume, 
				       int layer, 
				       int sensitive)
{
  return Acts::GeometryID().setVolume(volume)
                           .setLayer(layer)
                           .setSensitive(sensitive);
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

  m_hitIdClusKey = findNode::getClass<std::map<TrkrDefs::cluskey, unsigned int>>(topNode, "HitIDClusIDActsMap");
  if(!m_hitIdClusKey)
    {
      std::cout << PHWHERE << "ERROR: Can't find HitIdClusIdActsMap. Exiting."
		<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;

    }
  

  return Fun4AllReturnCodes::EVENT_OK;
}

