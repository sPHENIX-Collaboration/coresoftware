
#include "PHActsTrkProp.h"
#include "MakeActsGeometry.h"
#include "ActsTrack.h"
#include "ActsCovarianceRotater.h"

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
  , m_actsFitResults(nullptr)
  , m_hitIdClusKey(nullptr)
  , m_sourceLinks(nullptr)
{
  Verbosity(0);
}

 PHActsTrkProp::~PHActsTrkProp()
{
}

int PHActsTrkProp::Setup(PHCompositeNode* topNode)
{
  if(Verbosity() > 1)
    std::cout << "Setup PHActsTrkProp" << std::endl;
  createNodes(topNode);
  
  if (getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
    return Fun4AllReturnCodes::ABORTEVENT;

  /// Need to implement different chi2 criteria for the different layers
  /// and volumes to help the CKF out a little 
  /// First need to walk through layers/volumes to figure out appropriate
  /// numerical identifiers
  /// Leave these out for now once we return to track propagation
  /// m_sourceLinkSelectorConfig.layerMaxChi2 = {{2, {{2, 8}, {4, 7}}}};
  /// m_sourceLinkSelectorConfig.volumeMaxChi2 = {{2, 7}, {3, 8}};
  /// m_sourceLinkSelectorConfig.maxChi2 = 8;
  /// Set the allowed maximum number of source links to be large enough
  /// m_sourceLinkSelectorConfig.maxNumSourcelinksOnSurface = 100;
 
  findCfg.finder = FW::TrkrClusterFindingAlgorithm::makeFinderFunction(
                   m_tGeometry->tGeometry,
		   m_tGeometry->magField,
		   Acts::Logging::INFO);

  if(Verbosity() > 1)
    std::cout <<" Finish PHActsTrkProp setup" << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsTrkProp::Process()
{

  m_event++;

  if (Verbosity() > 1)
  {
    std::cout << PHWHERE << "Events processed: " << m_event << std::endl;
    std::cout << "Start PHActsTrkProp::process_event" << std::endl;
  }


  PerigeeSurface pSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>
    (Acts::Vector3D(0., 0., 0.));

  /// Collect all source links for the CKF
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

  ActsCovarianceRotater *rotater = new ActsCovarianceRotater();
  rotater->setVerbosity(Verbosity());

  for (SvtxTrackMap::Iter trackIter = m_trackMap->begin();
       trackIter != m_trackMap->end(); ++trackIter)
    {
      const unsigned int trackKey = trackIter->first;
      SvtxTrack *track = trackIter->second;
      
      if (!track)
	continue;
      
      if (Verbosity() > 1)
	{
	  std::cout << "found SvtxTrack " << trackIter->first << std::endl;
	  track->identify();
	}
      
      /// Get the necessary parameters and values for the TrackParameters
      const Acts::BoundSymMatrix seedCov = rotater->rotateSvtxTrackCovToActs(track);
      const Acts::Vector3D seedPos(track->get_x(),
				   track->get_y(),
				   track->get_z());
      const Acts::Vector3D seedMom(track->get_px(),
				   track->get_py(),
				   track->get_pz());
      
      // just set to 10 ns for now
      const double trackTime = 10 * Acts::UnitConstants::ns;
      const int trackQ = -1*track->get_charge();
      
      const FW::TrackParameters trackSeed(seedCov, seedPos,
					  seedMom, trackQ, trackTime);
      
      
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
      auto result = findCfg.finder(sourceLinks, trackSeed, ckfOptions);
      
      if(result.ok())
	{
	  const CKFFitResult& fitOutput = result.value();
	  auto parameterMap = fitOutput.fittedParameters;

	  Trajectory traj(fitOutput.fittedStates,
			  fitOutput.trackTips,
			  fitOutput.fittedParameters);
	
	  m_actsFitResults->insert(std::pair<const unsigned int, Trajectory>
				   (trackKey, traj));

	}
      else
	{
	  m_actsFitResults->insert(std::pair<const unsigned int, Trajectory>
				   (trackKey, FW::TrkrClusterMultiTrajectory()));

	  /// Track prop failed - set to arbitrary values
	  track->set_x(-9999.);
	  track->set_y(-9999.);
	  track->set_z(-9999.);
	  track->set_px(-9999.);
	  track->set_py(-9999.);
	  track->set_pz(-9999.);
	  
	  m_nBadFits++;
	}
	
    }

  /// Update the SvtxTrackMap by wiping it clean and adding the tracks from the CKF
  updateSvtxTrackMap();

  if(Verbosity() > 1)
    std::cout << "Finished process_event for PHActsTrkProp" << std::endl;


  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsTrkProp::ResetEvent(PHCompositeNode *topNode)
{
  m_actsFitResults->clear();

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsTrkProp::End()
{
  if (Verbosity() > 10)
  {
    std::cout << "Finished PHActsTrkProp" << std::endl;
  }

  std::cout << "PHActsTrkProp had " << m_nBadFits << " bad fits" << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

void PHActsTrkProp::updateSvtxTrackMap()
{
  
  /// Wipe the track map completely, since the CKF can return multiple tracks for a given
  /// track seed the trackseed track keys no longer have meaning
  m_trackMap->Reset();

  ActsCovarianceRotater *rotater = new ActsCovarianceRotater();
  rotater->setVerbosity(Verbosity());

  /// Iterate over the Trajectories and add them to the SvtxTrackMap
  std::map<const unsigned int, Trajectory>::iterator trackIter;
  int iTrack = 0;

  for (trackIter = m_actsFitResults->begin();
       trackIter != m_actsFitResults->end();
       ++trackIter)
  {
    Trajectory traj = trackIter->second;
 
    const auto& [trackTips, mj] = traj.trajectory();
    if(trackTips.empty()) 
      {
	if(Verbosity() > 5)
	  std::cout << "Empty multiTrajectory... continuing" << std::endl;
	continue;
      }
    
    for(const size_t& trackTip : trackTips) 
      {
	SvtxTrack_v1 track;
	track.set_id(iTrack);
    
	auto trajState = Acts::MultiTrajectoryHelpers::trajectoryState(mj, trackTip);
	
	if(! traj.hasTrackParameters(trackTip))
	  continue;
	
	const auto& fittedParameters = traj.trackParameters(trackTip);
	
	track.set_x(fittedParameters.position()(0) / Acts::UnitConstants::cm);
	track.set_y(fittedParameters.position()(1) / Acts::UnitConstants::cm);
	track.set_z(fittedParameters.position()(2) / Acts::UnitConstants::cm);
	track.set_px(fittedParameters.momentum()(0));
	track.set_py(fittedParameters.momentum()(1));
	track.set_pz(fittedParameters.momentum()(2));
	if(fittedParameters.covariance())
	  {
	    Acts::BoundSymMatrix rotatedCov = rotater->rotateActsCovToSvtxTrack(fittedParameters);
	    for(int i = 0; i < 6; ++i)
	      {
		for(int j = 0; j < 6; ++j)
		  {
		    track.set_error(i,j, rotatedCov(i,j));
		  }
	      }
	  }

	/// Loop over trajectory source links somehow ? TBD
	/// Logic will be: 
	/// SL->hitID()
	/// hitID -> clusterkey via map created in PHActsSourceLinks
	/// add cluster key to SvtxTrack
	/// Will want to match with clusterkey as determined in PHActsSourceLinks
	
	size_t nStates = trajState.nStates;
	size_t nMeasurements = trajState.nMeasurements;
	double chi2sum = trajState.chi2Sum;
	size_t NDF = trajState.NDF;
	
	++iTrack;
	m_trackMap->insert(&track);
      }

  }


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

  m_actsFitResults = findNode::getClass<std::map<const unsigned int, Trajectory>>(topNode, "ActsCKFResults");
  if(!m_actsFitResults)
    {
      m_actsFitResults = new std::map<const unsigned int, Trajectory>;
      PHDataNode<std::map<const unsigned int, Trajectory>> *fitNode = 
	new PHDataNode<std::map<const unsigned int, Trajectory>>(m_actsFitResults, "ActsCKFResults");
      
      svtxNode->addNode(fitNode);

    }


  return;
}

int PHActsTrkProp::getNodes(PHCompositeNode* topNode)
{

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

