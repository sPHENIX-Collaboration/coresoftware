
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
  , m_tGeometry(nullptr)
  , m_trackMap(nullptr)
  , m_actsProtoTracks(nullptr)
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
  createNodes(topNode);
  
  if (getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
    return Fun4AllReturnCodes::ABORTEVENT;

  /// Need to implement different chi2 criteria for the different layers
  /// and volumes to help the CKF out a little 
  /// First need to walk through layers/volumes to figure out appropriate
  /// numerical identifiers
  /// m_sourceLinkSelectorConfig.layerMaxChi2 = {{2, {{2, 8}, {4, 7}}}};
  /// m_sourceLinkSelectorConfig.volumeMaxChi2 = {{2, 7}, {3, 8}};
  m_sourceLinkSelectorConfig.maxChi2 = 8;
  // Set the allowed maximum number of source links to be large enough
  m_sourceLinkSelectorConfig.maxNumSourcelinksOnSurface = 100;
 
  findCfg.finder = FW::TrkrClusterFindingAlgorithm::makeFinderFunction(
                   m_tGeometry->tGeometry,
		   m_tGeometry->magField,
		   Acts::Logging::VERBOSE);



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
  

  for (SvtxTrackMap::Iter trackIter = m_trackMap->begin();
       trackIter != m_trackMap->end(); ++trackIter)
    {
      const SvtxTrack *track = trackIter->second;
      
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
      
      // just set to 40 ns for now?
      const double trackTime = 40 * Acts::UnitConstants::ns;
      const int trackQ = track->get_charge();
      
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
      auto result = findCfg.finder(sourceLinks, trackSeed, ckfOptions);
      
      if(result.ok())
	{
	  const auto& fitOutput = result.value();
	  auto parameterMap = fitOutput.fittedParameters;

	  /// how to get the associated source links from fit result?
	  std::vector<size_t> allSourceLinks = fitOutput.sourcelinkCandidateIndices;  
	  std::vector<SourceLink> trackSourceLinks;
	      
	  for(size_t i = 0; i < allSourceLinks.size(); ++i)
	    {
	      trackSourceLinks.push_back(sourceLinks.at(allSourceLinks.at(i)));
	    }
	
	  for(auto element : parameterMap)
	    {
	      const auto& params = element.second;
	      if(Verbosity() > 10)
		{
		  std::cout << "Fitted parameters for track finder" << std::endl;
		  std::cout << "Position : " << params.position().transpose()
			    << std::endl;
		  std::cout << "Momentum : " << params.momentum().transpose()
			    << std::endl;

		}


	      /// Get the finder results into a FW::TrackParameters 
	      const FW::TrackParameters trackSeed(element.second.covariance(),
						  element.second.position(),
						  element.second.momentum(),
						  element.second.charge(),
						  element.second.time());
	      
	      ActsTrack actsProtoTrack(trackSeed, trackSourceLinks);
	      m_actsProtoTracks->push_back(actsProtoTrack);
	    }
	  
	}
    }


  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsTrkProp::End()
{
  if (Verbosity() > 10)
  {
    std::cout << "Finished PHActsTrkProp" << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
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


  m_actsProtoTracks = findNode::getClass<std::vector<ActsTrack>>(topNode, "ActsProtoTracks");

  if(!m_actsProtoTracks)
    {
      m_actsProtoTracks = new std::vector<ActsTrack>;

      PHDataNode<std::vector<ActsTrack>> *protoTrackNode =
        new PHDataNode<std::vector<ActsTrack>>(m_actsProtoTracks, "ActsProtoTracks");
      
      svtxNode->addNode(protoTrackNode);

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

