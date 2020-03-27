/*!
 *  \file		PHActsTrkFitter.C
 *  \brief		Refit SvtxTracks with PHActs.
 *  \details	Refit SvtxTracks with PHActs.
 *  \author	        Tony Frawley <afrawley@fsu.edu>
 */

#include "PHActsTrkFitter.h"
#include "MakeActsGeometry.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <Acts/Surfaces/Surface.hpp>
#include <Acts/Surfaces/PlaneSurface.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>
#include <Acts/EventData/TrackParameters.hpp>

#include <ACTFW/EventData/Track.hpp>
#include <ACTFW/Framework/AlgorithmContext.hpp>
#include <ACTFW/Framework/IContextDecorator.hpp>
#include <ACTFW/Framework/WhiteBoard.hpp>
#include <ACTFW/Fitting/TrkrClusterFittingAlgorithm.hpp>
#include <ACTFW/Utilities/Options.hpp>


#include <cmath>                                                // for atoi
#include <iostream>
#include <vector>


PHActsTrkFitter::PHActsTrkFitter(const std::string& name)
  : PHTrackFitting(name)
  , m_event(0)
  , m_actsProtoTracks(nullptr)
  , m_fitCfgOptions(nullptr)
{
  Verbosity(0);
}

int PHActsTrkFitter::Setup(PHCompositeNode *topNode)
{

  if( getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
    return Fun4AllReturnCodes::ABORTEVENT;

 
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsTrkFitter::Process()
{
  m_event++;

  if (Verbosity() > 1)
    {
      std::cout << PHWHERE << "Events processed: " << m_event << std::endl;
      std::cout << "Start PHActsTrkfitter::process_event" << std::endl;
    }

  /// Construct a perigee surface as the target surface (?)
  auto pSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(
      Acts::Vector3D{0.,0.,0.});

  /*

      /// Call KF now. Have a vector of sourceLinks corresponding to clusters
      /// associated to this track and the corresponding track seed which corresponds
      /// to the PHGenFitTrkProp track seeds
      Acts::KalmanFitterOptions kfOptions(context.geoContext, 
					  context.magFieldContext,
					  context.calibContext,
					  &(*pSurface));

      

      /// Run the fitter
      auto result = fitCfg.fit(trackSourceLinks, trackSeed, kfOptions);



      /// Check that the result is okay
      if(result.ok()) {
	const auto& fitOutput = result.value();
	if(fitOutput.fittedParameters){
	  const auto& params = fitOutput.fittedParameters.value();
	  /// Get position, momentum from params
	  if(Verbosity() > 10){
	    std::cout<<"Fitted parameters for track"<<std::endl;
	    std::cout<<" position : " << params.position().transpose()<<std::endl;
	    std::cout<<" momentum : " << params.momentum().transpose()<<std::endl;
	    }
	  }

	}


      /// Add a new track to a container to put on the node tree


    
  */
  return 0;
}


 
int PHActsTrkFitter::End(PHCompositeNode* topNode)
{
  if(Verbosity() > 10)
    {
      std::cout<<"Finished PHActsTrkFitter"<<std::endl;
    }
  return Fun4AllReturnCodes::EVENT_OK;
}

PHActsTrkFitter::~PHActsTrkFitter()
{

}


int PHActsTrkFitter::createNodes(PHCompositeNode* topNode)
{
  
  return Fun4AllReturnCodes::EVENT_OK;
}



int PHActsTrkFitter::getNodes(PHCompositeNode* topNode)
{
  
  m_actsProtoTracks = findNode::getClass<std::vector<ActsTrack>>(topNode, "ActsProtoTracks");
  
  if(!m_actsProtoTracks)
    {
      std::cout << "Acts proto tracks not on node tree. Exiting." 
		<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  m_fitCfgOptions = findNode::getClass<FitCfgOptions>(topNode, "ActsFitCfg");
  
  if(!m_fitCfgOptions)
    {
      std::cout << "Acts FitCfgOptions not on node tree. Exiting."
		<< std::endl;
      
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  
  return Fun4AllReturnCodes::EVENT_OK;
}



