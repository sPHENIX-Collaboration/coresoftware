#include "PHActsVertexFitter.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <phool/phool.h>

/// Tracking includes
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxVertexMap.h>
#include <trackbase_historic/SvtxVertex.h>

#include <ActsExamples/Plugins/BField/BFieldOptions.hpp>
#include <ActsExamples/Plugins/BField/ScalableBField.hpp>

#include <Acts/MagneticField/ConstantBField.hpp>
#include <Acts/MagneticField/InterpolatedBFieldMap.hpp>
#include <Acts/MagneticField/SharedBField.hpp>
#include <Acts/EventData/TrackParameters.hpp>
#include <Acts/MagneticField/ConstantBField.hpp>
#include <Acts/Propagator/EigenStepper.hpp>
#include <Acts/Propagator/Propagator.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>
#include <Acts/Utilities/Definitions.hpp>
#include <Acts/Utilities/Helpers.hpp>
#include <Acts/Vertexing/FullBilloirVertexFitter.hpp>
#include <Acts/Vertexing/HelicalTrackLinearizer.hpp>
#include <Acts/Vertexing/LinearizedTrack.hpp>
#include <Acts/Vertexing/Vertex.hpp>
#include <Acts/Vertexing/VertexingOptions.hpp>

#include <iostream>

PHActsVertexFitter::PHActsVertexFitter(const std::string& name) 
: SubsysReco(name)
  , m_actsFitResults(nullptr)
  , m_event(0)
  , m_tGeometry(nullptr)
{
}


int PHActsVertexFitter::Init(PHCompositeNode *topNode)
{
  if(Verbosity() > 1)
    std::cout << "PHActsVertexFitter::Init" << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsVertexFitter::End(PHCompositeNode *topNode)
{
  if(Verbosity() > 1)
    std::cout << "PHActsVertexFitter::End " << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsVertexFitter::process_event(PHCompositeNode *topNode)
{

  if(getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
    return Fun4AllReturnCodes::ABORTRUN;

  auto logLevel = Acts::Logging::INFO;
  if(Verbosity() > 1)
    {
      std::cout << "Beginning PHActsVertexFitter::process_event number " 
		<< m_event << std::endl;
      logLevel = Acts::Logging::VERBOSE;
    }
  
  std::vector<const Acts::BoundParameters*> tracks = getTracks();
  
  auto logger = Acts::getDefaultLogger("PHActsVertexFitter", logLevel);

  /// Determine the input mag field type from the initial geometry created in
  /// MakeActsGeometry
  std::visit([&](auto& inputField) {

      /// Setup aliases
      using InputMagneticField = 
	typename std::decay_t<decltype(inputField)>::element_type;
      using MagneticField = Acts::SharedBField<InputMagneticField>;
      using Stepper = Acts::EigenStepper<MagneticField>;
      using Propagator = Acts::Propagator<Stepper>;
      using PropagatorOptions = Acts::PropagatorOptions<>;
      using TrackParameters = Acts::BoundParameters;
      using Linearizer = Acts::HelicalTrackLinearizer<Propagator>;
      using VertexFitter =
	Acts::FullBilloirVertexFitter<TrackParameters, Linearizer>;
      using VertexFitterOptions = Acts::VertexingOptions<TrackParameters>;
      
      /// Create necessary templated inputs for Acts vertex fitter
      MagneticField bField(std::move(inputField));
      auto propagator = std::make_shared<Propagator>(Stepper(bField));
      PropagatorOptions propagatorOpts(m_tGeometry->geoContext,
				       m_tGeometry->magFieldContext,
				       Acts::LoggerWrapper(*logger));
      
      typename VertexFitter::Config vertexFitterCfg;
      VertexFitter fitter(vertexFitterCfg);
      typename VertexFitter::State state(m_tGeometry->magFieldContext);
    
      typename Linearizer::Config linConfig(bField, propagator);
      Linearizer linearizer(linConfig);

      /// Can add a vertex fitting constraint as an option, if desired
      VertexFitterOptions vfOptions(m_tGeometry->geoContext,
				    m_tGeometry->magFieldContext);
      
      /// Call the fitter and get the result
      auto fitRes = fitter.fit(tracks,linearizer,
      		       vfOptions, state);

      if(fitRes.ok())
	{
	  Acts::Vertex<TrackParameters> fittedVertex;
	  fittedVertex = *fitRes;
	  if(Verbosity() > 3)
	    {
	      std::cout << "Fitted vertex position "
			<< fittedVertex.position().x() 
			<< ", " 
			<< fittedVertex.position().y() 
			<< ", "
			<< fittedVertex.position().z() 
			<< std::endl;

	    }
	}
      else
	{
	  if(Verbosity() > 3)
	    {
	      std::cout << "Acts vertex fit error: " 
			<< fitRes.error().message()
			<< std::endl;
	    }
	}
      
    }
    , m_tGeometry->magField
    ); /// end std::visit call

  if(Verbosity() > 1)
    std::cout << "Finished PHActsVertexFitter::process_event" 
	      << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}


std::vector<const Acts::BoundParameters*> PHActsVertexFitter::getTracks()
{
 
  std::vector<const Acts::BoundParameters*> trackPtrs;

  std::map<const unsigned int, Trajectory>::iterator trackIter;

  for (trackIter = m_actsFitResults->begin();
       trackIter != m_actsFitResults->end();
       ++trackIter)
  {
    const Trajectory traj = trackIter->second;
    const auto &[trackTips, mj] = traj.trajectory();
    
    for(const size_t &trackTip : trackTips)
      {
	if(traj.hasTrackParameters(trackTip))
	  {
	    const Acts::BoundParameters *param = new Acts::BoundParameters(traj.trackParameters(trackTip));
	 
	    trackPtrs.push_back(param);
	  }

      }
  }
  
  if(Verbosity() > 3)
    {
      std::cout << "Fitting a vertex for the following number of tracks "
		<< trackPtrs.size()
		<< std::endl;
      
      for(std::vector<const Acts::BoundParameters*>::iterator it = trackPtrs.begin();
	  it != trackPtrs.end(); ++it)
	{
	  const Acts::BoundParameters* param = *it;
	  std::cout << "Track position: (" 
		    << param->position(m_tGeometry->geoContext)(0)
		    <<", " << param->position(m_tGeometry->geoContext)(1) << ", "
		    << param->position(m_tGeometry->geoContext)(2) << ")" 
		    << std::endl;

	}
      
    }

  return trackPtrs;

}

int PHActsVertexFitter::getNodes(PHCompositeNode *topNode)
{
  m_actsFitResults = findNode::getClass<std::map<const unsigned int, Trajectory>>(topNode, "ActsTrajectories");
  if(!m_actsFitResults)
    {
      std::cout << PHWHERE << "Acts Trajectories not found on node tree, exiting."
		<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;

    }


  m_tGeometry = findNode::getClass<ActsTrackingGeometry>(topNode, "ActsTrackingGeometry");
  if(!m_tGeometry)
    {
      std::cout << PHWHERE << "ActsTrackingGeometry not on node tree. Exiting"
		<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  return Fun4AllReturnCodes::EVENT_OK;

}
