#include "PHActsVertexFinder.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <ActsExamples/Plugins/BField/BFieldOptions.hpp>
#include <ActsExamples/Plugins/BField/ScalableBField.hpp>

#include <Acts/EventData/TrackParameters.hpp>
#include <Acts/MagneticField/ConstantBField.hpp>
#include <Acts/MagneticField/InterpolatedBFieldMap.hpp>
#include <Acts/MagneticField/SharedBField.hpp>
#include <Acts/Propagator/EigenStepper.hpp>
#include <Acts/Propagator/Propagator.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>
#include <Acts/Utilities/Definitions.hpp>
#include <Acts/Utilities/Helpers.hpp>
#include <Acts/Utilities/Logger.hpp>
#include <Acts/Utilities/Units.hpp>
#include <Acts/Vertexing/FullBilloirVertexFitter.hpp>
#include <Acts/Vertexing/HelicalTrackLinearizer.hpp>
#include <Acts/Vertexing/ImpactPointEstimator.hpp>
#include <Acts/Vertexing/IterativeVertexFinder.hpp>
#include <Acts/Vertexing/LinearizedTrack.hpp>
#include <Acts/Vertexing/Vertex.hpp>
#include <Acts/Vertexing/VertexFinderConcept.hpp>
#include <Acts/Vertexing/VertexingOptions.hpp>
#include <Acts/Vertexing/ZScanVertexFinder.hpp>

#include <Acts/Geometry/GeometryContext.hpp>
#include <Acts/MagneticField/MagneticFieldContext.hpp>

#include <iostream>

PHActsVertexFinder::PHActsVertexFinder(const std::string &name)
  : PHInitVertexing(name)
  , m_actsFitResults(nullptr)
  , m_event(0)
  , m_maxVertices(50)
  , m_tGeometry(nullptr)
{

}


int PHActsVertexFinder::Setup(PHCompositeNode *topNode)
{
  int ret = createNodes(topNode);
  if(ret != Fun4AllReturnCodes::EVENT_OK)
    return ret;
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsVertexFinder::Process(PHCompositeNode *topNode)
{
  if(Verbosity() > 0)
    {
      std::cout << "Starting event " << m_event << " in PHActsVertexFinder"
		<< std::endl;
    }
  
  int ret = getNodes(topNode);
  if(ret != Fun4AllReturnCodes::EVENT_OK)
    return ret;

  auto logLevel = Acts::Logging::INFO;

  if(Verbosity() > 2)
    logLevel = Acts::Logging::VERBOSE;

  /// Get the list of tracks in Acts form
  std::vector<const Acts::BoundParameters*> trackPointers = getTracks();

  auto logger = Acts::getDefaultLogger("PHActsVertexFinder", logLevel);

  /// Determine the input mag field type from the initial geometry
  /// and run the vertex finding with the determined mag field
  std::visit([&](auto &inputField) {
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
	Acts::FullBilloirVertexFitter<TrackParameters,Linearizer>;
      using ImpactPointEstimator = 
	Acts::ImpactPointEstimator<TrackParameters, Propagator>;
      using VertexSeeder = Acts::ZScanVertexFinder<VertexFitter>;
      using VertexFinder = 
	Acts::IterativeVertexFinder<VertexFitter, VertexSeeder>;
      using VertexFinderOptions = Acts::VertexingOptions<TrackParameters>;

      static_assert(Acts::VertexFinderConcept<VertexSeeder>,
		    "VertexSeeder does not fulfill vertex finder concept.");
      static_assert(Acts::VertexFinderConcept<VertexFinder>,
		    "VertexFinder does not fulfill vertex finder concept.");

      MagneticField bField(std::move(inputField));
      auto propagator = std::make_shared<Propagator>(Stepper(bField));
      PropagatorOptions propagatorOpts(m_tGeometry->geoContext,
				       m_tGeometry->magFieldContext,
				       Acts::LoggerWrapper(*logger));
      
      /// Setup vertex finder now
      typename VertexFitter::Config vertexFitterConfig;
      VertexFitter vertexFitter(std::move(vertexFitterConfig));
      
      typename Linearizer::Config linearizerConfig(bField, propagator);
      Linearizer linearizer(std::move(linearizerConfig));
      
      typename ImpactPointEstimator::Config ipEstConfig(bField, propagator);
      ImpactPointEstimator ipEst(std::move(ipEstConfig));
      
      typename VertexSeeder::Config seederConfig(ipEst);
      VertexSeeder seeder(std::move(seederConfig));
      
      typename VertexFinder::Config finderConfig(std::move(vertexFitter), std::move(linearizer),
					std::move(seeder), ipEst);
      finderConfig.maxVertices = m_maxVertices;
      finderConfig.reassignTracksAfterFirstFit = true;
      VertexFinder finder(finderConfig);
      
      typename VertexFinder::State state(m_tGeometry->magFieldContext);
      VertexFinderOptions finderOptions(m_tGeometry->geoContext,
					m_tGeometry->magFieldContext);
      
      auto result = finder.find(trackPointers, finderOptions, state);
      
      if(result.ok())
	{
	  auto vertexCollection = *result;
	  
	  if(Verbosity() > 1)
	    {
	      std::cout << "Acts IVF found " << vertexCollection.size()
			<< " vertices in event" << std::endl;
	    }
	  
	  for(const auto& vertex : vertexCollection) 
	    {
	      if(Verbosity() > 1)
		{
		  std::cout << "Found vertex at (" << vertex.position().x()
			    << ", " << vertex.position().y() << ", " 
			    << vertex.position().z() << ")" << std::endl;
		}
	    }
	}
      else
	{
	  if(Verbosity() > 1)
	    {
	      std::cout << "Acts IVF returned error: " 
			<< result.error().message() << std::endl;
	    }
	}

    } /// end lambda
    , m_tGeometry->magField
    ); /// end std::visit call

  if(Verbosity() > 1)
    std::cout << "Finished PHActsVertexFinder::process_event" << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsVertexFinder::End(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}



std::vector<const Acts::BoundParameters*> PHActsVertexFinder::getTracks()
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

int PHActsVertexFinder::createNodes(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsVertexFinder::getNodes(PHCompositeNode *topNode)
{
  
  m_actsFitResults = findNode::getClass<std::map<const unsigned int, Trajectory>>
    (topNode, "ActsTrajectories");
  if(!m_actsFitResults)
    {
      std::cout << PHWHERE << "Acts Trajectories not found on node tree, exiting."
		<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;

    }
  
  m_tGeometry = findNode::getClass<ActsTrackingGeometry>(topNode, 
							 "ActsTrackingGeometry");
  if(!m_tGeometry)
    {
      std::cout << PHWHERE << "ActsTrackingGeometry not on node tree. Exiting"
		<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }


  return Fun4AllReturnCodes::EVENT_OK;
}
