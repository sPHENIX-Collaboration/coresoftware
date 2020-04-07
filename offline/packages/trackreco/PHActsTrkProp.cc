
#include "PHActsTrkProp.h"
#include "PHActsTracks.h"

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

#include <Acts/MagneticField/InterpolatedBFieldMap.hpp>
#include <Acts/MagneticField/SharedBField.hpp>
#include <Acts/MagneticField/ConstantBField.hpp>

#include <Acts/EventData/ChargePolicy.hpp>
#include <Acts/EventData/SingleCurvilinearTrackParameters.hpp>
#include <Acts/EventData/TrackParameters.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>

#include <Acts/Propagator/EigenStepper.hpp> // this include causes seg fault when put in header file for some reason
#include <Acts/Propagator/Propagator.hpp>
#include <Acts/Propagator/Navigator.hpp>
#include <Acts/Propagator/AbortList.hpp>
#include <Acts/Propagator/ActionList.hpp>
#include <Acts/Utilities/Helpers.hpp>
#include <Acts/Utilities/Units.hpp>

#include <ACTFW/Plugins/BField/BFieldOptions.hpp>
#include <ACTFW/Plugins/BField/ScalableBField.hpp>
#include <ACTFW/Framework/ProcessCode.hpp>
#include <ACTFW/Framework/WhiteBoard.hpp>
#include <ACTFW/EventData/Track.hpp>
#include <ACTFW/Framework/AlgorithmContext.hpp>

/// Setup aliases for creating propagator
/// For some reason putting these in the header file, with appropriate headers
/// causes seg fault. Propagator also must be instantiated immediately
using ConstantBField = Acts::ConstantBField;
using Stepper = Acts::EigenStepper<ConstantBField>;
using Propagator = Acts::Propagator<Stepper, Acts::Navigator>;
Propagator *propagator;

#include <cmath>
#include <iostream>
#include <vector>
#include <utility>

PHActsTrkProp::PHActsTrkProp(const std::string& name)
  : PHTrackPropagating(name)
  , m_event(0)
  , m_actsGeometry(nullptr)
  , m_minTrackPt(0.15)
  , m_maxStepSize(3.)
  , m_actsTracks(nullptr)
{
  Verbosity(0);
}

 PHActsTrkProp::~PHActsTrkProp()
{
}

int PHActsTrkProp::Setup(PHCompositeNode* topNode)
{
  if (getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
    return Fun4AllReturnCodes::ABORTEVENT;

  /// Get the magnetic field and tracking geometry to setup the Acts::Stepper
  //FW::Options::BFieldVariant bFieldVar = m_actsGeometry->magField;
   Acts::Navigator navigator(m_actsGeometry->tGeometry);

  /// Just use the default magnetic field for now. Can access BField
  /// from m_actsGeometry using std::visit, but can't figure out how to 
  /// get necessary information out of lambda function
  ConstantBField bField (0, 0, 1.4 * Acts::UnitConstants::T);
  Stepper stepper(bField);
  propagator = new Propagator(std::move(stepper), std::move(navigator));
  //Propagator propagator(std::move(stepper), std::move(navigator));

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

  PerigeeSurface surface = Acts::Surface::makeShared<Acts::PerigeeSurface>
    (Acts::Vector3D(0., 0., 0.));

  
  for(std::vector<ActsTrack>::iterator trackIter = m_actsTracks->begin();
      trackIter != m_actsTracks->end();
      ++trackIter)
	
    {
      ActsTrack track = *trackIter;
     
      const FW::TrackParameters actsTrack = track.trackParams;
      
      PropagationOutput pOutput = propagate(actsTrack);
      
      std::vector<Acts::detail::Step> steps = pOutput.first;
      for (auto& step : steps){
	/// Get kinematic information of steps
	/// global x,y,z
	float x = step.position.x();
	float y = step.position.y();
	float z = step.position.z();

	auto direction = step.momentum.normalized();
	///global direction x,y,yz
	float dx = direction.x();
	float dy = direction.y();
	float dz = direction.z();

	if(Verbosity() > 1)
	  {
	    std::cout << "Acts track propagation step : "
		      << x << ", " << y << ", " << z 
		      << " and momentum direction " << dx
		      << ", " << dy << ", " << dz << std::endl;
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


PropagationOutput PHActsTrkProp::propagate(FW::TrackParameters parameters)
{
  
  PropagationOutput pOutput;
  
  /// Set Acts namespace aliases to reduce code clutter
  using MaterialInteractor = Acts::MaterialInteractor;
  using SteppingLogger     = Acts::detail::SteppingLogger;
  using DebugOutput        = Acts::detail::DebugOutputActor;
  using EndOfWorld         = Acts::detail::EndOfWorldReached;
  
  using ActionList
    = Acts::ActionList<SteppingLogger, MaterialInteractor, DebugOutput>;
  using AbortList         = Acts::AbortList<EndOfWorld>;
  using PropagatorOptions = Acts::PropagatorOptions<ActionList, AbortList>;
  
  /// Setup propagator options to hold context, and propagator characteristics
  PropagatorOptions options(m_actsGeometry->geoContext, 
			    m_actsGeometry->magFieldContext);
  options.pathLimit = std::numeric_limits<double>::max();
  options.debug     = true;
  
  /// Activate loop protection at some pt value
  options.loopProtection
    = (Acts::VectorHelpers::perp(parameters.momentum())
       < m_minTrackPt);

  /// Switch the material interaction on/off & eventually into logging mode
  /// Should all of these switches be configurable from e.g. constructor?
  auto& mInteractor = options.actionList.get<MaterialInteractor>();
  mInteractor.multipleScattering = true;
  mInteractor.energyLoss         = true;
  mInteractor.recordInteractions = true;

  /// Set the maximum step size
  options.maxStepSize = m_maxStepSize * Acts::UnitConstants::mm;
  
  /// Propagate using Acts::Propagator
  const auto& result
    = propagator->propagate(parameters, options).value();
  auto steppingResults = result.template get<SteppingLogger::result_type>();
  
  /// Set the stepping result
  pOutput.first = std::move(steppingResults.steps);

  /// Also set the material recording result - if configured
  if (mInteractor.recordInteractions) {
    auto materialResult
      = result.template get<MaterialInteractor::result_type>();
    pOutput.second = std::move(materialResult);
  }
  
  if(Verbosity() > 1)
    {
      auto& debugResult = result.template get<DebugOutput::result_type>();
      std::cout << "Acts Propagator Debug: " << debugResult.debugString 
		<< std::endl;
    }
  
  return pOutput;
 
}




int PHActsTrkProp::createNodes(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsTrkProp::getNodes(PHCompositeNode* topNode)
{

  m_actsGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");

  if (!m_actsGeometry)
  {
    std::cout << "ActsGeometry not on node tree. Exiting."
              << std::endl;

    return Fun4AllReturnCodes::ABORTEVENT;
  }

  /// TODO - change the name of this to reflect the new node put on the node
  /// tree from PHActsTrkFitter. Right now for testing, using ActsProtoTracks
  m_actsTracks = findNode::getClass<std::vector<ActsTrack>>(topNode, "ActsProtoTracks");
  
  if (!m_actsTracks)
    {
      std::cout << "ActsTracks not on node tree. Exiting."

		<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }


  return Fun4AllReturnCodes::EVENT_OK;
}

