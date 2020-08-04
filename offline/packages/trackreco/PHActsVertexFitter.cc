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

#include <ACTFW/Plugins/BField/BFieldOptions.hpp>
#include <ACTFW/Plugins/BField/ScalableBField.hpp>

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

PHActsVertexFitter::PHActsVertexFitter()
{}


int PHActsVertexFitter::Init(PHCompositeNode *topNode)
{
  if(getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
    return Fun4AllReturnCodes::ABORTRUN;

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsVertexFitter::End(PHCompositeNode *topNode)
{

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsVertexFitter::process_event(PHCompositeNode *topNode)
{


  std::vector<Acts::BoundParameters> tracks =
    getTracks();

  /// Determine the input mag field type from the initial geometry created in
  /// MakeActsGeometry
  std::visit([&](auto& inputField) {
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
      
      MagneticField bField(std::move(inputField));
      auto propagator = std::make_shared<Propagator>(Stepper(bField));
      PropagatorOptions propagatorOpts(m_tGeometry->geoContext,
				       m_tGeometry->magFieldContext);
      
      typename VertexFitter::Config vertexFitterCfg;
      VertexFitter fitter(vertexFitterCfg);
      //VertexFitter::State(state(m_tGeometry->magFieldContext));
      /*
      Linearizer::Config linConfig(bField, propagator);
      Linearizer linearizer(linConfig);
      */
    }
    , m_tGeometry->magField
    );

  return Fun4AllReturnCodes::EVENT_OK;
}


std::vector<Acts::BoundParameters> PHActsVertexFitter::getTracks()
{

  std::vector<Acts::BoundParameters> tracks;
  std::map<const unsigned int, Trajectory>::iterator trackIter;

  for (trackIter = m_actsFitResults->begin();
       trackIter != m_actsFitResults->end();
       ++trackIter)
  {
    const Trajectory traj = trackIter->second;
    const auto &[trackTips, mj] = traj.trajectory();
    
    /// For the KF this iterates once. For the CKF it may iterate multiple
    /// times per trajectory
    for(const size_t &trackTip : trackTips)
      {
	if(traj.hasTrackParameters(trackTip))
	  {
	    tracks.push_back(traj.trackParameters(trackTip));
	  }

      }
  }

  return tracks;

}

int PHActsVertexFitter::getNodes(PHCompositeNode *topNode)
{
  m_actsFitResults = findNode::getClass<std::map<const unsigned int, Trajectory>>(topNode, "ActsFitResults");
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
