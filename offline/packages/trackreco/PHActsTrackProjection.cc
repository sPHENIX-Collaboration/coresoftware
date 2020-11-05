#include "PHActsTrackProjection.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/PHDataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/PHTimer.h>

#include <calobase/RawTowerGeomContainer.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTower.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawCluster.h>
#include <calobase/RawClusterUtility.h>
#include <phgeom/PHGeomUtility.h>

#include <Acts/Geometry/GeometryIdentifier.hpp>
#include <Acts/MagneticField/ConstantBField.hpp>
#include <Acts/MagneticField/InterpolatedBFieldMap.hpp>
#include <Acts/MagneticField/SharedBField.hpp>
#include <Acts/Propagator/EigenStepper.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>

#include <ActsExamples/Plugins/BField/ScalableBField.hpp>

PHActsTrackProjection::PHActsTrackProjection(const std::string& name)
  : SubsysReco(name)
{
  m_caloNames[0] = "CEMC";
  m_caloNames[1] = "HCALIN";
  m_caloNames[2] = "HCALOUT";
}
int PHActsTrackProjection::makeCaloSurfacePtrs(PHCompositeNode *topNode)
{
  for(int calLayer = 0; calLayer < m_nCaloLayers; calLayer++)
    {
      if(setCaloContainerNodes(topNode, calLayer) != Fun4AllReturnCodes::EVENT_OK)
	return Fun4AllReturnCodes::ABORTEVENT;
      
      float caloRadius = m_towerGeomContainer->get_radius();

      //std::shared_ptr<Acts::CylinderSurface> surf = 
      //Acts::Surface::makeShared<Acts::CylinderSurface>(args);
      /// Just put a placeholder surf in here for now while figuring
      /// out the cylinder args
      std::shared_ptr<Acts::PerigeeSurface> surf = 
	Acts::Surface::makeShared<Acts::PerigeeSurface>(
        Acts::Vector3D{0.,0.,caloRadius});
      m_caloSurfaces.insert(std::make_pair(m_caloNames[calLayer],
					   surf));

    }

  return Fun4AllReturnCodes::EVENT_OK;

}
int PHActsTrackProjection::Init(PHCompositeNode *topNode)
{
  int ret = makeCaloSurfacePtrs(topNode);
 

  return ret;
}

int PHActsTrackProjection::process_event(PHCompositeNode *topNode)
{
  if(Verbosity() > 1)
    std::cout << "PHActsTrackProjection : Starting process_event event "
	      << m_event << std::endl;

  for(int layer = 0; layer < m_nCaloLayers; layer++)
    {
      int ret = projectTracks(topNode, layer);
      if(ret != Fun4AllReturnCodes::EVENT_OK)
	return Fun4AllReturnCodes::ABORTEVENT;
    }


  if(Verbosity() > 1)
    std::cout << "PHActsTrackProjection : Finished process_event event "
	      << m_event << std::endl;

  m_event++;

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsTrackProjection::End(PHCompositeNode *topNode)
{
  
  return Fun4AllReturnCodes::EVENT_OK;
}


int PHActsTrackProjection::projectTracks(PHCompositeNode *topNode, 
					 const int calLayer)
{
  if(setCaloContainerNodes(topNode, calLayer) != Fun4AllReturnCodes::EVENT_OK)
    return Fun4AllReturnCodes::ABORTEVENT;
  
  std::map<const unsigned int, Trajectory>::iterator trajIter;
  for(trajIter = m_actsFitResults->begin();
      trajIter != m_actsFitResults->end();
      ++trajIter)
    {
      const auto trackKey = trajIter->first;
      const auto traj = trajIter->second;

      const auto &[trackTips, mj] = traj.trajectory();

      /// Skip failed track fits
      if(trackTips.empty())
	continue;

      for(const size_t& trackTip : trackTips)
	{
	  if(traj.hasTrackParameters(trackTip))
	    {
	      const auto &params = traj.trackParameters(trackTip);
	      auto cylSurf = 
		m_caloSurfaces.find(m_caloNames[calLayer])->second;
	      auto result = propagateTrack(params, cylSurf);

	      if(result.ok())
		{
		  auto trackStateParams = std::move(**result);
		  updateSvtxTrack(trackStateParams, trackKey);
		}

	    }
	}

    }

  return Fun4AllReturnCodes::EVENT_OK;
}

void PHActsTrackProjection::updateSvtxTrack(const Acts::BoundTrackParameters& params, const unsigned int trackKey)
{
  

}

BoundTrackParamPtrResult PHActsTrackProjection::propagateTrack(
        const FitParameters& params, 
	const SurfacePtr& targetSurf)
{
  
  if(Verbosity() > 1)
    std::cout << "Propagating silicon+MM fit params momentum: " 
	      << params.momentum() << " and position " 
	      << params.position(m_tGeometry->geoContext)
	      << std::endl;

  return std::visit([params, targetSurf, this]
		    (auto && inputField) -> BoundTrackParamPtrResult {
      using InputMagneticField = 
	typename std::decay_t<decltype(inputField)>::element_type;
      using MagneticField      = Acts::SharedBField<InputMagneticField>;
      using Stepper            = Acts::EigenStepper<MagneticField>;
      using Propagator         = Acts::Propagator<Stepper>;

      MagneticField field(inputField);
      Stepper stepper(field);
      Propagator propagator(stepper);

      Acts::Logging::Level logLevel = Acts::Logging::FATAL;
      if(Verbosity() > 3)
	logLevel = Acts::Logging::VERBOSE;

      auto logger = Acts::getDefaultLogger("PHTpcResiduals", logLevel);
      
      Acts::PropagatorOptions<> options(m_tGeometry->geoContext,
					m_tGeometry->magFieldContext,
					Acts::LoggerWrapper{*logger});
     
      auto result = propagator.propagate(params, *targetSurf, 
					 options);
   
      if(result.ok())
	return std::move((*result).endParameters);
      else
	return result.error();
   },
     std::move(m_tGeometry->magField));

}

int PHActsTrackProjection::setCaloContainerNodes(PHCompositeNode *topNode,
						  const int calLayer)
{
  std::string towerGeoNodeName = "TOWERGEOM_" + m_caloNames[calLayer];
  std::string towerNodeName    = "TOWER_CALIB_" + m_caloNames[calLayer];
  std::string clusterNodeName  = "CLUSTER_" + m_caloNames[calLayer];

  m_towerGeomContainer = findNode::getClass<RawTowerGeomContainer>
    (topNode, towerGeoNodeName.c_str());

  m_towerContainer = findNode::getClass<RawTowerContainer>
    (topNode, towerNodeName.c_str());
  
  m_clusterContainer = findNode::getClass<RawClusterContainer>
    (topNode, clusterNodeName.c_str());

  if(!m_towerGeomContainer or !m_towerContainer or !m_clusterContainer)
    {
      std::cout << PHWHERE 
		<< "Calo geometry and/or cluster container not found on node tree. Exiting"
		<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsTrackProjection::getNodes(PHCompositeNode *topNode)
{
  m_tGeometry = findNode::getClass<ActsTrackingGeometry>(topNode, "ActsTrackingGeometry");
  if(!m_tGeometry)
    {
      std::cout << "ActsTrackingGeometry not on node tree. Exiting."
		<< std::endl;
      
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  
    m_actsFitResults = findNode::getClass<std::map<const unsigned int, Trajectory>>
                     (topNode, "ActsFitResults");

  if (!m_actsFitResults)
  {
    std::cout << PHWHERE << "No Acts fit results on node tree. Bailing."
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
int PHActsTrackProjection::createNodes(PHCompositeNode *topNode)
{
  
  return Fun4AllReturnCodes::EVENT_OK;
}
