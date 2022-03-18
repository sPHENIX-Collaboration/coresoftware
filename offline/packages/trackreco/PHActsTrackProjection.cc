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

#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackState.h>
#include <trackbase_historic/SvtxVertexMap.h>
#include <trackbase_historic/SvtxVertex.h>

#include <calobase/RawTowerGeomContainer.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTower.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawCluster.h>
#include <calobase/RawClusterUtility.h>
#include <phgeom/PHGeomUtility.h>

#include <Acts/Geometry/GeometryIdentifier.hpp>
#include <Acts/Propagator/EigenStepper.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>
#include <Acts/MagneticField/MagneticFieldProvider.hpp>

#include <CLHEP/Vector/ThreeVector.h> 
#include <math.h>

PHActsTrackProjection::PHActsTrackProjection(const std::string& name)
  : SubsysReco(name)
  , m_trajectories(nullptr)
{
  m_caloNames.push_back("CEMC");
  m_caloNames.push_back("HCALIN");
  m_caloNames.push_back("HCALOUT");

  m_caloTypes.push_back(SvtxTrack::CEMC);
  m_caloTypes.push_back(SvtxTrack::HCALIN);
  m_caloTypes.push_back(SvtxTrack::HCALOUT);
}

int PHActsTrackProjection::InitRun(PHCompositeNode *topNode)
{
  if(Verbosity() > 1)
    std::cout << "PHActsTrackProjection begin Init" << std::endl;

  int ret = makeCaloSurfacePtrs(topNode);

  if(getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
    ret = Fun4AllReturnCodes::ABORTEVENT;

  if(ret == Fun4AllReturnCodes::ABORTEVENT)
    {
      /// If calos aren't available, set a flag so that job doesn't
      /// quit processing but the flag will skip process event
      m_calosAvailable = false;
      return Fun4AllReturnCodes::EVENT_OK;
    }

  if(Verbosity() > 1)
    std::cout << "PHActsTrackProjection finished Init" << std::endl;
  
  return ret;
}

int PHActsTrackProjection::process_event(PHCompositeNode *topNode)
{
  if(Verbosity() > 1)
    {
      std::cout << "PHActsTrackProjection : Starting process_event event "
	      << m_event << std::endl;
    }

  if(!m_calosAvailable) { return Fun4AllReturnCodes::EVENT_OK; }

  for(int layer = 0; layer < m_nCaloLayers; layer++)
    {
      if(Verbosity() > 1)
	std::cout << "Processing calo layer " 
		  << m_caloNames.at(layer) << std::endl;

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

int PHActsTrackProjection::Init(PHCompositeNode */*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsTrackProjection::End(PHCompositeNode */*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}


int PHActsTrackProjection::projectTracks(PHCompositeNode *topNode, 
					 const int caloLayer)
{
  if (setCaloContainerNodes(topNode, caloLayer) 
      != Fun4AllReturnCodes::EVENT_OK)
    return Fun4AllReturnCodes::ABORTEVENT;

  for(const auto& [trackKey, traj] : *m_trajectories)
    {
      const auto track = m_trackMap->get(trackKey);
      const auto& trackTips = traj.tips();
  
      if(trackTips.size() > 1 and Verbosity() > 0)
	{ 
	  std::cout << PHWHERE 
		    << "More than 1 track tip per track. Should never happen..."
		    << std::endl;
	}
      for(const auto& trackTip : trackTips)
	{
	  const auto params = traj.trackParameters(trackTip);
	  auto cylSurf = 
	    m_caloSurfaces.find(m_caloNames.at(caloLayer))->second;
	  
	  auto result = propagateTrack(params, cylSurf);
	  
	  if(result.ok())
	    {
	      auto trackStateParams = std::move(**result);
	      updateSvtxTrack(trackStateParams, 
			      track, caloLayer);
	    }
	}
    }

  return Fun4AllReturnCodes::EVENT_OK;
}

Acts::BoundTrackParameters 
PHActsTrackProjection::makeTrackParams(SvtxTrack* track)
{

  Acts::Vector3 momentum(track->get_px(), 
			 track->get_py(), 
			 track->get_pz());
  
  auto actsVertex = getVertex(track);
  auto perigee = 
    Acts::Surface::makeShared<Acts::PerigeeSurface>(actsVertex);
  auto actsFourPos = 
    Acts::Vector4(track->get_x() * Acts::UnitConstants::cm,
		  track->get_y() * Acts::UnitConstants::cm,
		  track->get_z() * Acts::UnitConstants::cm,
		  10 * Acts::UnitConstants::ns);

  Acts::BoundSymMatrix cov;
  for(int i = 0; i < 6; i++)
    for(int j = 0; j < 6; j++)
      { cov(i,j) = track->get_acts_covariance(i,j); }

  return ActsExamples::TrackParameters::create(perigee, m_tGeometry->geoContext,
					       actsFourPos, momentum,
					       track->get_charge() / track->get_p(),
					       cov).value();

}
Acts::Vector3 PHActsTrackProjection::getVertex(SvtxTrack *track)
{
  auto vertexId = track->get_vertex_id();
  const SvtxVertex* svtxVertex = m_vertexMap->get(vertexId);

  Acts::Vector3 vertex(svtxVertex->get_x() * Acts::UnitConstants::cm, 
		       svtxVertex->get_y() * Acts::UnitConstants::cm, 
		       svtxVertex->get_z() * Acts::UnitConstants::cm);
  return vertex;
}


void PHActsTrackProjection::updateSvtxTrack(
               const Acts::BoundTrackParameters& params, 
	       SvtxTrack* svtxTrack,
	       const int caloLayer)
{
  auto projectionPos = params.position(m_tGeometry->geoContext);

  auto projectionPhi = atan2(projectionPos(1), projectionPos(0));
  auto projectionEta = asinh(projectionPos(2) / 
			     sqrt(projectionPos(0) * projectionPos(0)
				  + projectionPos(1) * projectionPos(1)));


  if(Verbosity() > 2)
    {
      double radius = sqrt(projectionPos(0) * projectionPos(0)
			   + projectionPos(1) * projectionPos(1));
      std::cout << "Track projection phi/eta " << projectionPhi
		<< " and " << projectionEta << std::endl
		<< " projection position " << projectionPos.transpose()
		<< " and radius is " << radius << std::endl;
	}

  if(fabs(projectionEta) >= 1.1)
    return;

  auto phiBin = m_towerGeomContainer->get_phibin(projectionPhi);
  auto etaBin = m_towerGeomContainer->get_etabin(projectionEta);
  
  double energy3x3 = 0.0;
  double energy5x5 = 0.0;

  getSquareTowerEnergies(phiBin, etaBin, energy3x3, energy5x5);

  if(Verbosity() > 2)
    std::cout << "3x3/5x5 energy sums " << energy3x3 
	      << " and " << energy5x5 << std::endl;

  svtxTrack->set_cal_energy_3x3(m_caloTypes.at(caloLayer),
				energy3x3);
  svtxTrack->set_cal_energy_5x5(m_caloTypes.at(caloLayer),
				energy5x5);

  double minIndex = NAN;
  double minDphi = NAN;
  double minDeta = NAN;
  double minE = NAN;
  
  getClusterProperties(projectionPhi, projectionEta,
		       minIndex, minDphi, minDeta, minE);

  if(!std::isnan(minIndex))
    {
      svtxTrack->set_cal_dphi(m_caloTypes.at(caloLayer),
			      minDphi);
      svtxTrack->set_cal_deta(m_caloTypes.at(caloLayer),
			      minDeta);
      svtxTrack->set_cal_cluster_id(m_caloTypes.at(caloLayer),
				    minIndex);
      svtxTrack->set_cal_cluster_e(m_caloTypes.at(caloLayer),
				   minE);
      if(Verbosity() > 2)
	std::cout << "Calo cluster has dphi " << minDphi
		  << " and deta " << minDeta << " and clusID "
		  << minIndex << " and energy " << minE 
		  << std::endl;

    }

  return;
}

void PHActsTrackProjection::getClusterProperties(double phi,
						 double eta,
						 double& minIndex,
						 double& minDphi,
						 double& minDeta,
						 double& minE)
{
  double minR = DBL_MAX;
  auto clusterMap = m_clusterContainer->getClustersMap();
  for(const auto& [key, cluster] : clusterMap)
    {
      const auto clusterEta = 
	RawClusterUtility::GetPseudorapidity(*cluster, 
					     CLHEP::Hep3Vector(0,0,0));
      const auto dphi = deltaPhi(phi - cluster->get_phi());
      const auto deta = eta - clusterEta;
      const auto r = sqrt(pow(dphi,2) + pow(deta,2));

      if(r < minR)
	{
	  minIndex = key;
	  minR = r;
	  minDphi = dphi;
	  minDeta = deta;
	  minE = cluster->get_energy();
	}
    }

  return;
}

void PHActsTrackProjection::getSquareTowerEnergies(int phiBin, 
						   int etaBin,
						   double& energy3x3,
						   double& energy5x5)
{
  for(int iphi = phiBin - 2; iphi <= phiBin + 2; ++iphi) 
    {
      for(int ieta = etaBin - 2; ieta <= etaBin + 2; ++ieta)
	{
	  /// Check the phi periodic boundary conditions
	  int wrapPhi = iphi;
	  if(wrapPhi < 0)
	    wrapPhi += m_towerGeomContainer->get_phibins();
	  if(wrapPhi >= m_towerGeomContainer->get_phibins())
	    wrapPhi -= m_towerGeomContainer->get_phibins();
	  
	  /// Check the eta boundary conditions
	  if (ieta < 0 or ieta >= m_towerGeomContainer->get_etabins())
	    continue;
	  
	  auto tower = m_towerContainer->getTower(ieta, wrapPhi);
	  
	  if(!tower)
	    continue;
	  
	  energy5x5 += tower->get_energy();
	  if(abs(iphi - phiBin) <= 1 and abs(ieta - etaBin) <= 1)
	    energy3x3 += tower->get_energy();
	  
	}
    }

  return;
}

BoundTrackParamPtrResult PHActsTrackProjection::propagateTrack(
        const Acts::BoundTrackParameters& params, 
	const SurfacePtr& targetSurf)
{
  
  if(Verbosity() > 1) {
    std::cout << "Propagating final track fit with momentum: " 
	      << params.momentum() << " and position " 
	      << params.position(m_tGeometry->geoContext)
	      << std::endl
	      << "track fit phi/eta "
	      << atan2(params.momentum()(1), 
		       params.momentum()(0)) 
	      << " and " 
	      << atanh(params.momentum()(2) 
		       / params.momentum().norm())
	      << std::endl;
  }

  using Stepper = Acts::EigenStepper<>;
  using Propagator = Acts::Propagator<Stepper>;

  Stepper stepper(m_tGeometry->magField);
  Propagator propagator(stepper);

  
  Acts::Logging::Level logLevel = Acts::Logging::FATAL;
  if(Verbosity() > 3)
    { logLevel = Acts::Logging::VERBOSE; }
  
  auto logger = Acts::getDefaultLogger("PHActsTrackProjection", 
				       logLevel);
  
  Acts::PropagatorOptions<> options(m_tGeometry->geoContext,
				    m_tGeometry->magFieldContext,
				    Acts::LoggerWrapper{*logger});
  
  auto result = propagator.propagate(params, *targetSurf, 
				     options);

  if(result.ok())
    { return std::move((*result).endParameters); }
  
  return result.error();

}

int PHActsTrackProjection::setCaloContainerNodes(PHCompositeNode *topNode,
						  const int caloLayer)
{
  std::string towerGeoNodeName = "TOWERGEOM_" + m_caloNames.at(caloLayer);
  std::string towerNodeName    = "TOWER_CALIB_" + m_caloNames.at(caloLayer);
  std::string clusterNodeName  = "CLUSTER_" + m_caloNames.at(caloLayer);

  m_towerGeomContainer = findNode::getClass<RawTowerGeomContainer>
    (topNode, towerGeoNodeName.c_str());

  m_towerContainer = findNode::getClass<RawTowerContainer>
    (topNode, towerNodeName.c_str());
  
  m_clusterContainer = findNode::getClass<RawClusterContainer>
    (topNode, clusterNodeName.c_str());

  if(m_useCemcPosRecalib and 
     m_caloNames.at(caloLayer).compare("CEMC") == 0) 
    {
      std::string nodeName = "CLUSTER_POS_COR_" + m_caloNames.at(caloLayer);
      m_clusterContainer = findNode::getClass<RawClusterContainer>
	(topNode, nodeName.c_str());
    }
  
  if(!m_towerGeomContainer or !m_towerContainer or !m_clusterContainer)
    {
      std::cout << PHWHERE 
		<< "Calo geometry and/or cluster container not found on node tree. Track projections to calos won't be filled."
		<< std::endl;
      m_calosAvailable = false;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsTrackProjection::makeCaloSurfacePtrs(PHCompositeNode *topNode)
{
  for(int caloLayer = 0; caloLayer < m_nCaloLayers; caloLayer++)
    {
      if(setCaloContainerNodes(topNode, caloLayer) != Fun4AllReturnCodes::EVENT_OK)
	return Fun4AllReturnCodes::ABORTEVENT;
      
      /// Default to using calo radius
      double caloRadius = m_towerGeomContainer->get_radius() 
	* Acts::UnitConstants::cm;
      if(m_caloRadii.find(m_caloTypes.at(caloLayer)) != m_caloRadii.end())
	{ 
	  caloRadius = m_caloRadii.find(m_caloTypes.at(caloLayer))->second
	    * Acts::UnitConstants::cm; 
	}
    
      /// Extend farther so that there is at least surface there, for high
      /// curling tracks. Can always reject later
      const auto eta = 2.5;
      const auto theta = 2. * atan(exp(-eta));
      const auto halfZ = caloRadius / tan(theta) 
	* Acts::UnitConstants::cm;
      
      /// Make a cylindrical surface at (0,0,0) aligned along the z axis
      auto transform = Acts::Transform3::Identity();

      std::shared_ptr<Acts::CylinderSurface> surf = 
	Acts::Surface::makeShared<Acts::CylinderSurface>(transform,
							 caloRadius,
							 halfZ);
  
      m_caloSurfaces.insert(std::make_pair(m_caloNames.at(caloLayer),
					   surf));
    }

  if(Verbosity() > 1)
    {
      for(const auto& [name, surfPtr] : m_caloSurfaces)
	{
	  std::cout << "Cylinder " << name << " has center "
		    << surfPtr.get()->center(m_tGeometry->geoContext)
		    << std::endl;
	}
    }

  return Fun4AllReturnCodes::EVENT_OK;

}

int PHActsTrackProjection::getNodes(PHCompositeNode *topNode)
{
  m_trajectories = findNode::getClass<std::map<const unsigned int, Trajectory>>(topNode, "ActsTrajectories");
  if(!m_trajectories)
    {
      std::cout << PHWHERE << "No Acts trajectories on node tree, bailing."
		<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  m_vertexMap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
  if(!m_vertexMap)
    {
      std::cout << PHWHERE << "No vertex map on node tree, bailing."
		<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  m_tGeometry = findNode::getClass<ActsTrackingGeometry>(
			  topNode, "ActsTrackingGeometry");
  if(!m_tGeometry)
    {
      std::cout << "ActsTrackingGeometry not on node tree. Exiting."
		<< std::endl;
      
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  
  m_trackMap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if(!m_trackMap)
    {
      std::cout << PHWHERE << "No SvtxTrackMap on node tree. Bailing."
		<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  return Fun4AllReturnCodes::EVENT_OK;
}

double PHActsTrackProjection::deltaPhi(const double& phi)
{
  if (phi > M_PI) 
    return phi - 2. * M_PI;
  else if (phi <= -M_PI) 
    return phi + 2.* M_PI;
  else 
    return phi;
}
