#include "PHActsTrackProjection.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/PHTimer.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <trackbase/ActsTrackFittingAlgorithm.h>

#include <trackbase_historic/ActsTransformations.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackState.h>
#include <trackbase_historic/SvtxTrackState_v1.h>
#include <trackbase_historic/SvtxVertex.h>
#include <trackbase_historic/SvtxVertexMap.h>

#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawClusterUtility.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainerv1.h>
#include <calobase/RawTowerGeomContainer.h>
#include <phgeom/PHGeomUtility.h>

#include <Acts/Geometry/GeometryIdentifier.hpp>
#include <Acts/MagneticField/ConstantBField.hpp>
#include <Acts/MagneticField/MagneticFieldProvider.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>

#include <CLHEP/Vector/ThreeVector.h>
#include <math.h>

PHActsTrackProjection::PHActsTrackProjection(const std::string& name)
  : SubsysReco(name)
{
  m_caloNames.push_back("CEMC");
  m_caloNames.push_back("HCALIN");
  m_caloNames.push_back("HCALOUT");
  m_caloNames.push_back("OUTER_CEMC");
  m_caloNames.push_back("OUTER_HCALIN");
  m_caloNames.push_back("OUTER_HCALOUT");

  m_caloTypes.push_back(SvtxTrack::CEMC);
  m_caloTypes.push_back(SvtxTrack::HCALIN);
  m_caloTypes.push_back(SvtxTrack::HCALOUT);
  m_caloTypes.push_back(SvtxTrack::OUTER_CEMC);
  m_caloTypes.push_back(SvtxTrack::OUTER_HCALIN);
  m_caloTypes.push_back(SvtxTrack::OUTER_HCALOUT);
}

int PHActsTrackProjection::InitRun(PHCompositeNode* topNode)
{
  if (Verbosity() > 1)
  {
    std::cout << "PHActsTrackProjection begin Init" << std::endl;
  }

  int ret = makeCaloSurfacePtrs(topNode);

  /// only need to print warning once
  m_calosAvailable = false;

  if (getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
  {
    ret = Fun4AllReturnCodes::ABORTEVENT;
  }

  if (Verbosity() > 1)
  {
    std::cout << "PHActsTrackProjection finished Init" << std::endl;
  }

  return ret;
}

int PHActsTrackProjection::process_event(PHCompositeNode* topNode)
{
  if (Verbosity() > 1)
  {
    std::cout << "PHActsTrackProjection : Starting process_event event "
              << m_event << std::endl;
  }

  for (int layer = 0; layer < m_nCaloLayers; layer++)
  {
    if (Verbosity() > 1)
    {
      std::cout << "Processing calo layer "
                << m_caloNames.at(layer) << std::endl;
    }

    if (setCaloContainerNodes(topNode, layer) != Fun4AllReturnCodes::EVENT_OK)
    {
      continue;
    }

    int ret = projectTracks(layer);
    if (ret != Fun4AllReturnCodes::EVENT_OK)
    {
      return Fun4AllReturnCodes::ABORTEVENT;
    }
    ret = projectTracks(layer+m_nCaloLayers);
    if (ret != Fun4AllReturnCodes::EVENT_OK)
    {
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  if (Verbosity() > 1)
  {
    std::cout << "PHActsTrackProjection : Finished process_event event "
              << m_event << std::endl;
  }

  m_event++;

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsTrackProjection::Init(PHCompositeNode* /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsTrackProjection::End(PHCompositeNode* /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsTrackProjection::projectTracks(const int caloLayer)
{
  ActsPropagator prop(m_tGeometry);
  for (const auto& [key, track] : *m_trackMap)
  {
    const auto params = prop.makeTrackParams(track, m_vertexMap);
    auto cylSurf =
        m_caloSurfaces.find(m_caloNames.at(caloLayer))->second;

    auto result = propagateTrack(params, caloLayer, cylSurf);
    if (result.ok())
    {
      updateSvtxTrack(result.value(), track, caloLayer);
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void PHActsTrackProjection::updateSvtxTrack(
    const ActsPropagator::BoundTrackParamPair& parameters,
    SvtxTrack* svtxTrack,
    const int caloLayer)
{
  float pathlength = parameters.first / Acts::UnitConstants::cm;
  auto params = parameters.second;
  float calorad = m_caloRadii.find(m_caloTypes.at(caloLayer))->second;
  SvtxTrackState_v1 out(calorad);

  auto projectionPos = params.position(m_tGeometry->geometry().getGeoContext());
  const auto momentum = params.momentum();
  out.set_x(projectionPos.x() / Acts::UnitConstants::cm);
  out.set_y(projectionPos.y() / Acts::UnitConstants::cm);
  out.set_z(projectionPos.z() / Acts::UnitConstants::cm);
  out.set_px(momentum.x());
  out.set_py(momentum.y());
  out.set_pz(momentum.z());

  if (Verbosity() > 1)
  {
    std::cout << "Adding track state for caloLayer " << caloLayer
              << " at pathlength " << pathlength << " with position " << projectionPos.transpose() << std::endl;
  }

  ActsTransformations transformer;
  const auto globalCov = transformer.rotateActsCovToSvtxTrack(params);
  for (int i = 0; i < 6; ++i)
  {
    for (int j = 0; j < 6; ++j)
    {
      out.set_error(i, j, globalCov(i, j));
    }
  }

  svtxTrack->insert_state(&out);
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
  for (const auto& [key, cluster] : clusterMap)
  {
    const auto clusterEta =
        RawClusterUtility::GetPseudorapidity(*cluster,
                                             CLHEP::Hep3Vector(0, 0, 0));
    const auto dphi = deltaPhi(phi - cluster->get_phi());
    const auto deta = eta - clusterEta;
    const auto r = sqrt(pow(dphi, 2) + pow(deta, 2));

    if (r < minR)
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
  for (int iphi = phiBin - 2; iphi <= phiBin + 2; ++iphi)
  {
    for (int ieta = etaBin - 2; ieta <= etaBin + 2; ++ieta)
    {
      /// Check the phi periodic boundary conditions
      int wrapPhi = iphi;
      if (wrapPhi < 0)
        wrapPhi += m_towerGeomContainer->get_phibins();
      if (wrapPhi >= m_towerGeomContainer->get_phibins())
        wrapPhi -= m_towerGeomContainer->get_phibins();

      /// Check the eta boundary conditions
      if (ieta < 0 or ieta >= m_towerGeomContainer->get_etabins())
        continue;

      unsigned int towerkey = (ieta << 16U) + wrapPhi;
      auto tower = m_towerContainer->get_tower_at_key(towerkey);

      if (!tower)
        continue;

      energy5x5 += tower->get_energy();
      if (abs(iphi - phiBin) <= 1 and abs(ieta - etaBin) <= 1)
        energy3x3 += tower->get_energy();
    }
  }

  return;
}

PHActsTrackProjection::BoundTrackParamResult
PHActsTrackProjection::propagateTrack(
    const Acts::BoundTrackParameters& params,
    const int /*caloLayer*/,
    const SurfacePtr& targetSurf)
{
  ActsPropagator propagator(m_tGeometry);
  propagator.constField();
  propagator.verbosity(Verbosity());
  propagator.setConstFieldValue(1.4 * Acts::UnitConstants::T);

  return propagator.propagateTrackFast(params, targetSurf);
}

int PHActsTrackProjection::setCaloContainerNodes(PHCompositeNode* topNode,
                                                 const int caloLayer)
{
  std::string towerGeoNodeName = "TOWERGEOM_" + m_caloNames.at(caloLayer);
  std::string towerNodeName = "TOWERINFO_CALIB_" + m_caloNames.at(caloLayer);
  std::string clusterNodeName = "CLUSTER_" + m_caloNames.at(caloLayer);

  m_towerGeomContainer = findNode::getClass<RawTowerGeomContainer>(topNode, towerGeoNodeName.c_str());

  m_towerContainer = findNode::getClass<TowerInfoContainerv1>(topNode, towerNodeName.c_str());

  m_clusterContainer = findNode::getClass<RawClusterContainer>(topNode, clusterNodeName.c_str());

  if (m_useCemcPosRecalib and
      m_caloNames.at(caloLayer).compare("CEMC") == 0)
  {
    std::string nodeName = "CLUSTER_POS_COR_" + m_caloNames.at(caloLayer);
    m_clusterContainer = findNode::getClass<RawClusterContainer>(topNode, nodeName.c_str());
  }

  if (!m_towerGeomContainer or !m_towerContainer or !m_clusterContainer)
  {
    if (m_calosAvailable)
    {
      std::cout << PHWHERE
                << "Calo geometry and/or cluster container for " << m_caloNames.at(caloLayer)
                << "not found on node tree. Track projections to calos won't be filled."
                << std::endl;
    }

    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsTrackProjection::makeCaloSurfacePtrs(PHCompositeNode* topNode)
{
  for (int caloLayer = 0; caloLayer < m_nCaloLayers; caloLayer++)
  {
    if (setCaloContainerNodes(topNode, caloLayer) != Fun4AllReturnCodes::EVENT_OK)
    {
      continue;
    }

    /// Default to using calo radius
    double caloRadius = m_towerGeomContainer->get_radius();
    double caloOuterRadius = m_towerGeomContainer->get_radius() + m_towerGeomContainer->get_thickness();

    if (m_caloRadii.find(m_caloTypes.at(caloLayer)) != m_caloRadii.end())
    {
      caloRadius = m_caloRadii.find(m_caloTypes.at(caloLayer))->second;
    }
    else
    {
      m_caloRadii.insert(std::make_pair(m_caloTypes.at(caloLayer), caloRadius));
    }

    caloRadius *= Acts::UnitConstants::cm;

    if (m_caloRadii.find(m_caloTypes.at(caloLayer+m_nCaloLayers)) != m_caloRadii.end())
    {
      caloOuterRadius = m_caloRadii.find(m_caloTypes.at(caloLayer+m_nCaloLayers))->second;
    }
    else
    {
      m_caloRadii.insert(std::make_pair(m_caloTypes.at(caloLayer+m_nCaloLayers), caloOuterRadius));
    }

    caloOuterRadius *= Acts::UnitConstants::cm;

    /// Extend farther so that there is at least surface there, for high
    /// curling tracks. Can always reject later
    const auto eta = 2.5;
    const auto theta = 2. * atan(exp(-eta));
    const auto halfZ = caloRadius / tan(theta) * Acts::UnitConstants::cm;
    const auto halfZOuter = caloOuterRadius / tan(theta) * Acts::UnitConstants::cm;

    /// Make a cylindrical surface at (0,0,0) aligned along the z axis
    auto transform = Acts::Transform3::Identity();

    std::shared_ptr<Acts::CylinderSurface> surf =
        Acts::Surface::makeShared<Acts::CylinderSurface>(transform,
                                                         caloRadius,
                                                         halfZ);
    std::shared_ptr<Acts::CylinderSurface> outer_surf = 
        Acts::Surface::makeShared<Acts::CylinderSurface>(transform, 
                                                         caloOuterRadius, 
                                                         halfZOuter);
    if (Verbosity() > 1)
    {
      std::cout << "Creating  cylindrical surface at " << caloRadius << std::endl;
      std::cout << "Creating  cylindrical surface at " << caloOuterRadius << std::endl;
    }
    m_caloSurfaces.insert(std::make_pair(m_caloNames.at(caloLayer),surf));
    m_caloSurfaces.insert(std::make_pair(m_caloNames.at(caloLayer+m_nCaloLayers), outer_surf));
  }

  if (Verbosity() > 1)
  {
    for (const auto& [name, surfPtr] : m_caloSurfaces)
    {
      std::cout << "Cylinder " << name << " has center "
                << surfPtr.get()->center(m_tGeometry->geometry().getGeoContext()).transpose()
                << std::endl;
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsTrackProjection::getNodes(PHCompositeNode* topNode)
{
  m_vertexMap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
  if (!m_vertexMap)
  {
    std::cout << PHWHERE << "No vertex map on node tree, bailing."
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_tGeometry = findNode::getClass<ActsGeometry>(
      topNode, "ActsGeometry");
  if (!m_tGeometry)
  {
    std::cout << "ActsTrackingGeometry not on node tree. Exiting."
              << std::endl;

    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_trackMap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if (!m_trackMap)
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
    return phi + 2. * M_PI;
  else
    return phi;
}
