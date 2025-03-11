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
#include <globalvertex/SvtxVertex.h>
#include <globalvertex/SvtxVertexMap.h>

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

namespace
{
  static const std::map<SvtxTrack::CAL_LAYER,std::string> m_caloNames =
  {
    {SvtxTrack::CEMC, "CEMC"},
    {SvtxTrack::HCALIN, "HCALIN"},
    {SvtxTrack::HCALOUT, "HCALOUT"},
    {SvtxTrack::OUTER_CEMC, "OUTER_CEMC"},
    {SvtxTrack::OUTER_HCALIN, "OUTER_HCALIN"},
    {SvtxTrack::OUTER_HCALOUT, "OUTER_HCALOUT"}
  };
}

PHActsTrackProjection::PHActsTrackProjection(const std::string& name)
  : SubsysReco(name)
{}

int PHActsTrackProjection::InitRun(PHCompositeNode* topNode)
{
  if (Verbosity() > 1)
  {
    std::cout << "PHActsTrackProjection begin Init" << std::endl;
  }

  int ret = makeCaloSurfacePtrs(topNode);

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

int PHActsTrackProjection::process_event(PHCompositeNode* /*topNode*/)
{

  for( const auto& [layer, name]:m_caloNames )
  {
    if (Verbosity())
    {
      std::cout << "Processing calo layer " << name << std::endl;
    }
    int ret = projectTracks(layer);
    if (ret != Fun4AllReturnCodes::EVENT_OK)
    {
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

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

int PHActsTrackProjection::projectTracks(SvtxTrack::CAL_LAYER caloLayer)
{

  // make sure caloSurface is valid
  const auto surface_iter = m_caloSurfaces.find(caloLayer);
  if( surface_iter == m_caloSurfaces.end() ) return Fun4AllReturnCodes::EVENT_OK;
  const auto& cylSurf = surface_iter->second;

  // create propagator
  ActsPropagator prop(m_tGeometry);

  // loop over tracks
  for (const auto& [key, track] : *m_trackMap)
  {
    auto params = prop.makeTrackParams(track, m_vertexMap);
    if(!params.ok())
    {
      continue;
    }

    // propagate
    const auto result = propagateTrack(params.value(), cylSurf);
    if (result.ok())
    {
      // update track
      updateSvtxTrack(result.value(), track, caloLayer);
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void PHActsTrackProjection::updateSvtxTrack(
    const ActsPropagator::BoundTrackParamPair& parameters,
    SvtxTrack* svtxTrack,
    SvtxTrack::CAL_LAYER caloLayer)
{
  const float pathlength = parameters.first / Acts::UnitConstants::cm;
  const auto params = parameters.second;
  const float calorad = m_caloRadii.at(caloLayer);
  SvtxTrackState_v1 out(calorad);

  const auto projectionPos = params.position(m_tGeometry->geometry().getGeoContext());
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

PHActsTrackProjection::BoundTrackParamResult
PHActsTrackProjection::propagateTrack(
    const Acts::BoundTrackParameters& params,
    const SurfacePtr& targetSurf)
{
  ActsPropagator propagator(m_tGeometry);
  propagator.constField();
  propagator.verbosity(Verbosity());
  propagator.setConstFieldValue(m_constFieldVal * Acts::UnitConstants::T);

  return propagator.propagateTrackFast(params, targetSurf);
}

int PHActsTrackProjection::makeCaloSurfacePtrs(PHCompositeNode* topNode)
{
  using calo_pair_t = std::pair<SvtxTrack::CAL_LAYER,SvtxTrack::CAL_LAYER>;

  for( const auto& [first, second]:std::initializer_list<calo_pair_t>{
    {SvtxTrack::CEMC,SvtxTrack::OUTER_CEMC},
    {SvtxTrack::HCALIN,SvtxTrack::OUTER_HCALIN },
    {SvtxTrack::HCALOUT,SvtxTrack::OUTER_HCALOUT }} )
  {
    const auto& caloname( m_caloNames.at(first) );

    // get tower geometry node name
    const std::string towerGeoNodeName = "TOWERGEOM_" + caloname;

    // get tower geometry container
    const auto towerGeomContainer = findNode::getClass<RawTowerGeomContainer>(topNode, towerGeoNodeName.c_str());

    if( !towerGeomContainer )
    {
      std::cout << PHWHERE << "-"
        << " Calo tower geometry container for " << caloname
        << " not found on node tree. Track projections to calos won't be filled."
        << std::endl;
      continue;
    }

    // get calorimeter inner radius and store
    double caloRadius = towerGeomContainer->get_radius();
    {
      const auto iter = m_caloRadii.find(first);
      if( iter == m_caloRadii.end() ) m_caloRadii.emplace(first,caloRadius);
      else caloRadius = iter->second;
    }

    // get calorimeter outer radius and store
    double caloOuterRadius = towerGeomContainer->get_radius() + towerGeomContainer->get_thickness();
    {
      const auto iter = m_caloRadii.find(second);
      if( iter == m_caloRadii.end() ) m_caloRadii.emplace(second,caloOuterRadius);
      else caloOuterRadius = iter->second;
    }

    // convert to ACTS units
    caloRadius *= Acts::UnitConstants::cm;
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
    m_caloSurfaces.emplace(first,surf);
    m_caloSurfaces.emplace(second,outer_surf);
  }

  if (Verbosity() > 1)
  {
    for (const auto& [layer, surfPtr] : m_caloSurfaces)
    {
      std::cout << "Cylinder " << m_caloNames.at(layer) << " has center "
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
