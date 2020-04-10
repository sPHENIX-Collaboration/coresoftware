/// See header file for comments on class, functions, member variables
#include "PHActsSourceLinks.h"
#include "MakeActsGeometry.h"

/// Tracking includes
#include <g4detectors/PHG4CylinderCellGeom.h>
#include <g4detectors/PHG4CylinderCellGeomContainer.h>
#include <g4detectors/PHG4CylinderGeom.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <intt/CylinderGeomIntt.h>
#include <intt/InttDefs.h>
#include <mvtx/CylinderGeom_Mvtx.h>
#include <mvtx/MvtxDefs.h>
#include <tpc/TpcDefs.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>

/// Fun4All includes
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <phool/phool.h>

/// Acts includes
#include <ACTFW/Utilities/Options.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>
#include <Acts/Surfaces/PlaneSurface.hpp>
#include <Acts/Surfaces/Surface.hpp>

/// std (and the like) includes
#include <cmath>
#include <iostream>
#include <memory>
#include <utility>

/// Root includes
#include <TGeoNode.h>
#include <TMatrixDSym.h>
#include <TMatrixT.h>
#include <TVector3.h>

PHActsSourceLinks::PHActsSourceLinks(const std::string &name)
  : SubsysReco(name)
  , m_clusterMap(nullptr)
  , m_hitIdClusKey(nullptr)
  , m_sourceLinks(nullptr)
  , m_actsGeometry(nullptr)
  , m_geomContainerMvtx(nullptr)
  , m_geomContainerIntt(nullptr)
  , m_geomContainerTpc(nullptr)
  , m_minSurfZ(0.0)  /// Initialize these to 0, get from MakeActsGeometry
  , m_maxSurfZ(0.0)
  , m_nSurfZ(0)
  , m_nSurfPhi(0)
  , m_surfStepPhi(0.0)
  , m_surfStepZ(0.0)
  , m_moduleStepPhi(0.0)
  , m_modulePhiStart(0.0)
{
  Verbosity(0);
}

int PHActsSourceLinks::Init(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsSourceLinks::InitRun(PHCompositeNode *topNode)
{
  if (Verbosity() > 10)
  {
    std::cout << "Starting PHActsSourceLinks::InitRun" << std::endl;
  }

  /// Check and create nodes that this module will build
  createNodes(topNode);

  /// Check if Acts geometry has been built and is on the node tree
  MakeActsGeometry *actsGeometry = new MakeActsGeometry();
  actsGeometry->BuildAllGeometry(topNode);

  m_minSurfZ = actsGeometry->getMinSurfZ();
  m_maxSurfZ = actsGeometry->getMaxSurfZ();
  m_nSurfZ = actsGeometry->getNSurfZ();
  m_nSurfPhi = actsGeometry->getNSurfPhi();

  /// These are arbitrary TPC Surface subdivisions, set in MakeActsGeometry
  m_surfStepZ = actsGeometry->getSurfStepZ();
  m_moduleStepPhi = actsGeometry->getModuleStepPhi();
  m_modulePhiStart = actsGeometry->getModulePhiStart();
  m_surfStepPhi = actsGeometry->getSurfStepPhi();

  m_clusterNodeMap = actsGeometry->getNodeMap();
  m_clusterSurfaceMap = actsGeometry->getSurfaceMapSilicon();
  m_clusterSurfaceMapTpc = actsGeometry->getSurfaceMapTpc();
  m_actsGeometry->tGeometry = actsGeometry->getTGeometry();
  m_actsGeometry->magField = actsGeometry->getMagField();
  m_actsGeometry->geoContext = actsGeometry->getGeoContext();
  m_actsGeometry->magFieldContext = actsGeometry->getMagFieldContext();
  m_actsGeometry->calibContext = actsGeometry->getCalibContext();

  if (Verbosity() > 10)
  {
    std::cout << "Finished PHActsSourceLinks::InitRun" << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsSourceLinks::process_event(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
  {
    std::cout << "Start PHActsSourceLinks process_event" << std::endl;
  }

  /// Get the nodes from the node tree
  if (getNodes(topNode) == Fun4AllReturnCodes::ABORTEVENT)
    return Fun4AllReturnCodes::ABORTEVENT;

  /// Arbitrary hitId that is used to mape between cluster key and an
  /// unsigned int which Acts can take
  unsigned int hitId = 0;

  TrkrClusterContainer::ConstRange clusRange = m_clusterMap->getClusters();
  TrkrClusterContainer::ConstIterator clusIter;

  for (clusIter = clusRange.first; clusIter != clusRange.second; ++clusIter)
  {
    const TrkrDefs::cluskey clusKey = clusIter->first;
    const TrkrCluster *cluster = clusIter->second;

    const unsigned int layer = TrkrDefs::getLayer(clusKey);

    /// Create the clusKey hitId pair to insert into the map
    const unsigned int trkrId = TrkrDefs::getTrkrId(clusKey);
    m_hitIdClusKey->insert(std::pair<TrkrDefs::cluskey, unsigned int>(clusKey, hitId));

    /// Local coordinates and surface to be set by the correct tracking
    /// detector function below
    TMatrixD localErr(3, 3);
    double local2D[2] = {0};
    Surface surface;

    /// Run the detector specific function for getting the local coordinates
    /// of the cluster, as well as the corresponding Acts::Surface
    if (trkrId == TrkrDefs::mvtxId)
    {
      surface = getMvtxLocalCoords(local2D, localErr, cluster, clusKey);
      /// Make sure things returned appropriately
      if (!surface)
      {
        /// Error message already printed in function, return abort
        return Fun4AllReturnCodes::ABORTEVENT;
      }
    }
    else if (trkrId == TrkrDefs::inttId)
    {
      surface = getInttLocalCoords(local2D, localErr, cluster, clusKey);
      if (!surface)
      {
        /// Error message already printed in function, return abort
        return Fun4AllReturnCodes::ABORTEVENT;
      }
    }
    else if (trkrId == TrkrDefs::tpcId)
    {
      surface = getTpcLocalCoords(local2D, localErr, cluster, clusKey);

      if (!surface)
      {
        /// Error message already printed in function, return abort
        return Fun4AllReturnCodes::ABORTEVENT;
      }
    }
    else
    {
      std::cout << "Invalid trkrId found in " << PHWHERE
                << std::endl
                << "Skipping this cluster"
                << std::endl;
    }

    /// ====================================================
    /// Finished with detector specific cluster stuff
    /// We have the data needed to construct an Acts  measurement
    /// for this cluster
    /// ====================================================

    /// Get the 2D location covariance uncertainty for the cluster
    Acts::BoundMatrix cov = Acts::BoundMatrix::Zero();
    cov(Acts::eLOC_0, Acts::eLOC_0) = localErr[0][0];
    cov(Acts::eLOC_1, Acts::eLOC_0) = localErr[1][0];
    cov(Acts::eLOC_0, Acts::eLOC_1) = localErr[0][1];
    cov(Acts::eLOC_1, Acts::eLOC_1) = localErr[1][1];

    /// local and localErr contain the position and covariance
    /// matrix in local coords
    if (Verbosity() > 0)
    {
      std::cout << "    layer " << layer << std::endl;
      for (int i = 0; i < 2; ++i)
      {
        std::cout << " i " << i << "   local 2D position " << local2D[i]
                  << std::endl;
      }

      std::cout << "    local covariance matrix:" << std::endl;
      std::cout << cov << std::endl;
    }

    /// Cluster positions on GeoObject/Surface
    Acts::BoundVector loc = Acts::BoundVector::Zero();
    loc[Acts::eLOC_0] = local2D[0];
    loc[Acts::eLOC_1] = local2D[1];

    if (Verbosity() > 0)
    {
      std::cout << "Layer " << layer
                << " create measurement for trkrid " << trkrId
                << " surface " << surface->name() << " surface type "
                << surface->type() << " local x " << loc[Acts::eLOC_0]
                << " local y " << loc[Acts::eLOC_1] << std::endl;
    }

    /// TrkrClusterSourceLink creates an Acts::FittableMeasurement
    SourceLink sourceLink(hitId, surface, loc, cov);

    /// Add the sourceLink to the container
    m_sourceLinks->insert(std::pair<unsigned int, SourceLink>(hitId, sourceLink));

    hitId++;
  }

  if (Verbosity() > 10)
  {
    //m_hitIdClusKey
    std::map<TrkrDefs::cluskey, unsigned int>::iterator it = m_hitIdClusKey->begin();
    while (it != m_hitIdClusKey->end())
    {
      std::cout << "cluskey " << it->first << " has hitid " << it->second
                << std::endl;
      ++it;
    }
  }

  /// Add the data to the nodes that were previously created (?)

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsSourceLinks::End(PHCompositeNode *topNode)
{
  if (Verbosity() > 10)
  {
    std::cout << "Finished PHActsSourceLinks" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

Surface PHActsSourceLinks::getTpcLocalCoords(double (&local2D)[2],
                                             TMatrixD &localErr,
                                             const TrkrCluster *cluster,
                                             const TrkrDefs::cluskey clusKey)
{
  const float x = cluster->getPosition(0);
  const float y = cluster->getPosition(1);
  const float z = cluster->getPosition(2);

  // In local coords the covariances are in the  r*phi vs z frame
  // They have been rotated into global coordinates in TrkrCluster
  TMatrixD worldErr(3, 3);
  for (int i = 0; i < 3; ++i)
  {
    for (int j = 0; j < 3; j++)
    {
      worldErr[i][j] = cluster->getError(i, j);
    }
  }

  /// Extract detector element IDs to access the correct Surface
  TVector3 world(x, y, z);

  /// Get some geometry values
  const double clusPhi = atan2(world[1], world[0]);
  const double radius = sqrt(x * x + y * y);
  const double rClusPhi = radius * clusPhi;
  const double zTpc = world[2];

  const unsigned int layer = TrkrDefs::getLayer(clusKey);
  const unsigned int sectorId = TpcDefs::getSectorId(clusKey);
  const unsigned int side = TpcDefs::getSide(clusKey);

  const double modulePhiLow = m_modulePhiStart + (double) sectorId * m_moduleStepPhi;

  const unsigned int iPhi = (clusPhi - modulePhiLow) / m_surfStepPhi;
  const unsigned int iZ = fabs(zTpc) / m_surfStepZ;
  const unsigned int iPhiZ = iPhi + 100. * iZ;  /// for making a map key

  if (Verbosity() > 0)
  {
    const double checkSurfRphiCenter = radius * (modulePhiLow + (double) iPhi * m_surfStepPhi + m_surfStepPhi / 2.0);
    double checkSurfZCenter = (double) iZ * m_surfStepZ + m_surfStepZ / 2.0;
    if (side == 0)
      checkSurfZCenter = -checkSurfZCenter;

    std::cout << " xworld " << world[0] << " yworld "
              << world[1] << " clusphi " << clusPhi << std::endl;
    std::cout << " sectorid " << sectorId << " side "
              << side << " iphi " << iPhi << " iz "
              << iZ << " i_phi_z " << iPhiZ << std::endl;
    std::cout << "r_clusphi " << rClusPhi
              << " check_surf_rphi center " << checkSurfRphiCenter
              << " ztpc " << zTpc << " check_surf_z_center "
              << checkSurfZCenter << std::endl;
  }

  /// Get the surface key to find the surface from the map
  const TrkrDefs::cluskey surfkey = TpcDefs::genClusKey(layer, sectorId,
                                                        side, iPhiZ);
  std::map<TrkrDefs::cluskey, Surface>::iterator surfIter;

  surfIter = m_clusterSurfaceMapTpc->find(surfkey);
  if (surfIter == m_clusterSurfaceMapTpc->end())
  {
    std::cout << PHWHERE << "Failed to find surface, should be impossible!" << std::endl;
    return nullptr;
  }

  Surface surface = surfIter->second->getSharedPtr();

  /// Transformation of cluster to local surface coords
  /// Coords are r*phi relative to surface r-phi center, and z
  /// relative to surface z center
  Acts::Vector3D center = surface->center(m_actsGeometry->geoContext);
  double surfRphiCenter = atan2(center[1], center[0]) * radius;
  double surfZCenter = center[2];
  if (Verbosity() > 0)
  {
    std::cout << "surface center readback:   x " << center[0]
              << " y " << center[1] << " phi "
              << atan2(center[1], center[0]) << " z " << center[2]
              << " surf_rphi_center " << surfRphiCenter
              << " surf_z_center " << surfZCenter
              << std::endl;
  }

  local2D[0] = rClusPhi - surfRphiCenter;
  local2D[1] = zTpc - surfZCenter;

  if (Verbosity() > 0)
  {
    std::cout << "r_clusphi " << rClusPhi << " surf_rphi center "
              << surfRphiCenter
              << " ztpc " << zTpc << " surf_z_center " << surfZCenter
              << " local rphi " << local2D[0] << " local z " << local2D[1]
              << std::endl;
  }

  localErr = transformCovarToLocal(clusPhi, worldErr);

  return surface;
}

Surface PHActsSourceLinks::getInttLocalCoords(double (&local2D)[2],
                                              TMatrixD &localErr,
                                              const TrkrCluster *cluster,
                                              const TrkrDefs::cluskey clusKey)
{
  TVector3 local(0, 0, 0);

  const float x = cluster->getPosition(0);
  const float y = cluster->getPosition(1);
  const float z = cluster->getPosition(2);

  // In local coords the covariances are in the  r*phi vs z frame
  // They have been rotated into global coordinates in TrkrCluster
  TMatrixD worldErr(3, 3);
  for (int i = 0; i < 3; ++i)
  {
    for (int j = 0; j < 3; j++)
    {
      worldErr[i][j] = cluster->getError(i, j);
    }
  }

  /// Extract detector element IDs to access the correct Surface
  const TVector3 world(x, y, z);

  /// Get the INTT geometry
  const unsigned int ladderZId = InttDefs::getLadderZId(clusKey);
  const unsigned int ladderPhiId = InttDefs::getLadderPhiId(clusKey);
  const unsigned int layer = TrkrDefs::getLayer(clusKey);

  const TrkrDefs::hitsetkey hitSetKey = InttDefs::genHitSetKey(layer, ladderZId, ladderPhiId);
  if (Verbosity() > 10)
  {
    std::cout << "Intt cluster with ladderzid " << ladderZId
              << " ladderphid " << ladderPhiId << std::endl;
  }
  /// Get the TGeoNode
  const TGeoNode *sensorNode = getNodeFromClusterMap(hitSetKey);

  if (!sensorNode)
  {
    std::cout << PHWHERE << "No entry in TGeo map for cluster: layer "
              << layer << " ladderZId " << ladderZId
              << " ladderPhiId " << ladderPhiId
              << " - should be impossible!" << std::endl;
    return nullptr;
  }

  /// Now we have the geo node, so find the corresponding Acts::Surface
  Surface surface = getSurfaceFromClusterMap(hitSetKey);
  if (!surface)
  {
    std::cout << PHWHERE
              << "Failed to find associated surface element - should be impossible!"
              << std::endl;
    return nullptr;
  }

  // transform position back to local coords on sensor
  CylinderGeomIntt *layerGeom =
      dynamic_cast<CylinderGeomIntt *>(m_geomContainerIntt->GetLayerGeom(layer));
  local = layerGeom->get_local_from_world_coords(ladderZId,
                                                 ladderPhiId,
                                                 world);
  local2D[0] = local[1];  // r*phi
  local2D[1] = local[2];  // z

  if (Verbosity() > 10)
  {
    double segcent[3];
    layerGeom->find_segment_center(ladderZId, ladderPhiId, segcent);
    std::cout << "   segment center: " << segcent[0]
              << " " << segcent[1] << " " << segcent[2] << std::endl;
    std::cout << "   world; " << world[0] << " " << world[1]
              << " " << world[2] << std::endl;
    std::cout << "   local; " << local[0] << " " << local[1]
              << " " << local[2] << std::endl;
  }

  /// Get the local covariance error
  localErr = getInttCovarLocal(layer, ladderZId, ladderPhiId, worldErr);

  return surface;
}

Surface PHActsSourceLinks::getMvtxLocalCoords(double (&local2D)[2],
                                              TMatrixD &localErr,
                                              const TrkrCluster *cluster,
                                              const TrkrDefs::cluskey clusKey)
{
  TVector3 local(0, 0, 0);

  const float x = cluster->getPosition(0);
  const float y = cluster->getPosition(1);
  const float z = cluster->getPosition(2);

  // In local coords the covariances are in the  r*phi vs z frame
  // They have been rotated into global coordinates in TrkrCluster
  TMatrixD worldErr(3, 3);
  for (int i = 0; i < 3; ++i)
  {
    for (int j = 0; j < 3; j++)
    {
      worldErr[i][j] = cluster->getError(i, j);
    }
  }

  /// Extract detector element IDs to access the correct Surface
  const TVector3 world(x, y, z);

  /// Get the Mvtx geometry
  const unsigned int staveId = MvtxDefs::getStaveId(clusKey);
  const unsigned int chipId = MvtxDefs::getChipId(clusKey);
  const unsigned int layer = TrkrDefs::getLayer(clusKey);

  if (Verbosity() > 10)
  {
    std::cout << "MVTX cluster with stave id: " << staveId
              << " and chip Id " << chipId << std::endl;
  }

  /// Generate the hitsetkey
  TrkrDefs::hitsetkey hitSetKey = MvtxDefs::genHitSetKey(layer,
                                                         staveId,
                                                         chipId);

  /// Get the TGeoNode
  const TGeoNode *sensorNode = getNodeFromClusterMap(hitSetKey);

  if (!sensorNode)
  {
    std::cout << PHWHERE << "No entry in TGeo map for cluster: layer "
              << layer << " staveid " << staveId
              << " chipid " << chipId
              << " - should be impossible!" << std::endl;
    return nullptr;
  }

  /// Now we have the geo node, so find the corresponding Acts::Surface
  Surface surface = getSurfaceFromClusterMap(hitSetKey);
  if (!surface)
  {
    std::cout << PHWHERE
              << "Failed to find associated surface element - should be impossible!"
              << std::endl;
    return nullptr;
  }

  CylinderGeom_Mvtx *layerGeom = dynamic_cast<CylinderGeom_Mvtx *>(m_geomContainerMvtx->GetLayerGeom(layer));
  local = layerGeom->get_local_from_world_coords(staveId, 0, 0,
                                                 chipId, world);
  local2D[0] = local[0];
  local2D[1] = local[2];

  if (Verbosity() > 10)
  {
    double segcent[3];
    layerGeom->find_sensor_center(staveId, 0, 0, chipId, segcent);
    std::cout << "   segment center: " << segcent[0] << " "
              << segcent[1] << " " << segcent[2] << std::endl;
    std::cout << "   world; " << world[0] << " "
              << world[1] << " " << world[2] << std::endl;
    std::cout << "   local; " << local[0] << " "
              << local[1] << " " << local[2] << std::endl;
  }

  // transform covariance matrix back to local coords on chip
  localErr = getMvtxCovarLocal(layer, staveId, chipId, worldErr);

  return surface;
}

int PHActsSourceLinks::getNodes(PHCompositeNode *topNode)
{
  m_clusterMap = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!m_clusterMap)
  {
    std::cout << PHWHERE
              << "TRKR_CLUSTER node not found on node tree. Exiting."
              << std::endl;

    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_geomContainerMvtx = findNode::getClass<
      PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_MVTX");
  if (!m_geomContainerMvtx)
  {
    std::cout << PHWHERE << " CYLINDERGEOM_MVTX  node not found on node tree"
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_geomContainerTpc =
      findNode::getClass<PHG4CylinderCellGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  if (!m_geomContainerTpc)
  {
    std::cout << PHWHERE << "ERROR: Can't find node CYLINDERCELLGEOM_SVTX"
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_geomContainerIntt = findNode::getClass<
      PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_INTT");
  if (!m_geomContainerIntt)
  {
    std::cout << PHWHERE << " CYLINDERGEOM_INTT  node not found on node tree"
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void PHActsSourceLinks::createNodes(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  /// Get the DST Node
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));

  /// Check that it is there
  if (!dstNode)
  {
    std::cerr << "DST Node missing, quitting" << std::endl;
    throw std::runtime_error("failed to find DST node in PHActsSourceLinks::createNodes");
  }

  /// Get the tracking subnode
  PHCompositeNode *svtxNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "SVTX"));

  /// Check that it is there
  if (!svtxNode)
  {
    svtxNode = new PHCompositeNode("SVTX");
    dstNode->addNode(svtxNode);
  }

  /// See if the map is already on the node tree
  m_hitIdClusKey = findNode::getClass<std::map<TrkrDefs::cluskey, unsigned int>>(topNode, "HitIDClusIDActsMap");

  /// If not add it
  if (!m_hitIdClusKey)
  {
    m_hitIdClusKey = new std::map<TrkrDefs::cluskey, unsigned int>;
    PHDataNode<std::map<TrkrDefs::cluskey, unsigned int>> *hitMapNode =
        new PHDataNode<std::map<TrkrDefs::cluskey, unsigned int>>(m_hitIdClusKey, "HitIDClusIDActsMap");
    svtxNode->addNode(hitMapNode);
  }

  m_clusterSurfaceMap = findNode::getClass<std::map<TrkrDefs::hitsetkey, Surface>>(topNode, "ClusterSurfaceActsMap");
  
  if(!m_clusterSurfaceMap)
    {
      m_clusterSurfaceMap = new std::map<TrkrDefs::hitsetkey, Surface>;
      PHDataNode<std::map<TrkrDefs::hitsetkey, Surface>> *keySurfaceNode = 
	new PHDataNode<std::map<TrkrDefs::hitsetkey, Surface>>(m_clusterSurfaceMap, "HitSetKeySurfaceActsMap");
      svtxNode->addNode(keySurfaceNode);
    }

    m_clusterSurfaceMapTpc = findNode::getClass<std::map<TrkrDefs::cluskey, Surface>>(topNode, "ClusterTpcSurfaceActsMap");
  
  if(!m_clusterSurfaceMapTpc)
    {
      m_clusterSurfaceMapTpc = new std::map<TrkrDefs::cluskey, Surface>;
      PHDataNode<std::map<TrkrDefs::cluskey, Surface>> *keyTpcSurfaceNode = 
	new PHDataNode<std::map<TrkrDefs::cluskey, Surface>>(m_clusterSurfaceMapTpc, "ClusKeyTpcSurfaceActsMap");
      svtxNode->addNode(keyTpcSurfaceNode);
    }

  m_actsGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsFitCfg");

  if (!m_actsGeometry)
  {
    m_actsGeometry = new ActsGeometry();
    PHDataNode<ActsGeometry> *fitNode = new PHDataNode<ActsGeometry>(m_actsGeometry, "ActsGeometry");
    svtxNode->addNode(fitNode);
  }

  /// Do the same for the SourceLink container
  m_sourceLinks = findNode::getClass<std::map<unsigned int, SourceLink>>(topNode, "TrkrClusterSourceLinks");

  if (!m_sourceLinks)
  {
    m_sourceLinks = new std::map<unsigned int, SourceLink>;
    PHDataNode<std::map<unsigned int, SourceLink>>
        *sourceLinkNode = new PHDataNode<std::map<unsigned int, SourceLink>>(m_sourceLinks, "TrkrClusterSourceLinks");

    svtxNode->addNode(sourceLinkNode);
  }

  return;
}

TGeoNode *PHActsSourceLinks::getNodeFromClusterMap(TrkrDefs::hitsetkey hitSetKey)
{
  /// Get the TGeoNode for this hit set key
  std::map<TrkrDefs::hitsetkey, TGeoNode *>::iterator mapIter;
  mapIter = m_clusterNodeMap.find(hitSetKey);
  TGeoNode *sensorNode;

  /// Make sure we found it
  if (mapIter != m_clusterNodeMap.end())
  {
    sensorNode = mapIter->second;
    if (Verbosity() > 0)
    {
      std::cout << "Found TGeoNode in m_clusterNodeMap for hitsetkey "
                << hitSetKey
                << " node " << sensorNode->GetName()
                << std::endl;
    }
  }
  else
  {
    /// If it couldn't be found then return nullptr
    return nullptr;
  }

  // Otherwise return node
  return sensorNode;
}

Surface PHActsSourceLinks::getSurfaceFromClusterMap(TrkrDefs::hitsetkey hitSetKey)
{
  Surface surface;

  std::map<TrkrDefs::hitsetkey, Surface>::iterator
      surfaceIter;

  surfaceIter = m_clusterSurfaceMap->find(hitSetKey);

  /// Check to make sure we found the surface in the map
  if (surfaceIter != m_clusterSurfaceMap->end())
  {
    surface = surfaceIter->second;
    if (Verbosity() > 0)
    {
      std::cout << "Got surface pair " << surface->name()
                << " surface type " << surface->type()
                << std::endl;
    }
  }
  else
  {
    /// If it doesn't exit, return nullptr
    return nullptr;
  }

  return surface;
}

TMatrixD PHActsSourceLinks::getMvtxCovarLocal(const unsigned int layer, const unsigned int staveId, const unsigned int chipId, TMatrixD worldErr)
{
  TMatrixD localErr(3, 3);

  // rotate errors back to local coords
  double ladderLocation[3] = {0.0, 0.0, 0.0};

  // returns the center of the sensor in world coordinates - used to get the ladder phi location
  CylinderGeom_Mvtx *layergeom = dynamic_cast<CylinderGeom_Mvtx *>(m_geomContainerMvtx->GetLayerGeom(layer));
  layergeom->find_sensor_center(staveId, 0, 0, chipId, ladderLocation);
  double ladderPhi = atan2(ladderLocation[1], ladderLocation[0]);
  ladderPhi += layergeom->get_stave_phi_tilt();

  localErr = transformCovarToLocal(ladderPhi, worldErr);

  if (Verbosity() > 10)
  {
    for (int i = 0; i < 3; ++i)
    {
      for (int j = 0; j < 3; ++j)
      {
        std::cout << "  " << i << "    " << j
                  << " local_err " << localErr[i][j]
                  << std::endl;
      }
    }
  }

  return localErr;
}

TMatrixD PHActsSourceLinks::getInttCovarLocal(const unsigned int layer, const unsigned int ladderZId, const unsigned int ladderPhiId, TMatrixD worldErr)
{
  TMatrixD localErr(3, 3);

  // rotate errors back to local coords
  double ladderLocation[3] = {0.0, 0.0, 0.0};

  // rotate errors back to local coords
  CylinderGeomIntt *layerGeom = dynamic_cast<CylinderGeomIntt *>(m_geomContainerIntt->GetLayerGeom(layer));
  layerGeom->find_segment_center(ladderZId, ladderPhiId, ladderLocation);
  double ladderPhi = atan2(ladderLocation[1], ladderLocation[0]);

  localErr = transformCovarToLocal(ladderPhi, worldErr);

  if (Verbosity() > 10)
  {
    for (int i = 0; i < 3; ++i)
    {
      for (int j = 0; j < 3; ++j)
      {
        std::cout << "  INTT: " << i << "    " << j
                  << " local_err " << localErr[i][j] << std::endl;
      }
    }
  }
  return localErr;
}

TMatrixD PHActsSourceLinks::transformCovarToLocal(const double ladderPhi,
                                                  TMatrixD worldErr)
{
  TMatrixD localErr(3, 3);

  // this is the matrix that was used to rotate from local to global coords
  TMatrixD ROT(3, 3);
  ROT[0][0] = cos(ladderPhi);
  ROT[0][1] = -1.0 * sin(ladderPhi);
  ROT[0][2] = 0.0;
  ROT[1][0] = sin(ladderPhi);
  ROT[1][1] = cos(ladderPhi);
  ROT[1][2] = 0.0;
  ROT[2][0] = 0.0;
  ROT[2][1] = 0.0;
  ROT[2][2] = 1.0;
  // we want the inverse rotation
  ROT.Invert();

  TMatrixD ROT_T(3, 3);
  ROT_T.Transpose(ROT);

  localErr = ROT * worldErr * ROT_T;

  return localErr;
}
