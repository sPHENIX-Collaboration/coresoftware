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
#include <trackbase_historic/SvtxVertexMap.h>
#include <trackbase_historic/SvtxVertex.h>

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
#include <Acts/Surfaces/PerigeeSurface.hpp>
#include <Acts/Surfaces/PlaneSurface.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/Utilities/Units.hpp>

#include <ActsExamples/Utilities/Options.hpp>

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
  , m_useVertexMeasurement(false)
  , m_clusterMap(nullptr)
  , m_actsGeometry(nullptr)
  , m_hitIdClusKey(nullptr)
  , m_sourceLinks(nullptr)
  , m_magField("1.4")
  , m_magFieldRescale(-1.0)
  , m_geomContainerMvtx(nullptr)
  , m_geomContainerIntt(nullptr)
  , m_geomContainerTpc(nullptr)
  , m_tGeometry(nullptr)
 
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
  m_actsGeometry = new MakeActsGeometry();
  
  m_actsGeometry->setVerbosity(Verbosity());
  m_actsGeometry->setMagField(m_magField);
  m_actsGeometry->setMagFieldRescale(m_magFieldRescale);
  m_actsGeometry->buildAllGeometry(topNode);

  /// Set the tGeometry struct to be put on the node tree
  m_tGeometry->tGeometry = m_actsGeometry->getTGeometry();
  m_tGeometry->magField = m_actsGeometry->getMagField();
  m_tGeometry->calibContext = m_actsGeometry->getCalibContext();
  m_tGeometry->magFieldContext = m_actsGeometry->getMagFieldContext();
  m_tGeometry->geoContext = m_actsGeometry->getGeoContext();

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

if(m_useVertexMeasurement){
  addVerticesAsSourceLinks(topNode, hitId);

 }
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
    Acts::Vector2D local2D(0,0);
    
    Surface surface;
    Acts::BoundMatrix cov = Acts::BoundMatrix::Zero();
    
    /// Run the detector specific function for getting the local coordinates
    /// of the cluster, as well as the corresponding Acts::Surface
    if (trkrId == TrkrDefs::mvtxId)
    {
      surface = getMvtxLocalCoords(local2D, cov, cluster, clusKey);
      /// Make sure things returned appropriately
      if (!surface)
      {
        /// if we couldn't find the surface (shouldn't happen) just skip this hit
	continue;
      }
    }
    else if (trkrId == TrkrDefs::inttId)
    {
      surface = getInttLocalCoords(local2D, cov, cluster, clusKey);
      if (!surface)
      {
	/// if we couldn't find the surface (shouldn't happen) just skip this hit
	continue;
      }
    }
    else if (trkrId == TrkrDefs::tpcId)
    {
      surface = getTpcLocalCoords(local2D, cov, cluster, clusKey);

      if (!surface)
      {
        /// if we couldn't find the surface (shouldn't happen) just skip this hit
	continue;
      }
    }
    else if (trkrId == TrkrDefs::micromegasId)
      {
	// skip micromegas for now
	continue;
      }
    else
    {
      std::cout << "Invalid trkrId found in " << PHWHERE
                << std::endl
                << "Skipping this cluster"
                << std::endl;
	continue;
    }

    /// ====================================================
    /// Finished with detector specific cluster stuff
    /// We have the data needed to construct an Acts  measurement
    /// for this cluster
    /// ====================================================


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
                << " local y " << loc[Acts::eLOC_1] << std::endl << std::endl;
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


  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsSourceLinks::ResetEvent(PHCompositeNode *topNode)
{
  /// Clear out the maps for the next event
  m_hitIdClusKey->clear();
  m_sourceLinks->clear();
  

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

Surface PHActsSourceLinks::getTpcLocalCoords(Acts::Vector2D &local2D,
                                             Acts::BoundMatrix &localErr,
                                             const TrkrCluster *cluster,
                                             const TrkrDefs::cluskey clusKey)
{
  // cm
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
      worldErr[i][j] = cluster->getError(i, j);  // this is the cov error squared
    }
  }

  /// Extract detector element IDs to access the correct Surface
  TVector3 world(x, y, z);

  /// Get some geometry values (lengths in mm)
  const double clusPhi = atan2(world[1], world[0]);
  const double radius = sqrt(x * x + y * y) * Acts::UnitConstants::cm;
  const double rClusPhi = radius * clusPhi;
  const double zTpc = world[2] * Acts::UnitConstants::cm;

  const unsigned int layer = TrkrDefs::getLayer(clusKey);
  const unsigned int sectorId = TpcDefs::getSectorId(clusKey);
  const unsigned int side = TpcDefs::getSide(clusKey);
 
  /// Get the surface key to find the surface from the map
  TrkrDefs::hitsetkey tpcHitSetKey = TpcDefs::genHitSetKey(layer, sectorId, side);
  std::vector<double> worldVec = {world[0], world[1], world[2]};
  /// MakeActsGeometry has a helper function since many surfaces can exist on
  /// a given readout module
  Surface surface = m_actsGeometry->getTpcSurfaceFromCoords(tpcHitSetKey,
							    worldVec);

  /// If surface can't be found (shouldn't happen) return nullptr and skip this cluster
  if(!surface)
    {
      std::cout << PHWHERE
		<< "Failed to find associated surface element - should be impossible! Skipping measurement."
		<< std::endl;
      return nullptr;
    }
  if(Verbosity() > 0)
    {
      std::cout << "Stream of found TPC surface: " << std::endl;
      surface->toStream(m_tGeometry->geoContext, std::cout);
    }

  /// Transformation of cluster to local surface coords
  /// Coords are r*phi relative to surface r-phi center, and z
  /// relative to surface z center
  /// lengths in mm
  Acts::Vector3D center = surface->center(m_actsGeometry->getGeoContext());
  Acts::Vector3D normal = surface->normal(m_actsGeometry->getGeoContext());
  
  double surfRadius = sqrt(center[0]*center[0] + center[1]*center[1]);
  double surfPhiCenter = atan2(center[1], center[0]);
  double surfRphiCenter = atan2(center[1], center[0]) * surfRadius;
  double surfZCenter = center[2];
  
  if (Verbosity() > 0)
  {
    std::cout << std::endl << "surface center readback:   x " << center[0]
              << " y " << center[1]  << " z " << center[2] << " radius " << surfRadius << std::endl;
    std::cout << "Surface normal vector : "<< normal(0) << ", " 
	      << normal(1) << ", " << normal(2) << std::endl;
    std::cout << " surface center  phi " << atan2(center[1], center[0]) 
              << " surface center r*phi " << surfRphiCenter
              << " surface center z  " << surfZCenter
              << std::endl;
  }

  Acts::Vector3D globalPos(x * Acts::UnitConstants::cm,
			   y * Acts::UnitConstants::cm,
			   z * Acts::UnitConstants::cm);
  
  surface->globalToLocal(m_actsGeometry->getGeoContext(), globalPos,
			 surface->normal(m_actsGeometry->getGeoContext()),
			 local2D);
  
  /// Test that Acts surface transforms correctly back
  Acts::Vector3D actsGlobal(0,0,0);
  surface->localToGlobal(m_actsGeometry->getGeoContext(), local2D, 
			 Acts::Vector3D(1,1,1), actsGlobal);

  if (Verbosity() > 0)
  {
    std::cout << "cluster readback (mm):  x " << x*Acts::UnitConstants::cm <<  " y " << y*Acts::UnitConstants::cm << " z " << z*Acts::UnitConstants::cm 
	      << " radius " << radius << std::endl;
    std::cout << " cluster phi " << clusPhi << " cluster z " << zTpc << " r*clusphi " << rClusPhi << std::endl;
    std::cout << " local phi " << clusPhi - surfPhiCenter
              << " local rphi " << rClusPhi-surfRphiCenter 
 	      << " local z " << zTpc - surfZCenter  << std::endl;
    std::cout << " acts local : " <<local2D(0) <<"  "<<local2D(1) << std::endl;
    std::cout << " sPHENIX global : " << x * 10 << "  " << y * 10 << "  " 
	      << z * 10 << "  " << std::endl;
    std::cout << " acts global : " << actsGlobal(0) << "  " << actsGlobal(1) 
	      << "  " << actsGlobal(2) << std::endl;
  }

  TMatrixD sPhenixLocalErr = transformCovarToLocal(clusPhi, worldErr);

  /// Get the 2D location covariance uncertainty for the cluster (y and z)
  localErr(Acts::eLOC_0, Acts::eLOC_0) = 
    sPhenixLocalErr[1][1] * Acts::UnitConstants::cm2;
  localErr(Acts::eLOC_1, Acts::eLOC_0) = 
    sPhenixLocalErr[2][1] * Acts::UnitConstants::cm2;
  localErr(Acts::eLOC_0, Acts::eLOC_1) = 
    sPhenixLocalErr[1][2] * Acts::UnitConstants::cm2;
  localErr(Acts::eLOC_1, Acts::eLOC_1) = 
    sPhenixLocalErr[2][2] * Acts::UnitConstants::cm2;


  if(Verbosity() > 0)
    {
      for (int i = 0; i < 3; ++i)
	{
	  for (int j = 0; j < 3; j++)
	    {
	      std::cout << "  " << i << " "  << j << " worldErr " << worldErr[i][j]  << " localErr " << sPhenixLocalErr[i][j] << std::endl;
	    }
	}
    }

  return surface;
}

Surface PHActsSourceLinks::getInttLocalCoords(Acts::Vector2D &local2D,
                                              Acts::BoundMatrix &localErr,
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
 
  CylinderGeomIntt *layerGeom =
    dynamic_cast<CylinderGeomIntt *>(m_geomContainerIntt->GetLayerGeom(layer));
 
  local = layerGeom->get_local_from_world_coords(ladderZId,
						 ladderPhiId,
						 world);
 
  Acts::Vector3D globalPos(x * Acts::UnitConstants::cm,
			   y * Acts::UnitConstants::cm,
			   z * Acts::UnitConstants::cm);
  
  surface->globalToLocal(m_actsGeometry->getGeoContext(), globalPos,
			 surface->normal(m_actsGeometry->getGeoContext()),
			 local2D);

  /// Test that Acts surface transforms correctly back
  Acts::Vector3D actsGlobal(0,0,0);
  surface->localToGlobal(m_actsGeometry->getGeoContext(), local2D, 
			 Acts::Vector3D(1,1,1), actsGlobal);

  Acts::Vector3D normal = surface->normal(m_actsGeometry->getGeoContext());
  
  if (Verbosity() > 0)
  {
    double segcent[3];
    layerGeom->find_segment_center(ladderZId, ladderPhiId, segcent);
    std::cout << "Acts surface normal vector: " << normal(0) << ", "
	      << normal(1) << ", " <<normal(2)<<std::endl;
    std::cout << "   segment center: " << segcent[0]
              << " " << segcent[1] << " " << segcent[2] << std::endl;
    std::cout << "   world; " << world[0] << " " << world[1]
              << " " << world[2] << std::endl;
    std::cout << "   local; " << local[0] << " " << local[1]
              << " " << local[2] << std::endl;
    std::cout << " acts local " << local2D(0) << "  " << local2D(1) << std::endl;
    std::cout << " sPHENIX global : " << x * 10 << "  " << y * 10 << "  " 
	      << z * 10 << "  " << std::endl;
    std::cout << " acts global : " << actsGlobal(0) << "  " << actsGlobal(1) 
	      << "  " << actsGlobal(2) << std::endl;
  }

  /// Get the local covariance error
  localErr = getInttCovarLocal(layer, ladderZId, ladderPhiId, worldErr);


  return surface;
}

Surface PHActsSourceLinks::getMvtxLocalCoords(Acts::Vector2D &local2D,
                                              Acts::BoundMatrix &localErr,
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

  CylinderGeom_Mvtx *layerGeom = 
    dynamic_cast<CylinderGeom_Mvtx *>(m_geomContainerMvtx->GetLayerGeom(layer));
 
  // this is just for checking - however it disagrees in x by much more than the pixel size - why? tilt not right?
  local = layerGeom->get_local_from_world_coords(staveId, 0, 0,
						 chipId,
						 world);


  Acts::Vector3D globalPos(x * Acts::UnitConstants::cm,
			   y * Acts::UnitConstants::cm,
			   z * Acts::UnitConstants::cm);
  
  surface->globalToLocal(m_actsGeometry->getGeoContext(), globalPos,
			 surface->normal(m_actsGeometry->getGeoContext()),
			 local2D);

  /// Test that Acts surface transforms correctly back
  Acts::Vector3D actsGlobal(0,0,0);
  surface->localToGlobal(m_actsGeometry->getGeoContext(), local2D, 
			 Acts::Vector3D(1,1,1), actsGlobal);

  Acts::Vector3D normal = surface->normal(m_actsGeometry->getGeoContext());

  if (Verbosity() > 0)
  {
    double segcent[3];
    std::cout << "Acts normal vector: "<<normal(0) << ", " << normal(1) 
	      << ", " << normal(2) << std::endl;
    layerGeom->find_sensor_center(staveId, 0, 0, chipId, segcent);
    std::cout << "   segment center: " << segcent[0] << " "
              << segcent[1] << " " << segcent[2] << std::endl;
    std::cout << "   world; " << world[0] << " "
              << world[1] << " " << world[2] << std::endl;
    std::cout << "   our local; " << local[0] * 10 << " "
              << local[1] * 10 << " " << local[2] * 10 << std::endl;
    std::cout << " acts local " << local2D(0) << "  " << local2D(1) << std::endl;
    std::cout << "Our global : " << x * 10 << "  " << y * 10 << "   " << z * 10 
	      << std::endl;
    std::cout << "Acts transform global : " << actsGlobal(0) << "  " << actsGlobal(1)
	      << "   " << actsGlobal(2) << std::endl;
  }

  // transform covariance matrix back to local coords on chip
  localErr = getMvtxCovarLocal(layer, staveId, chipId, worldErr);

  return surface;
}

void PHActsSourceLinks::addVerticesAsSourceLinks(PHCompositeNode *topNode,
						 unsigned int &hitId)
{
  SvtxVertexMap *vertexMap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");

  if(!vertexMap)
    {
      std::cout << PHWHERE << "Can't get vertex map to add source links. Vertices won't be added as source links"
		<< std::endl;
      return;
    }

  for(SvtxVertexMap::Iter vertexIter = vertexMap->begin();
      vertexIter != vertexMap->end();
      ++vertexIter)
    {

      const SvtxVertex *vertex = vertexIter->second;

      const Acts::Vector3D globalPos(vertex->get_x() * Acts::UnitConstants::cm,
				     vertex->get_y() * Acts::UnitConstants::cm,
				     vertex->get_z() * Acts::UnitConstants::cm);

      /// Make a perigee surface corresponding to the vertex for the 
      /// "measurement" to live on
      auto pSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(
				     globalPos);

      Acts::Vector2D loc(0,0); 
      Acts::BoundMatrix cov = Acts::BoundMatrix::Zero();

      pSurface->globalToLocal(m_actsGeometry->getGeoContext(),
			      globalPos,
			      pSurface->normal(m_actsGeometry->getGeoContext()),
			      loc);

      TMatrixD worldErr(3,3);
      for(int i = 0; i < 3; i++)
	for(int j = 0; j < 3; j++)
	  worldErr(i,j) = vertex->get_error(i, j);

      if(Verbosity() > 1)
	{
	  std::cout << "global cov xx, yy, zz: " << worldErr(0,0) 
		    << ", " << worldErr(1,1) << ", " << worldErr(2,2)
		    << std::endl;
	}

      const float phi = atan2(vertex->get_y(), vertex->get_x());
      
      TMatrixD localErr(3,3);
      localErr = transformCovarToLocal(phi, worldErr);

      if(Verbosity() > 1)
	{
	  std::cout << "local cov rphi, z " << localErr[0][0]
		    << ", " << localErr[2][2] << std::endl;
	}

      Acts::BoundVector location = Acts::BoundVector::Zero();
      location[Acts::eLOC_0] = loc(0);
      location[Acts::eLOC_1] = loc(1);
      
      cov(Acts::eLOC_0, Acts::eLOC_0) =
	localErr[0][0] * Acts::UnitConstants::cm2;
      cov(Acts::eLOC_1, Acts::eLOC_0) =
	localErr[2][0] * Acts::UnitConstants::cm2;
      cov(Acts::eLOC_0, Acts::eLOC_1) =
	localErr[0][2] * Acts::UnitConstants::cm2;
      cov(Acts::eLOC_1, Acts::eLOC_1) =
	localErr[2][2] * Acts::UnitConstants::cm2;
      
      if(Verbosity() > 1)
	{
	  Acts::Vector3D center = pSurface->center(m_actsGeometry->getGeoContext());
	  Acts::Vector3D normal = pSurface->normal(m_actsGeometry->getGeoContext());

	  std::cout << "Adding vertex as a measurement in track fitting"
		    << std::endl
		    << "global loc: (" << globalPos(0) << ", " << globalPos(1)
		    << ", " << globalPos(2) << ")" << std::endl
		    << "local pos (" << loc(0) << ", " << loc(1) << ") " 
		    << std::endl << "surface center: (" << center(0)
		    << ", " << center(1) << ", " << center(2) << ")" << std::endl
		    << "surface normal : (" << normal(0) << ", " << normal(1)
		    << ", " << normal(2) << ") " << std::endl;
	}

      SourceLink vertexSL(hitId, pSurface, location, cov);
      
      m_sourceLinks->insert(std::pair<unsigned int, SourceLink>(hitId,vertexSL));

      hitId++;
    }

  return;

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

  /// Do the same for the SourceLink container
  m_sourceLinks = findNode::getClass<std::map<unsigned int, SourceLink>>(topNode, "TrkrClusterSourceLinks");

  if (!m_sourceLinks)
  {
    m_sourceLinks = new std::map<unsigned int, SourceLink>;
    PHDataNode<std::map<unsigned int, SourceLink>>
        *sourceLinkNode = new PHDataNode<std::map<unsigned int, SourceLink>>(m_sourceLinks, "TrkrClusterSourceLinks");

    svtxNode->addNode(sourceLinkNode);
  }

  m_tGeometry = findNode::getClass<ActsTrackingGeometry>(topNode,"ActsTrackingGeometry");
  if(!m_tGeometry)
    {
      m_tGeometry = new ActsTrackingGeometry();
      PHDataNode<ActsTrackingGeometry> *tGeoNode = new PHDataNode<ActsTrackingGeometry>(m_tGeometry,
								 "ActsTrackingGeometry");
      svtxNode->addNode(tGeoNode);
    }


  return;
}

TGeoNode *PHActsSourceLinks::getNodeFromClusterMap(TrkrDefs::hitsetkey hitSetKey)
{
  std::map<TrkrDefs::hitsetkey, TGeoNode*> clusterNodeMap = m_actsGeometry->getTGeoNodeMap();

  /// Get the TGeoNode for this hit set key
  std::map<TrkrDefs::hitsetkey, TGeoNode *>::iterator mapIter;
  mapIter = clusterNodeMap.find(hitSetKey);
  TGeoNode *sensorNode;

  /// Make sure we found it
  if (mapIter != clusterNodeMap.end())
  {
    sensorNode = mapIter->second;
    if (Verbosity() > 0)
    {
      std::cout << "Found TGeoNode in clusterNodeMap for hitsetkey "
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
  std::map<TrkrDefs::hitsetkey, Surface> clusterSurfaceMap = 
    m_actsGeometry->getSurfaceMapSilicon();

  std::map<TrkrDefs::hitsetkey, Surface>::iterator
      surfaceIter;

  surfaceIter = clusterSurfaceMap.find(hitSetKey);

  /// Check to make sure we found the surface in the map
  if (surfaceIter != clusterSurfaceMap.end())
  {
    surface = surfaceIter->second;
    if (Verbosity() > 0)
    {
      std::cout << "Got surface pair " << surface->name()
                << " surface type " << surface->type()
                << std::endl;
      surface->toStream(m_tGeometry->geoContext, std::cout);
      std::cout<<std::endl;
    }
  }
  else
  {
    /// If it doesn't exit, return nullptr
    return nullptr;
  }

  return surface;
}

Acts::BoundMatrix PHActsSourceLinks::getMvtxCovarLocal(const unsigned int layer, 
						       const unsigned int staveId,
						       const unsigned int chipId, 
						       TMatrixD worldErr)
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

  if (Verbosity() > 0)
  {
    for (int i = 0; i < 3; ++i)
    {
      for (int j = 0; j < 3; ++j)
      {
        std::cout << "  " << i << "    " << j
                  << " world_err " << worldErr[i][j] << " local_err " << localErr[i][j]
                  << std::endl;
      }
    }
  }

  Acts::BoundMatrix matrix = Acts::BoundMatrix::Zero();

  /// Get the 2D location covariance uncertainty for the cluster (x and z) 
  matrix(Acts::eLOC_0, Acts::eLOC_0) = 
    localErr[0][0] * Acts::UnitConstants::cm2;
  matrix(Acts::eLOC_1, Acts::eLOC_0) = 
    localErr[2][0] * Acts::UnitConstants::cm2;
  matrix(Acts::eLOC_0, Acts::eLOC_1) = 
    localErr[0][2] * Acts::UnitConstants::cm2;
  matrix(Acts::eLOC_1, Acts::eLOC_1) = 
    localErr[2][2] * Acts::UnitConstants::cm2;

  return matrix;
}

Acts::BoundMatrix PHActsSourceLinks::getInttCovarLocal(const unsigned int layer, const unsigned int ladderZId, const unsigned int ladderPhiId, TMatrixD worldErr)
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

  Acts::BoundMatrix matrix = Acts::BoundMatrix::Zero();
  /// Get the 2D location covariance uncertainty for the cluster (y and z)
  matrix(Acts::eLOC_0, Acts::eLOC_0) = 
    localErr[1][1] * Acts::UnitConstants::cm2;
  matrix(Acts::eLOC_1, Acts::eLOC_0) = 
    localErr[2][1] * Acts::UnitConstants::cm2;
  matrix(Acts::eLOC_0, Acts::eLOC_1) = 
    localErr[1][2] * Acts::UnitConstants::cm2;
  matrix(Acts::eLOC_1, Acts::eLOC_1) = 
    localErr[2][2] * Acts::UnitConstants::cm2;

  return matrix;
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
