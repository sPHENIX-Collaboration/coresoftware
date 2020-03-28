/*!
 *  \file		MakeActsGeometry.h
 *  \brief		Make Acts geometry from sPHENIX TGeometry
 *  \details	Make Acts geometry from sPHENIX TGeometry
 *  \author		Tony Frawley <afrawley@fsu.edu>
 */

#ifndef MAKE_ACTS_GEOMETRY_H
#define MAKE_ACTS_GEOMETRY_H

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>

#include <Acts/Utilities/Definitions.hpp>
#include <Acts/Utilities/BinnedArray.hpp>                      
#include <Acts/Utilities/Logger.hpp>                           
#include <Acts/EventData/MeasurementHelpers.hpp> 
#include <Acts/Geometry/TrackingGeometry.hpp>
#include <Acts/MagneticField/MagneticFieldContext.hpp>
#include <Acts/Utilities/CalibrationContext.hpp>

#include <ACTFW/TGeoDetector/TGeoDetector.hpp>
#include <ACTFW/Fitting/TrkrClusterFittingAlgorithm.hpp>
#include <ACTFW/Plugins/BField/BFieldOptions.hpp>

#include <map>
#include <memory>            
#include <string>
#include <vector>

class PHCompositeNode;
class PHG4CylinderGeomContainer;
class PHG4CylinderCellGeomContainer;
class TGeoManager;
class TGeoNode;

namespace FW {
  class IBaseDetector;
  class IContextDecorator;
}

namespace Acts {
  class Surface;
}

using Surface = std::shared_ptr<const Acts::Surface>;
using TrackingGeometry = std::shared_ptr<const Acts::TrackingGeometry>;


/**
 * This class is responsible for building the ActsGeometry from the sPHENIX
 * TGeometry. The code largely follows examples within the ACTFW code, 
 * specifically GeometryExampleBase.cpp. 
 * This class puts several nodes on the node tree which relate Acts::Surfaces
 * to sPHENIX TGeo objects, for use in building Acts SourceLinks.
 */
class MakeActsGeometry
{
 public:

  //! Default constructor
  MakeActsGeometry(const std::string& name = "MakeActsGeometry");

  //! Destructor
  ~MakeActsGeometry();

  int BuildAllGeometry(PHCompositeNode *topNode);
  
  void SetVerbosity(int verbosity)
  { m_verbosity = verbosity; }

  std::map<TrkrDefs::hitsetkey,Surface> getSurfaceMapSilicon()
    { return m_clusterSurfaceMapSilicon; }
  
  std::map<TrkrDefs::cluskey, Surface> getSurfaceMapTpc()
    { return m_clusterSurfaceMapTpc; }
  
  std::map<TrkrDefs::hitsetkey, TGeoNode*> getNodeMap()
    { return m_clusterNodeMap; }
  
  double getSurfStepZ() { return m_surfStepZ; }
  double getSurfStepPhi() { return m_surfStepPhi; }
  double getModuleStepPhi() { return m_moduleStepPhi; }
  double getModulePhiStart() { return m_modulePhiStart; }
  Acts::GeometryContext getGeoContext() { return m_geoCtxt; }
  
  std::vector<std::shared_ptr<FW::IContextDecorator>> getContextDecorators()
    { return m_contextDecorators; }
  
  double getMinSurfZ() { return m_minSurfZ; }
  double getMaxSurfZ() { return m_maxSurfZ; }
  double getNSurfZ() { return m_nSurfZ; }
  double getNSurfPhi() { return m_nSurfPhi; }     
  TrackingGeometry getTGeometry(){ return m_tGeometry; }
  FW::Options::BFieldVariant getMagField(){ return m_magneticField; }
  Acts::MagneticFieldContext getMagFieldContext() { return m_magFieldContext; }
  Acts::CalibrationContext getCalibContext() { return m_calibContext; }

 private:
  
  //! Get all the nodes
  int GetNodes(PHCompositeNode*);
  
  //!Create New nodes
  int CreateNodes(PHCompositeNode*);
  
  // silicon layers made by BuildSiliconLayers and its helper functions
  void BuildSiliconLayers();

  int MakeSiliconGeometry(int argc, char* argv[], FW::IBaseDetector& detector);
  
  void getInttKeyFromNode(TGeoNode *gnode);
  
  void getMvtxKeyFromNode(TGeoNode *gnode);
  
  TrkrDefs::hitsetkey GetMvtxHitSetKeyFromCoords(unsigned int layer, 
						 std::vector<double> &world);
  
  TrkrDefs::hitsetkey GetInttHitSetKeyFromCoords(unsigned int layer,
						 std::vector<double> &world);
  
  void isActive(TGeoNode *gnode);

  void BuildTpcSurfaceMap();

  void MakeTGeoNodeMap(PHCompositeNode*);

  PHG4CylinderGeomContainer* m_geomContainerMvtx;
  
  PHG4CylinderGeomContainer* m_geomContainerIntt;
  
  PHG4CylinderCellGeomContainer* m_geomContainerTpc;

  TGeoManager* m_geoManager;

  Acts::GeometryContext  m_geoCtxt;

  std::vector<std::shared_ptr<FW::IContextDecorator> > m_contextDecorators;

  /// Several maps that connect Acts world to sPHENIX G4 world 
  std::map<TrkrDefs::hitsetkey, TGeoNode*> m_clusterNodeMap;
  std::map<TrkrDefs::hitsetkey,Surface> m_clusterSurfaceMapSilicon;
  std::map<TrkrDefs::cluskey, Surface> m_clusterSurfaceMapTpc;

  // TPC surface subdivisions
  double m_minSurfZ;
  double m_maxSurfZ;
  unsigned int m_nSurfZ;
  unsigned int m_nSurfPhi;
  double m_surfStepPhi;
  double m_surfStepZ;
  double m_moduleStepPhi;
  double m_modulePhiStart;

  // these don't change, we are building the tpc this way!
  const unsigned int m_nTpcLayers = 48;
  const unsigned int m_nTpcModulesPerLayer = 12;
  const unsigned int m_nTpcSides = 2;

  // The acts geometry object
  TGeoDetector m_detector;

  TrackingGeometry m_tGeometry;
  FW::Options::BFieldVariant m_magneticField;
  Acts::CalibrationContext m_calibContext;
  Acts::MagneticFieldContext m_magFieldContext;

  int m_verbosity;
};

#endif
