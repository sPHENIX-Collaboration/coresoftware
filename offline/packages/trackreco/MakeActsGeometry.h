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

#include <ActsExamples/TGeoDetector/TGeoDetector.hpp>
#include <ActsExamples/Fitting/TrkrClusterFittingAlgorithm.hpp>
#include <ActsExamples/Plugins/BField/BFieldOptions.hpp>

#include <map>
#include <memory>            
#include <string>
#include <vector>

class PHCompositeNode;
class PHG4CylinderGeomContainer;
class PHG4CylinderCellGeomContainer;
class TGeoManager;
class TGeoNode;
class TGeoVolume;

namespace ActsExamples {
  class IBaseDetector;
  class IContextDecorator;
}

namespace Acts {
  class Surface;
}

using Surface = std::shared_ptr<const Acts::Surface>;
using TrackingGeometry = std::shared_ptr<const Acts::TrackingGeometry>;
using TrackingVolumePtr = std::shared_ptr<const Acts::TrackingVolume>;

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

  /// Main function to build all acts geometry for use in the fitting modules
  int buildAllGeometry(PHCompositeNode *topNode);
  
  /// Functions to edit TGeoManager to include TPC boxes
  void editTPCGeometry(PHCompositeNode *topNode);
  void addActsTpcSurfaces(TGeoVolume *tpc_gas_vol, 
			  TGeoManager *geoManager);
  void addActsMicromegasSurfaces(TGeoVolume *micromegas_vol, 
				 TGeoManager *geoManager);

  void setVerbosity(int verbosity)
    { m_verbosity = verbosity; }

  std::map<TrkrDefs::hitsetkey,Surface> getSurfaceMapSilicon()
    { return m_clusterSurfaceMapSilicon; }
  
  std::map<TrkrDefs::hitsetkey, std::vector<Surface>> getSurfaceMapTpc()
    { return m_clusterSurfaceMapTpcEdit; }
  
  std::map<TrkrDefs::hitsetkey, TGeoNode*> getTGeoNodeMap()
    { return m_clusterNodeMap; }
   
  std::vector<std::shared_ptr<ActsExamples::IContextDecorator>> getContextDecorators()
    { return m_contextDecorators; }
  
  /// Getters for acts geometry that is needed by fitter functions
  TrackingGeometry getTGeometry(){ return m_tGeometry; }
  ActsExamples::Options::BFieldVariant getMagField()
    { return m_magneticField; }
  Acts::MagneticFieldContext getMagFieldContext() 
    { return m_magFieldContext; }
  Acts::CalibrationContext getCalibContext() 
    { return m_calibContext; }
  Acts::GeometryContext getGeoContext() 
    { return m_geoCtxt; }
  
  /// Gets tpc surface from a cluster coordinate and hitsetkey. Necessary
  /// since there are many tpc surfaces per read out module
  Surface getTpcSurfaceFromCoords(TrkrDefs::hitsetkey hitsetkey, 
    std::vector<double> &world);

  void setMagField(std::string magField)
    {m_magField = magField;}
  void setMagFieldRescale(double magFieldRescale)
    {m_magFieldRescale = magFieldRescale;}

 private:
  
  //! Get all the nodes
  int getNodes(PHCompositeNode*);
  
  //!Create New nodes
  int createNodes(PHCompositeNode*);
  
  /// Silicon layers made by BuildSiliconLayers and its helper functions
  void buildActsSurfaces();

  /// Function that mimics ActsFW::GeometryExampleBase
  void makeGeometry(int argc, char* argv[], 
		    ActsExamples::IBaseDetector& detector);
  
  /// Get hitsetkey from TGeoNode for each detector geometry
  void getInttKeyFromNode(TGeoNode *gnode);
  void getMvtxKeyFromNode(TGeoNode *gnode);
  void getTpcKeyFromNode(TGeoNode *gnode);
  
  /// Make the Surface<-->TrkrDef::hitsetkey map pairs for each of 
  /// the various subdetectors
  void makeMvtxMapPairs(TrackingVolumePtr &mvtxVolume);
  void makeInttMapPairs(TrackingVolumePtr &inttVolume);
  void makeTpcMapPairs(TrackingVolumePtr &tpcVolume);
  
  /// Get subdetector hitsetkey from the local sensor unit coordinates
  TrkrDefs::hitsetkey getMvtxHitSetKeyFromCoords(unsigned int layer, 
						 std::vector<double> &world);
  TrkrDefs::hitsetkey getInttHitSetKeyFromCoords(unsigned int layer,
						 std::vector<double> &world);
  TrkrDefs::hitsetkey getTpcHitSetKeyFromCoords(std::vector<double> &world);

  /// Helper diagnostic function for identifying active layers in subdetectors
  void isActive(TGeoNode *gnode, int nmax_print);

  /// Makes map of TrkrHitSetKey<-->TGeoNode
  void makeTGeoNodeMap(PHCompositeNode *topNode);

  /// Subdetector geometry containers for getting layer information
  PHG4CylinderGeomContainer* m_geomContainerMvtx;  
  PHG4CylinderGeomContainer* m_geomContainerIntt;
  PHG4CylinderCellGeomContainer* m_geomContainerTpc;

  TGeoManager* m_geoManager; 

  /// Acts Context decorators, which may contain e.g. calibration information
  std::vector<std::shared_ptr<ActsExamples::IContextDecorator> > 
    m_contextDecorators;

  /// Several maps that connect Acts world to sPHENIX G4 world 
  std::map<TrkrDefs::hitsetkey, TGeoNode*> m_clusterNodeMap;
  std::map<TrkrDefs::hitsetkey, Surface> m_clusterSurfaceMapSilicon;
  std::map<TrkrDefs::hitsetkey, std::vector<Surface>> m_clusterSurfaceMapTpcEdit;
  std::map<TrkrDefs::cluskey, Surface> m_clusterSurfaceMapTpc;
  
  /// These don't change, we are building the tpc this way!
  const static unsigned int m_nTpcLayers = 48;
  const unsigned int m_nTpcModulesPerLayer = 12;
  const unsigned int m_nTpcSides = 2;

  /// TPC Acts::Surface subdivisions
  double m_minSurfZ;
  double m_maxSurfZ;
  unsigned int m_nSurfZ;
  unsigned int m_nSurfPhi;
  double m_surfStepPhi;
  double m_surfStepZ;
  double m_moduleStepPhi;
  double m_modulePhiStart;

  /// Debugger for printing out tpc active volumes
  int nprint_tpc;

  /// TPC TGeoManager editing box surfaces subdivisions
  const static int m_nTpcSectors = 3;
  const double m_minRadius[m_nTpcSectors] = {30.0, 40.0, 60.0};
  const double m_maxRadius[m_nTpcSectors] = {40.0, 60.0, 77.0};
  double layer_thickness_sector[m_nTpcSectors] = {0};
  double m_layerRadius[m_nTpcLayers] = {0};
  double m_layerThickness[m_nTpcLayers] = {0};

  // Micromegas box surfacrs use same phi and z segmentation as TPC, but layer details are different
  const static int m_nMmLayers = 2;
  double m_mmLayerRadius[m_nMmLayers] = {82.2565, 82.6998};
  double m_mmLayerThickness[m_nMmLayers] = {0.3};

  // Spaces to prevent boxes from touching when placed
  const double half_width_clearance_thick = 0.4999;
  const double half_width_clearance_phi = 0.4999;
  const double half_width_clearance_z = 0.4999;

  /// The acts geometry object
  TGeoDetector m_detector;

  /// Acts geometry objects that are needed to create (for example) the fitter
  TrackingGeometry m_tGeometry;
  ActsExamples::Options::BFieldVariant m_magneticField;
  Acts::GeometryContext  m_geoCtxt;  
  Acts::CalibrationContext m_calibContext;
  Acts::MagneticFieldContext m_magFieldContext;

  /// Verbosity value handed from PHActsSourceLinks
  int m_verbosity;

  /// Magnetic field components to set Acts magnetic field
  std::string m_magField;
  double m_magFieldRescale;
};

#endif
