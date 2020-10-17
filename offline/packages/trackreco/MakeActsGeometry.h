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

#include "ActsTrackingGeometry.h"
#include "ActsSurfaceMaps.h"

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

namespace ActsExamples 
{
  class IBaseDetector;
  class IContextDecorator;
}

namespace Acts 
{
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
class MakeActsGeometry : public SubsysReco
{
 public:

  //! Default constructor
  MakeActsGeometry(const std::string& name = "MakeActsGeometry");

  //! Destructor
  ~MakeActsGeometry();

  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

  std::vector<std::shared_ptr<ActsExamples::IContextDecorator>> getContextDecorators()
    { return m_contextDecorators; }

  void setMagField(const std::string &magField)
    {m_magField = magField;}
  void setMagFieldRescale(double magFieldRescale)
    {m_magFieldRescale = magFieldRescale;}

  double getSurfStepPhi() {return m_surfStepPhi;}
  double getSurfStepZ() {return m_surfStepZ;}

 private:
  /// Main function to build all acts geometry for use in the fitting modules
  int buildAllGeometry(PHCompositeNode *topNode);
 
  //! Get all the nodes
  int getNodes(PHCompositeNode*);
  
  //!Create New nodes
  int createNodes(PHCompositeNode*);
  
  /// Functions to edit TGeoManager to include TPC boxes
  void setPlanarSurfaceDivisions();
  void editTPCGeometry(PHCompositeNode *topNode);
  void addActsTpcSurfaces(TGeoVolume *tpc_gas_vol, 
			  TGeoManager *geoManager);
  void addActsMicromegasSurfaces(int mm_layer, TGeoVolume *micromegas_vol, 
				 TGeoManager *geoManager);


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
  void makeMmMapPairs(TrackingVolumePtr &tpcVolume);
  
  /// Get subdetector hitsetkey from the local sensor unit coordinates
  TrkrDefs::hitsetkey getMvtxHitSetKeyFromCoords(unsigned int layer, 
						 std::vector<double> &world);
  TrkrDefs::hitsetkey getInttHitSetKeyFromCoords(unsigned int layer,
						 std::vector<double> &world);
  TrkrDefs::hitsetkey getTpcHitSetKeyFromCoords(std::vector<double> &world);
  TrkrDefs::hitsetkey getMmHitSetKeyFromCoords(std::vector<double> &world);

  /// Helper diagnostic function for identifying active layers in subdetectors
  void isActive(TGeoNode *gnode, int nmax_print);

  /// Makes map of TrkrHitSetKey<-->TGeoNode
  void makeTGeoNodeMap(PHCompositeNode *topNode);
  
  void unpackVolumes();

  /// Subdetector geometry containers for getting layer information
  PHG4CylinderGeomContainer* m_geomContainerMvtx = nullptr;
  PHG4CylinderGeomContainer* m_geomContainerInt = nullptrt;
  PHG4CylinderCellGeomContainer* m_geomContainerTpc = nullptr;

  TGeoManager* m_geoManager = nullptr;

  /// Acts Context decorators, which may contain e.g. calibration information
  std::vector<std::shared_ptr<ActsExamples::IContextDecorator> > 
    m_contextDecorators;

  /// Several maps that connect Acts world to sPHENIX G4 world 
  std::map<TrkrDefs::hitsetkey, TGeoNode*> m_clusterNodeMap;
  std::map<TrkrDefs::hitsetkey, Surface> m_clusterSurfaceMapSilicon;
  std::map<TrkrDefs::hitsetkey, std::vector<Surface>> m_clusterSurfaceMapTpcEdit;
  std::map<TrkrDefs::hitsetkey, std::vector<Surface>> m_clusterSurfaceMapMmEdit;
  
  /// These don't change, we are building the tpc this way!
  const static unsigned int m_nTpcLayers = 48;
  const unsigned int m_nTpcModulesPerLayer = 12;
  const unsigned int m_nTpcSides = 2;

  /// TPC Acts::Surface subdivisions
  double m_minSurfZ = 0.;
  double m_maxSurfZ = 105.5;
  unsigned int m_nSurfZ = 1;
  unsigned int m_nSurfPhi = 12;
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

  // Micromegas box surfaces use same phi and z segmentation as TPC, but layer details are different
  const static int m_nMmLayers = 2;
  const unsigned int m_mmLayerNumber[m_nMmLayers] = {55, 56};
  double m_mmLayerRadius[m_nMmLayers] = {82.2565, 82.6998};
  double m_mmLayerThickness[m_nMmLayers] = {0.3, 0.3};  // cm

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

  /// Structs to put on the node tree which carry around ActsGeom info
  ActsTrackingGeometry *m_actsGeometry = nullptr;
  ActsSurfaceMaps *m_surfMaps = nullptr;

  /// Verbosity value handed from PHActsSourceLinks
  int m_verbosity = 0;

  /// Magnetic field components to set Acts magnetic field
  std::string m_magField ="1.4" ;
  double m_magFieldRescale = -1.;

  bool m_buildMMs = false;

};

#endif
