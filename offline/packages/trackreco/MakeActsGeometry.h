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
#include <Acts/Utilities/BinnedArray.hpp>                       // for Binne...
#include <Acts/Utilities/Logger.hpp>                            // for getDe...
#include <Acts/EventData/MeasurementHelpers.hpp>  // for GeometryContext
#include <Acts/Geometry/TrackingGeometry.hpp>

#include <ACTFW/TGeoDetector/TGeoDetector.hpp>
#include <ACTFW/Fitting/TrkrClusterFittingAlgorithm.hpp>
#include <ACTFW/Plugins/BField/BFieldOptions.hpp>

#include <map>
#include <memory>                // for shared_ptr
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
  class TrackingGeometry;
}


/**
 * A struct that contains the necessary geometry objects that the fitter
 * needs in PHActsTrkFitter. To be put on the node tree
 */
struct FitCfgOptions 
{
  /// Two constructor options
  FitCfgOptions(){}
  FitCfgOptions(std::shared_ptr<const Acts::TrackingGeometry> tGeo,
	        FW::Options::BFieldVariant mag,
		Acts::CalibrationContext calib,
		Acts::GeometryContext geo,
		Acts::MagneticFieldContext magField)
  : tGeometry(tGeo)
  , magField (mag)
  , calibContext(calib)
  , geoContext(geo)
  , magFieldContext(magField)
  {}

  /// Acts tracking geometry
  std::shared_ptr<const Acts::TrackingGeometry> tGeometry;

  /// Acts magnetic field
  FW::Options::BFieldVariant magField;

  /// Acts calibration context, grabbed from geometry building
  Acts::CalibrationContext calibContext;
  
  /// Acts geometry context, grabbed from geometry building
  Acts::GeometryContext geoContext;

  /// Acts magnetic field context, grabbed from geometry building
  Acts::MagneticFieldContext magFieldContext;
};


class MakeActsGeometry
{
 public:

  //! Default constructor
  MakeActsGeometry(const std::string& name = "MakeActsGeometry");

  //! dtor
  ~MakeActsGeometry();

  int BuildAllGeometry(PHCompositeNode *topNode);
  void SetVerbosity(int verbosity)
  {
    _verbosity = verbosity;
  }

  std::map<TrkrDefs::hitsetkey,std::shared_ptr<const Acts::Surface>> getSurfaceMapSilicon()
    {
      return _cluster_surface_map_silicon;
    }
  
  std::map<TrkrDefs::cluskey, std::shared_ptr<const Acts::Surface>> getSurfaceMapTpc()
    {
      return _cluster_surface_map_tpc;
    }
  
  std::map<TrkrDefs::hitsetkey, TGeoNode*> getNodeMap()
    {
      return _cluster_node_map;
    }
  
  double getSurfStepZ()
  {
    return SurfStepZ;
  }
  double getSurfStepPhi()
  {
    return SurfStepPhi;
  }
  
  double getModuleStepPhi()
  {
    return ModuleStepPhi;
  }
  
  double getModulePhiStart()
  {
    return ModulePhiStart;
  }
  
  Acts::GeometryContext  getGeoContext()
    {
      return geo_ctxt;
    }
  
  std::vector<std::shared_ptr<FW::IContextDecorator> > getContextDecorators()
    {
      return contextDecorators;
    }
  
  double getMinSurfZ() { return MinSurfZ; }
  double getMaxSurfZ() { return MaxSurfZ; }
  double getNSurfZ() { return NSurfZ; }
  double getNSurfPhi() { return NSurfPhi; }
  
  //Flags of different kinds of outputs
  enum Flag
  {
    //all disabled
    NONE = 0,
  };
    
  FitCfgOptions* getFitCfgOptions() { return m_fitCfgOptions; }

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
  TrkrDefs::hitsetkey GetMvtxHitSetKeyFromCoords(unsigned int layer, std::vector<double> &world);
  TrkrDefs::hitsetkey GetInttHitSetKeyFromCoords(unsigned int layer, std::vector<double> &world);
  void isActive(TGeoNode *gnode);

  void BuildTpcSurfaceMap();

  void MakeTGeoNodeMap(PHCompositeNode*);

  PHG4CylinderGeomContainer* _geom_container_mvtx;
  PHG4CylinderGeomContainer* _geom_container_intt;
  PHG4CylinderCellGeomContainer* _geom_container_tpc;

  TGeoManager* _geomanager;

  Acts::GeometryContext  geo_ctxt;

  std::vector<std::shared_ptr<FW::IContextDecorator> > contextDecorators;

  /// Several maps that connect Acts world to sPHENIX G4 world 
  std::map<TrkrDefs::hitsetkey, TGeoNode*> _cluster_node_map;
  std::map<TrkrDefs::hitsetkey,std::shared_ptr<const Acts::Surface>> _cluster_surface_map_silicon;
  std::map<TrkrDefs::cluskey, std::shared_ptr<const Acts::Surface>> _cluster_surface_map_tpc;

  // TPC surface subdivisions
  double MinSurfZ;
  double MaxSurfZ;
  unsigned int NSurfZ;
  unsigned int NSurfPhi;
  double SurfStepPhi;
  double SurfStepZ;
  double ModuleStepPhi;
  double ModulePhiStart;

  // these don't change, we are building the tpc this way!
  const unsigned int NTpcLayers = 48;
  const unsigned int NTpcModulesPerLayer = 12;
  const unsigned int NTpcSides = 2;


  Acts::Logging::Level logLevel;
  // The acts geometry object
  TGeoDetector detector;

 
  FitCfgOptions *m_fitCfgOptions;
  int _verbosity;
};

#endif
