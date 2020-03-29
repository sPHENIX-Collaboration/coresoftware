
#ifndef TRACKRECO_PHACTSSOURCELINKS_H
#define TRACKRECO_PHACTSSOURCELINKS_H

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>

#include <TMatrixDfwd.h>
#include <map>
#include <string>
#include <vector>

/// Acts includes to create all necessary definitions
#include <Acts/Utilities/BinnedArray.hpp>
#include <Acts/Utilities/Definitions.hpp>
#include <Acts/Utilities/Logger.hpp>

#include <Acts/EventData/MeasurementHelpers.hpp>
#include <Acts/Geometry/TrackingGeometry.hpp>
#include <Acts/MagneticField/MagneticFieldContext.hpp>
#include <Acts/Utilities/CalibrationContext.hpp>


#include <ACTFW/EventData/Track.hpp>
#include <ACTFW/EventData/TrkrClusterSourceLink.hpp>
#include <ACTFW/Plugins/BField/BFieldOptions.hpp>

class PHCompositeNode;
class TrkrClusterContainer;
class TrkrCluster;
class TGeoNode;
class PHG4CylinderGeomContainer;
class PHG4CylinderCellGeomContainer;

namespace FW
{
class IBaseDetector;
}

namespace Acts
{
class Surface;
}

using Surface = std::shared_ptr<const Acts::Surface>;
using SourceLink = FW::Data::TrkrClusterSourceLink;

/**
 * A struct that contains the necessary geometry objects that the fitter
 * needs in PHActsTrkFitter. To be put on the node tree
 */
struct FitCfgOptions
{
  /// Two constructor options
  FitCfgOptions() {}
  FitCfgOptions(std::shared_ptr<const Acts::TrackingGeometry> tGeo,
                FW::Options::BFieldVariant mag,
                Acts::CalibrationContext calib,
                Acts::GeometryContext geo,
                Acts::MagneticFieldContext magField)
    : tGeometry(tGeo)
    , magField(mag)
    , calibContext(calib)
    , geoContext(geo)
    , magFieldContext(magField)
  {
  }

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

/**
 * This class is responsible for creating Acts TrkrClusterSourceLinks from
 * the SvtxClusters. The class creates a node of TrkrClusterSourceLinks and
 * places it on the node tree for other Acts modules to use (e.g. for use of
 * track propagation or track fitting). SourceLinks are created from the 
 * SvtxClusterContainer contents (SvtxClusters) and are put into a special
 * SourceLinkContainer that is defined in TrkrClusterSourceLink.hpp
 */
class PHActsSourceLinks : public SubsysReco
{
 public:
  /// Default constructor
  PHActsSourceLinks(const std::string &name = "PHActsSourceLinks");

  /// Destructor
  virtual ~PHActsSourceLinks() {}

  /// SubsysReco inherited functions
  int End(PHCompositeNode *topNode);
  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);

 private:
  /**
   * Functions
   */
  /// Create sourceLink related nodes on tree if they don't yet exist
  void createNodes(PHCompositeNode *topNode);

  /// Function to get necessary nodes off node tree
  int getNodes(PHCompositeNode *topNode);

  /// Get a TGeoNode from the m_clusterNodeMap
  TGeoNode *getNodeFromClusterMap(TrkrDefs::hitsetkey hitSetKey);

  /// Get a Surface from the m_surfaceNodeMap;
  Surface getSurfaceFromClusterMap(TrkrDefs::hitsetkey hitSetKey);

  /// Get the local covariance matrix for the Mvtx cluster in a form Acts
  /// can accept
  TMatrixD getMvtxCovarLocal(const unsigned int layer,
                             const unsigned int staveid,
                             const unsigned int chipid,
                             TMatrixD world_err);

  /// Get the local covariance matrix for the Intt cluster in a form Acts
  /// can accept
  TMatrixD getInttCovarLocal(const unsigned int layer,
                             const unsigned int ladderZId,
                             const unsigned int ladderPhiId,
                             TMatrixD worldErr);

  /// Transform the covariance to the local coordinate frame
  TMatrixD transformCovarToLocal(const double ladderPhi, TMatrixD worldErr);

  /// Function which returns MVTX local coordinates and error, as well as
  /// corresponding surface
  Surface getMvtxLocalCoords(double (&local2D)[2], TMatrixD &localErr,
                             const TrkrCluster *cluster,
                             const TrkrDefs::cluskey clusKey);

  /// Same as above, except for INTT
  Surface getInttLocalCoords(double (&local2D)[2], TMatrixD &localErr,
                             const TrkrCluster *cluster,
                             const TrkrDefs::cluskey clusKey);

  /// Same as above, except for TPC
  Surface getTpcLocalCoords(double (&local2D)[2], TMatrixD &localErr,
                            const TrkrCluster *cluster,
                            const TrkrDefs::cluskey clusKey);

  /**
   * Member variables
   */

  /// SvtxCluster node
  TrkrClusterContainer *m_clusterMap;

  /// Map relating arbitrary hitid to TrkrDef::cluskey for SourceLink, to be put
  /// on node tree by this module
  std::map<TrkrDefs::cluskey, unsigned int> *m_hitIdClusKey;

  /// Map for source hitid:sourcelink, to be put on node tree by this module
  std::map<unsigned int, SourceLink> *m_sourceLinks;

  /// Map relating hit set keys to TGeoNodes
  std::map<TrkrDefs::hitsetkey, TGeoNode *> m_clusterNodeMap;

  /// Map relating hit set keys to Acts::Surfaces
  std::map<TrkrDefs::hitsetkey, Surface> m_clusterSurfaceMap;

  /// Get the TPC surface cluster map
  std::map<TrkrDefs::cluskey, Surface> m_clusterSurfaceMapTpc;

  /// The fit cfg options, created in MakeActsGeometry, to be put on node tree
  /// for PHActsTrkFitter
  FitCfgOptions *m_fitCfgOptions;

  /// Tracking geometry objects
  PHG4CylinderGeomContainer *m_geomContainerMvtx;
  PHG4CylinderGeomContainer *m_geomContainerIntt;
  PHG4CylinderCellGeomContainer *m_geomContainerTpc;

  /// These don't change, we are building the tpc this way!
  const unsigned int m_nTpcLayers = 48;
  const unsigned int m_nTpcModulesPerLayer = 12;
  const unsigned int m_nTpcSides = 2;

  /// TPC surface subdivisions. These come from MakeActsGeometry
  double m_minSurfZ;
  double m_maxSurfZ;
  unsigned int m_nSurfZ;
  unsigned int m_nSurfPhi;
  double m_surfStepPhi;
  double m_surfStepZ;
  double m_moduleStepPhi;
  double m_modulePhiStart;
};

#endif
