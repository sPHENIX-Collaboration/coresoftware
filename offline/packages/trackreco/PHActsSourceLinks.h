#ifndef TRACKRECO_PHACTSSOURCELINKS_H
#define TRACKRECO_PHACTSSOURCELINKS_H

#include <trackbase/ActsTrackingGeometry.h>
#include "ActsSurfaceMaps.h"

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>

#include <TMatrixDfwd.h>
#include <map>
#include <string>
#include <vector>

#include <boost/bimap.hpp>

/// Acts includes to create all necessary definitions
#include <Acts/Utilities/BinnedArray.hpp>
#include <Acts/Utilities/Definitions.hpp>
#include <Acts/Utilities/Logger.hpp>

#include <Acts/EventData/MeasurementHelpers.hpp>
#include <Acts/Geometry/TrackingGeometry.hpp>
#include <Acts/MagneticField/MagneticFieldContext.hpp>
#include <Acts/Utilities/CalibrationContext.hpp>

#include <ActsExamples/EventData/Track.hpp>
#include <ActsExamples/EventData/TrkrClusterSourceLink.hpp>
#include <ActsExamples/Plugins/BField/BFieldOptions.hpp>


class PHCompositeNode;
class TrkrClusterContainer;
class TrkrCluster;
class TGeoNode;
class PHG4CylinderGeomContainer;
class PHG4CylinderCellGeomContainer;


namespace ActsExamples
{
  class IBaseDetector;
}

namespace Acts
{
  class Surface;
}

using Surface = std::shared_ptr<const Acts::Surface>;
using SourceLink = ActsExamples::TrkrClusterSourceLink;

typedef boost::bimap<TrkrDefs::cluskey, unsigned int> CluskeyBimap;

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
  int ResetEvent(PHCompositeNode *topNode);
  void useVertexAsMeasurement(bool useVertexMeasurement)
    {m_useVertexMeasurement = useVertexMeasurement;}
  void setMagField(const std::string &magField)
    {m_magField = magField;}
  void setMagFieldRescale(double magFieldRescale)
    {m_magFieldRescale = magFieldRescale;}

  void SetUseTruthClusters(bool setit){_use_truth_clusters = setit;}
 
 private:
  /**
   * Functions
   */
  /// Create sourceLink related nodes on tree if they don't yet exist
  void createNodes(PHCompositeNode *topNode);

  /// Function to get necessary nodes off node tree
  int getNodes(PHCompositeNode *topNode);

  /// Get a TGeoNode from the m_clusterNodeMap
  TGeoNode* getNodeFromClusterMap(TrkrDefs::hitsetkey hitSetKey);

  /// Get a Surface from the m_surfaceNodeMap;
  Surface getSurfaceFromClusterMap(TrkrDefs::hitsetkey hitSetKey);

  /// Get the local covariance matrix for the Mvtx cluster in a form Acts
  /// can accept
  Acts::BoundMatrix getMvtxCovarLocal(const unsigned int layer,
                             const unsigned int staveid,
                             const unsigned int chipid,
                             TMatrixD world_err);

  /// Get the local covariance matrix for the Intt cluster in a form Acts
  /// can accept
  Acts::BoundMatrix getInttCovarLocal(const unsigned int layer,
                             const unsigned int ladderZId,
                             const unsigned int ladderPhiId,
                             TMatrixD worldErr);

  /// Transform the covariance to the local coordinate frame
  TMatrixD transformCovarToLocal(const double ladderPhi, TMatrixD worldErr);

  /// Function which returns MVTX local coordinates and error, as well as
  /// corresponding surface
  Surface getMvtxLocalCoords(Acts::Vector2D &local2D, 
			     Acts::BoundMatrix &localErr,
                             const TrkrCluster *cluster,
                             const TrkrDefs::cluskey clusKey);

  /// Same as above, except for INTT
  Surface getInttLocalCoords(Acts::Vector2D &local2D, 
			     Acts::BoundMatrix &localErr,
                             const TrkrCluster *cluster,
                             const TrkrDefs::cluskey clusKey);

  /// Same as above, except for TPC
  Surface getTpcLocalCoords(Acts::Vector2D &local2D, 
			    Acts::BoundMatrix &localErr,
                            const TrkrCluster *cluster,
                            const TrkrDefs::cluskey clusKey);

  Surface getMmLocalCoords(Acts::Vector2D &local2D,
                                             Acts::BoundMatrix &localErr,
                                             const TrkrCluster *cluster,
			   const TrkrDefs::cluskey clusKey);

  void addVerticesAsSourceLinks(PHCompositeNode *topNode,
				unsigned int &hitId);

  /// Gets tpc surface from a cluster coordinate and hitsetkey. Necessary
  /// since there are many tpc surfaces per read out module
  Surface getTpcSurfaceFromCoords(TrkrDefs::hitsetkey hitsetkey, 
    std::vector<double> &world);

  Surface getMmSurfaceFromCoords(TrkrDefs::hitsetkey hitsetkey, 
    std::vector<double> &world);


  /**
   * Member variables
   */

  bool m_useVertexMeasurement = false;
  bool _use_truth_clusters = false;

  /// SvtxCluster node
  TrkrClusterContainer *m_clusterMap = nullptr;

  /// Map relating arbitrary hitid to TrkrDef::cluskey for SourceLink, to be put
  /// on node tree by this module
  CluskeyBimap *m_hitIdClusKey;

  /// Map for source hitid:sourcelink, to be put on node tree by this module
  std::map<unsigned int, SourceLink> *m_sourceLinks;

  /// Magnetic field components to set Acts magnetic field
  std::string m_magField = "1.4";
  double m_magFieldRescale = -1.;

  /// Tracking geometry objects
  PHG4CylinderGeomContainer *m_geomContainerMvtx = nullptr;
  PHG4CylinderGeomContainer *m_geomContainerIntt = nullptr;
  PHG4CylinderCellGeomContainer *m_geomContainerTpc = nullptr;

  ActsTrackingGeometry *m_tGeometry = nullptr;
  ActsSurfaceMaps *m_surfMaps = nullptr;
};

#endif
