
#ifndef TRACKRECO_PHACTSSOURCELINKS_H
#define TRACKRECO_PHACTSSOURCELINKS_H

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>

#include <string>
#include <map>
#include <vector>
#include <TMatrixDfwd.h>  

/// Acts includes to create all necessary definitions
#include <Acts/Utilities/Definitions.hpp>
#include <Acts/Utilities/BinnedArray.hpp>                    
#include <Acts/Utilities/Logger.hpp>   
#include <ACTFW/EventData/Track.hpp>       
#include "TrkrClusterSourceLink.hpp"

class PHCompositeNode;
class TrkrClusterContainer;
class TGeoNode;
class PHG4CylinderGeomContainer;
class PHG4CylinderCellGeomContainer;

namespace FW {
  class IBaseDetector;
}

namespace Acts {
  class Surface;
}


/**
 * This class is responsible for creating Acts TrkrClusterSourceLinks from
 * the SvtxClusters. The class creates a node of TrkrClusterSourceLinks and
 * places it on the node tree for other Acts modules to use (e.g. for use of
 * track propagation or track fitting)
 */
class PHActsSourceLinks : public SubsysReco
{

 public:
   /// Default constructor
   PHActsSourceLinks(const std::string& name = "PHActsSourceLinks");

  /// Destructor
  virtual ~PHActsSourceLinks() {}

  /// SubsysReco inherited functions
  int End(PHCompositeNode*);
  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);



 private:

  /**
   * Functions
   */
   /// Create nodes on tree if they don't yet exist
  void createNodeTree(PHCompositeNode *topNode);

  /// Get a TGeoNode from the m_clusterNodeMap
  TGeoNode* getNodeFromClusterMap(TrkrDefs::hitsetkey hitSetKey);

  /// Get a Surface from the m_surfaceNodeMap;
  std::shared_ptr<const Acts::Surface> getSurfaceFromClusterMap(TrkrDefs::hitsetkey hitSetKey);

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
  
 

  /**
   * Member variables
   */

  /// SvtxCluster node
  TrkrClusterContainer *m_clusterMap;

  /// Function to get necessary nodes off node tree
  int GetNodes(PHCompositeNode *topNode);

  /// Map relating arbitrary hitid to TrkrDef::cluskey for SourceLink
  std::map<TrkrDefs::cluskey, unsigned int> m_hitIdClusKey;

  /// Container for source links
  TrkrClusterSourceLinkContainer m_sourceLinks;

  /// Map relating hit set keys to TGeoNodes
  std::map<TrkrDefs::hitsetkey, TGeoNode*> *m_clusterNodeMap;

  /// Map relating hit set keys to Acts::Surfaces
  std::map<TrkrDefs::hitsetkey, std::shared_ptr<const Acts::Surface>> *m_clusterSurfaceMap;

  /// Get the TPC surface cluster map
  std::map<TrkrDefs::cluskey, std::shared_ptr<const Acts::Surface>> m_clusterSurfaceMapTpc;
  /// Tracking geometry objects
  PHG4CylinderGeomContainer* m_geomContainerMvtx;
  PHG4CylinderGeomContainer* m_geomContainerIntt;
  PHG4CylinderCellGeomContainer* m_geomContainerTpc;

  /// These don't change, we are building the tpc this way!
  const unsigned int m_nTpcLayers = 48;
  const unsigned int m_nTpcModulesPerLayer = 12;
  const unsigned int m_nTpcSides = 2;

  /// TPC surface subdivisions. These will come from the ActsGeometry helper class once it is finished
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
