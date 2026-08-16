#ifndef TPC_LASERCLUSTERHELPER_H
#define TPC_LASERCLUSTERHELPER_H

#include <trackbase/ActsGeometry.h>
#include <trackbase/TrkrDefs.h>

#include <array>

class ActsGeometry;
class LaserCluster;
class PHCompositeNode;
class PHG4TpcGeomContainer;
class PHGarfield;
class TPolyLine3D;

class LaserClusterHelper
{
  public:
    LaserClusterHelper () = default;

    void loadNodes(PHCompositeNode *topNode);
  
    Acts::Vector3 getHitPosition(TrkrDefs::hitsetkey, TrkrDefs::hitkey) const;
    Acts::Vector3 getClusterCentroid(LaserCluster*) const;
    std::array<double, 3> getClusterHardwareCentroid(LaserCluster*) const;
    Acts::Vector3 getClusterCentroidWithPHGarfield(LaserCluster*) const;

    void set_useZ(bool use) { m_useZ = use; }
    void set_useGlobal(bool use) { m_useGlobal = use; }
  private:

    bool InterpolateOrClosestAtZ(TPolyLine3D* path, double zTarget, double& rOut, double& phiOut) const;

    PHGarfield* m_phgarfield{nullptr};
    ActsGeometry *m_tGeometry{nullptr};
    PHG4TpcGeomContainer *m_geom_container{nullptr};

    bool m_useZ{false};
    bool m_useGlobal{true};

};

#endif
