#ifndef TPC_LASERCLUSTERHELPER_H
#define TPC_LASERCLUSTERHELPER_H

#include <trackbase/ActsGeometry.h>
#include <trackbase/TrkrDefs.h>

class ActsGeometry;
class LaserCluster;
class PHCompositeNode;
class PHG4TpcGeomContainer;

class LaserClusterHelper
{
  public:
    LaserClusterHelper () = default;

    void loadNodes(PHCompositeNode *topNode);
  
    Acts::Vector3 getHitPosition(TrkrDefs::hitsetkey, TrkrDefs::hitkey) const;
    Acts::Vector3 getClusterCentroid(LaserCluster*) const;

    void set_useZ(bool use) { m_useZ = use; }
    void set_useGlobal(bool use) { m_useGlobal = use; }
  private:

    ActsGeometry *m_tGeometry{nullptr};
    PHG4TpcGeomContainer *m_geom_container{nullptr};

    bool m_useZ{false};
    bool m_useGlobal{true};

};

#endif
