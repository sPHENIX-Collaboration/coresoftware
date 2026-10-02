#ifndef TPC_LASERCLUSTERHELPER_H
#define TPC_LASERCLUSTERHELPER_H

#include <fun4all/SubsysReco.h>


#include <trackbase/ActsGeometry.h>
#include <trackbase/TrkrDefs.h>

#include <array>

class ActsGeometry;
class LaserCluster;
class PHCompositeNode;
class PHG4TpcGeomContainer;

class LaserClusterHelper : public SubsysReco
{
  public:
    LaserClusterHelper () = default;

    void loadNodes(PHCompositeNode *topNode);
  
    Acts::Vector3 getHitPosition(TrkrDefs::hitsetkey, TrkrDefs::hitkey) const;
    Acts::Vector3 getClusterCentroid(LaserCluster*) const;
    std::array<double, 3> getClusterHardwareCentroid(LaserCluster*) const;

    void set_useZ(bool use) { m_useZ = use; }
    void set_useGlobal(bool use) { m_useGlobal = use; }
    void set_useDouble(bool use) { m_useDouble = use; }
  private:

    ActsGeometry *m_tGeometry{nullptr};
    PHG4TpcGeomContainer *m_geom_container{nullptr};

    bool m_useZ{false};
    bool m_useGlobal{true};
    bool m_useDouble{false};

};

#endif
