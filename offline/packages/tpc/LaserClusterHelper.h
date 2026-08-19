#ifndef TPC_LASERCLUSTERHELPER_H
#define TPC_LASERCLUSTERHELPER_H

#include <trackbase/ActsGeometry.h>
#include <trackbase/TrkrDefs.h>

#include <array>
#include <memory>

class ActsGeometry;
class LaserCluster;
class PHCompositeNode;
class PHG4TpcGeomContainer;
class PHGarfield;
class TPolyLine3D;

class LaserClusterHelper
{
  public:
    LaserClusterHelper ();
    ~LaserClusterHelper();

    void loadNodes(PHCompositeNode *topNode);
  
    Acts::Vector3 getHitPosition(TrkrDefs::hitsetkey, TrkrDefs::hitkey) const;
    Acts::Vector3 getClusterCentroid(LaserCluster*) const;
    std::array<double, 3> getClusterHardwareCentroid(LaserCluster*) const;
    Acts::Vector3 getClusterCentroidWithPHGarfield(LaserCluster*) const;

    void set_useZ(bool use) { m_useZ = use; }
    void set_useGlobal(bool use) { m_useGlobal = use; }
    void set_garfield_cmvoltage(double use) { m_garfield_cmvoltage = use; }
    void set_garfield_zerofield(bool use) { m_garfield_zerofield = use; }
    void set_garfield_spacechargescaleside0(double use) { m_garfield_spacechargescaleside0 = use; }
    void set_garfield_spacechargescaleside1(double use) { m_garfield_spacechargescaleside1 = use; }
    void set_garfield_stepns(double use) { m_garfield_stepns = use; }
    void set_garfield_spacechargefieldmap(std::string use) { m_spacechargefieldmap = use; }
  private:
    ActsGeometry *m_tGeometry{nullptr};
    PHG4TpcGeomContainer *m_geom_container{nullptr};
    std::unique_ptr<PHGarfield> m_phgarfield;

    bool m_useZ{false};
    bool m_useGlobal{true};

    double m_garfield_cmvoltage{380.0};
    bool m_garfield_zerofield{false};
    double m_garfield_spacechargescaleside0{1.0};
    double m_garfield_spacechargescaleside1{1.0};
    double m_garfield_stepns{50.0};
    std::string m_spacechargefieldmap{"/sphenix/user/dloomis/Distortions/garfield_fields/sphenix_rossegger_garfield_field.root"};


};

#endif
