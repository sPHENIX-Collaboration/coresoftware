#include "LaserClusterHelper.h"
 
#include <phgarfield/PHGarfield.h>

#include <trackbase/LaserCluster.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrDefs.h>
 
#include <g4detectors/PHG4TpcGeom.h>
#include <g4detectors/PHG4TpcGeomContainer.h>
 
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <TPolyLine3D.h>
 
#include <cmath>
#include <iostream>
#include <limits>

namespace
{
  Acts::Vector3 invalid(std::numeric_limits<double>::quiet_NaN(),
                        std::numeric_limits<double>::quiet_NaN(),
                        std::numeric_limits<double>::quiet_NaN());

  std::array<double, 3> invalidArr = {std::numeric_limits<double>::quiet_NaN(),
                        std::numeric_limits<double>::quiet_NaN(),
                        std::numeric_limits<double>::quiet_NaN()};
}

//____________________________________________________________________________
void LaserClusterHelper::loadNodes(PHCompositeNode* topNode)
{
    m_tGeometry = findNode::getClass<ActsGeometry>(topNode,"ActsGeometry");
    if(!m_tGeometry)
    {
        std::cout << "LaserClusterHelper::loadNodes - ActsGeometry not found on node tree" << std::endl;
    }

    m_geom_container = findNode::getClass<PHG4TpcGeomContainer>(topNode, "TPCGEOMCONTAINER");
    if(!m_geom_container)
    {
        std::cout << "LaserClusterHelper::loadNodes - TPCGEOMCONTAINER not found on node tree" << std::endl;
    }
}

//____________________________________________________________________________
Acts::Vector3 LaserClusterHelper::getHitPosition(TrkrDefs::hitsetkey hitsetkey, TrkrDefs::hitkey hitkey) const
{
    //const Acts::Vector3 invalid(std::numeric_limits<double>::quiet_NaN(),
    //                            std::numeric_limits<double>::quiet_NaN(),
    //                            std::numeric_limits<double>::quiet_NaN());

    if(!m_tGeometry || !m_geom_container)
    {
        return invalid;
    }

    const int layer = TrkrDefs::getLayer(hitsetkey);
    const int side = TpcDefs::getSide(hitsetkey);

    PHG4TpcGeom *layer_geom = m_geom_container->GetLayerCellGeom(layer);
    if(!layer_geom)
    {
        return invalid;
    }

    const int iphi = TpcDefs::getPad(hitkey);
    const int it = TpcDefs::getTBin(hitkey);

    const double radius = layer_geom->get_radius();
    const double phi = layer_geom->get_phi(iphi, side);
    
    const double env_x = radius * cos(phi);
    const double env_y = radius * sin(phi);
    double env_z = 0.0;
    //hard code at 0 until better z coordinate calibration is determined
    if(m_useZ)
    {
        double vdrift = m_tGeometry->get_drift_velocity();
        double tdriftmax = layer_geom->get_max_driftlength() / vdrift;

        double zdriftlength = layer_geom->get_zcenter(it) * vdrift;
        // convert z drift length to z position in the TPC
        env_z = tdriftmax * vdrift - zdriftlength;
        if (side == 0)
        {
            env_z = -env_z;
        }
    }

    Acts::Vector3 env_global(env_x, env_y, env_z);
    if(!m_useGlobal)
    {
        return env_global;
    }

    return m_tGeometry->transformTpcEnvelopeToWorld(env_global);
}

//____________________________________________________________________________
Acts::Vector3 LaserClusterHelper::getClusterCentroid(LaserCluster* cluster) const
{
    //const Acts::Vector3 invalid(std::numeric_limits<double>::quiet_NaN(),
    //                            std::numeric_limits<double>::quiet_NaN(),
    //                            std::numeric_limits<double>::quiet_NaN());
    
    if(!cluster)
    {
        return invalid;
    }

    Acts::Vector3 weightedSum(0.0, 0.0, 0.0);
    double adcSum = 0.0;

    const unsigned int nhits = cluster->getNhits();
    for(unsigned int i=0; i<nhits; ++i)
    {
        const LaserClusterHitInfo hit= cluster->getHit(i);
        const Acts::Vector3 hitCoords = getHitPosition(hit.hitsetkey, hit.hitkey);
        if(hitCoords.hasNaN())
        {
            continue;
        }

        weightedSum += hit.adc * hitCoords;
        adcSum += hit.adc;
    }

    if(adcSum <= 0.0)
    {
        return invalid;
    }

    return weightedSum / adcSum;
}

//____________________________________________________________________________
std::array<double, 3> LaserClusterHelper::getClusterHardwareCentroid(LaserCluster* cluster) const
{
    //const Acts::Vector3 invalid(std::numeric_limits<double>::quiet_NaN(),
    //                            std::numeric_limits<double>::quiet_NaN(),
    //                            std::numeric_limits<double>::quiet_NaN());
    
    if(!cluster)
    {
        return invalidArr;
    }

    Acts::Vector3 weightedSum(0.0, 0.0, 0.0);
    double adcSum = 0.0;
    double layerSum = 0.0;
    double iphiSum = 0.0;
    double itSum = 0.0;

    const unsigned int nhits = cluster->getNhits();
    for(unsigned int i=0; i<nhits; ++i)
    {
        const LaserClusterHitInfo hit= cluster->getHit(i);

        const int layer = TrkrDefs::getLayer(hit.hitsetkey);
        const int iphi = TpcDefs::getPad(hit.hitkey);
        const int it = TpcDefs::getTBin(hit.hitkey);

        adcSum += hit.adc;
        layerSum += layer * hit.adc;
        iphiSum += iphi * hit.adc;
        itSum += it * hit.adc;
    }

    if(adcSum <= 0.0)
    {
        return invalidArr;
    }

    return {layerSum / adcSum, iphiSum / adcSum, itSum / adcSum};
}


Acts::Vector3 LaserClusterHelper::getClusterCentroidWithPHGarfield(LaserCluster* cluster) const
{
    if(!cluster)
    {
        return invalid;
    }

    if (cluster->getNhits() < 1)
    {
        return invalid;
    }

    Acts::Vector3 centroid = getClusterCentroid(cluster);
    if (std::isnan(centroid[0]) || std::isnan(centroid[1]) || std::isnan(centroid[2]))
    {
        return invalid;
    }

    // transform back to local if necessary so that z position can be put at readout plane
    if (m_useGlobal)
    {
        centroid = m_tGeometry->transformTpcWorldToEnvelope(centroid);
    }

    const LaserClusterHitInfo hit = cluster->getHit(0);
    const int side = TpcDefs::getSide(hit.hitsetkey);

    Acts::Vector3 readoutPos(centroid[0], centroid[1], side == 1 ? 102 : -102);
    if (m_useGlobal)
    {
        readoutPos = m_tGeometry->transformTpcEnvelopeToWorld(readoutPos);
    }

    std::unique_ptr<TPolyLine3D> path;
    if (!m_useGlobal)
    {
        path.reset(m_phgarfield->ReverseDrift(readoutPos[0], readoutPos[1], readoutPos[2]));
    }
    else
    {
        path.reset(m_phgarfield->ReverseDriftGlobalCoords(readoutPos[0], readoutPos[1], readoutPos[2]));
    }
    
    // just get the last point on the polyline which should be on the central membrane
    if (!path || path->GetN() < 1)
    {
        return invalid;
    }
    const int nPoints = path->GetN();
    const float* points = path->GetP();
    const double xCM = points[3 * (nPoints-1) + 0];
    const double yCM = points[3 * (nPoints-1) + 1];
    const double zCM = points[3 * (nPoints-1) + 2];
    Acts::Vector3 centroidatCM(xCM,yCM,zCM);

    return centroidatCM;
}