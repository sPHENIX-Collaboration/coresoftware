#include "LaserClusterHelper.h"
 
#include <trackbase/LaserCluster.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrDefs.h>
 
#include <g4detectors/PHG4TpcGeom.h>
#include <g4detectors/PHG4TpcGeomContainer.h>
 
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
 
#include <cmath>
#include <iostream>
#include <limits>

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
Acts::Vector3 LaserClusterHelper::getHitGlobalPosition(TrkrDefs::hitsetkey hitsetkey, TrkrDefs::hitkey hitkey) const
{
    const Acts::Vector3 invalid(std::numeric_limits<double>::quiet_NaN(),
                                std::numeric_limits<double>::quiet_NaN(),
                                std::numeric_limits<double>::quiet_NaN());

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
    return m_tGeometry->transformTpcEnvelopeToWorld(env_global);
}

//____________________________________________________________________________
Acts::Vector3 LaserClusterHelper::getClusterCentroid(LaserCluster* cluster) const
{
    const Acts::Vector3 invalid(std::numeric_limits<double>::quiet_NaN(),
                                std::numeric_limits<double>::quiet_NaN(),
                                std::numeric_limits<double>::quiet_NaN());
    
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
        const Acts::Vector3 global = getHitGlobalPosition(hit.hitsetkey, hit.hitkey);
        if(global.hasNaN())
        {
            continue;
        }

        weightedSum += hit.adc * global;
        adcSum += hit.adc;
    }

    if(adcSum <= 0.0)
    {
        return invalid;
    }

    return weightedSum / adcSum;
}