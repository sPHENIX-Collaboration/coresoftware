// This file is really -*- C++ -*-.
#ifndef ClusterIso_h 
#define ClusterIso_h

#include <fun4all/SubsysReco.h>

#include <phool/getClass.h>         
#include <phool/PHCompositeNode.h>  

#include <calobase/RawCluster.h>          
#include <calobase/RawClusterUtility.h>   
#include <calobase/RawClusterContainer.h> 

#include <calobase/RawTower.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerContainer.h> 
#include <calobase/RawTowerGeomContainer.h>  

#include <g4vertex/GlobalVertex.h>    
#include <g4vertex/GlobalVertexMap.h> 

#include <cmath>
class PHCompositeNode;

/** \Brief Tool to find isolation energy of each EMCal cluster.
 * 
 * This tool finds isoET of clusters by summing towers energy
 * in a cone of radius R around the cluster and subtracting 
 * the cluster from the sum
 */

class ClusterIso: public SubsysReco
{
public:
/**
 * Constructor for ClusterIso Class
 */
ClusterIso(const std::string& ,float eTCut, float coneSize);

virtual int Init(PHCompositeNode*);
virtual int process_event(PHCompositeNode*);
virtual int End(PHCompositeNode*);

void seteTCut(float x);
void setConeSize(float x);
const float geteTCut();
const float getConeSize();
const CLHEP::Hep3Vector getVertex();


private:
  float m_eTCut; ///< cluster must be over this energy set in constructor
  float m_coneSize; ///< delta R around cluster that is considered
  float m_vx; ///< new vertex x value  
  float m_vy; ///< new vertex y value 
  float m_vz; ///< new vertex z value
};

/** \Brief Function to get correct tower eta
 *
 * Each tower is calculated using the vertex (0,0,0)
 * which is incorrect in many collisions. This function 
 * uses geometry to find eta using correct vertex.
 */
inline double getTowerEta(RawTowerGeom* tower_geom, double vx, double vy, double vz) 
  {
    if(vx==0&&vy==0&&vz==0){
      return tower_geom->get_eta();
    }
    else{
     double r= tower_geom->get_center_radius();
     double x = r*cos(tower_geom->get_phi())-vx;
     double y = r*sin(tower_geom->get_phi())-vy;
     double z = r/tan(2*atan2(exp(-1*tower_geom->get_eta()),1))-vz;
     r= sqrt(x*x+y*y);
     return -log(tan(atan2(r,z)/2.));
    }
}

/** \Brief Function to find delta R between 2 objects
 *
 * Takes eta and phi of each object and returns the difference 
 * of the etas and phis added in quadrature. Used to find towers
 * inside a cone of delta R around a cluster.
 */
inline const float deltaR( float eta1, float eta2, float phi1, float phi2 ) {

    float deta = eta1 - eta2;
    float dphi = phi1 - phi2;
    if ( dphi > 3.14159 ) dphi -= 2 * 3.14159; //corrects to keep range -pi to pi
    if ( dphi < -3.14159 ) dphi += 2 * 3.14159; //corrects to keep range -pi to pi

    return sqrt( deta*deta + dphi*dphi);
}

#endif
