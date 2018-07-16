// This file is really -*- C++ -*-.
#ifndef ClusterIso_h 
#define ClusterIso_h

#include <fun4all/SubsysReco.h>
#include <CLHEP/Vector/ThreeVector.h>
#include <cmath>

class PHCompositeNode;
class Hep3Vector;
class RawTowerGeom;

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
  ClusterIso(const std::string& ,float eTCut, int coneSize);

  virtual int Init(PHCompositeNode*);
  virtual int process_event(PHCompositeNode*);
  virtual int End(PHCompositeNode*);

  void seteTCut(float x);
  void setConeSize(int x);
  const float geteTCut();
  const int getConeSize();
  const CLHEP::Hep3Vector getVertex();

private:
  double getTowerEta(RawTowerGeom* tower_geom, double vx, double vy, double vz); 
  float m_eTCut; ///< cluster must be over this energy set in constructor
  float m_coneSize; ///< delta R around cluster that is considered
  float m_vx; ///< new vertex x value  
  float m_vy; ///< new vertex y value 
  float m_vz; ///< new vertex z value
};



/** \Brief Function to find delta R between 2 objects
 *
 * Takes eta and phi of each object and returns the difference 
 * of the etas and phis added in quadrature. Used to find towers
 * inside a cone of delta R around a cluster.
 */
inline const float deltaR( float eta1, float eta2, float phi1, float phi2 ) {
  float deta = eta1 - eta2;
  float dphi = phi1 - phi2;
  if ( dphi > M_PI ) dphi -= 2 * M_PI ; //corrects to keep range -pi to pi
  if ( dphi < -1*M_PI ) dphi += 2 * M_PI ; //corrects to keep range -pi to pi
  return sqrt( deta*deta + dphi*dphi);
}

#endif
