// This file is really -*- C++ -*-.
#ifndef CLUSTERISO_CLUSTERISO_H
#define CLUSTERISO_CLUSTERISO_H

#include <fun4all/SubsysReco.h>

#include <CLHEP/Vector/ThreeVector.h>

#include <cmath>
#include <string>

class PHCompositeNode;
class RawTowerGeom;

/** \Brief Tool to find isolation energy of each EMCal cluster.
 * 
 * This tool finds isoET of clusters by summing towers energy
 * in a cone of radius R around the cluster and subtracting 
 * the cluster from the sum
 */

class ClusterIso : public SubsysReco
{
 public:
  /**
   * Constructor for ClusterIso Class
   * the coneSize is taken in as an integer multiple of .1 ie if you want R=.2 pass 2
   */
  ClusterIso(const std::string&, float eTCut, int coneSize, bool do_subtracted, bool do_unsubtracted);

  int Init(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;
  int End(PHCompositeNode*) override;

  void seteTCut(float x);
  void setConeSize(int x);
  /*const*/ float geteTCut();
  //! returns coneSize*10 as an int
  /*const*/ int getConeSize();
  /*const*/ CLHEP::Hep3Vector getVertex();

 private:
  double getTowerEta(RawTowerGeom* tower_geom, double vx, double vy, double vz);
  float m_eTCut;     ///< The minimum required transverse energy in a cluster for ClusterIso to be run
  float m_coneSize;  ///< Size of the cone used to isolate a given cluster
  float m_vx;        ///< Correct vertex x coordinate
  float m_vy;        ///< Correct vertex y coordinate
  float m_vz;        ///< Correct vertex z coordinate
  bool m_do_subtracted;
  bool m_do_unsubtracted;
};

/** \Brief Function to find delta R between 2 objects
 *
 * Takes the eta and phi of each object and returns the difference 
 * of the etas and phis added in quadrature. Used to find towers
 * inside a cone of delta R around a cluster.
 */
inline /*const*/ float deltaR(float eta1, float eta2, float phi1, float phi2)
{
  float deta = eta1 - eta2;
  float dphi = phi1 - phi2;
  if (dphi > M_PI) dphi -= 2 * M_PI;       //corrects to keep range -pi to pi
  if (dphi < -1 * M_PI) dphi += 2 * M_PI;  //corrects to keep range -pi to pi
  return sqrt(deta * deta + dphi * dphi);
}

#endif
