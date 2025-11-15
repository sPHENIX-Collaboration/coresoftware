// This file is really -*- C++ -*-.
#ifndef CLUSTERISO_CLUSTERISO_H
#define CLUSTERISO_CLUSTERISO_H

#include <fun4all/SubsysReco.h>

#include <CLHEP/Vector/ThreeVector.h>

#include <cmath>
#include <string>

class PHCompositeNode;
class RawTowerGeom;
class TowerInfo;

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

  void seteTCut(float eTCut);
  void setConeSize(int coneSize);
  float geteTCut() const;
  //! returns coneSize*10 as an int
  int getConeSize() const;
  CLHEP::Hep3Vector getVertex() const;
  void set_use_towerinfo(bool usetowerinfo)
  {
    m_use_towerinfo = usetowerinfo;
  };

  // Minimum tower energy (GeV) required for a tower to contribute to isolation sum
  void setMinTowerEnergy(float emin) { m_minTowerEnergy = emin; }
  float getMinTowerEnergy() const { return m_minTowerEnergy; }

  void set_cluster_node_name(const std::string& name)
  {
    m_cluster_node_name = name;
  }

 private:
  double getTowerEta(RawTowerGeom* tower_geom, double vx, double vy, double vz);
  bool IsAcceptableTower(TowerInfo* tower);
  float m_eTCut{};     ///< The minimum required transverse energy in a cluster for ClusterIso to be run
  float m_coneSize{};  ///< Size of the cone used to isolate a given cluster
  float m_vx;          ///< Correct vertex x coordinate
  float m_vy;          ///< Correct vertex y coordinate
  float m_vz;          ///< Correct vertex z coordinate
  bool m_do_subtracted;
  bool m_do_unsubtracted;
  bool m_use_towerinfo{true};
  std::string m_cluster_node_name{"CLUSTERINFO_CEMC"};
  float m_minTowerEnergy{-100};  ///< Minimum tower energy for inclusion in isolation calculation
};

#endif
