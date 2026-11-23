#ifndef KFPARTICLESPHENIX_KFPARTICLECALOINFO_H
#define KFPARTICLESPHENIX_KFPARTICLECALOINFO_H

#include "KFParticle_Config.h"
#include "KFParticle_utilities.h"

// Cluster
#include <calobase/RawClusterDefs.h>

#include <cmath>
#include <limits>
#include <string>
#include <vector>

class KFParticle;
class PHCompositeNode;
class RawClusterContainer;
class RawTowerContainer;
class RawTowerGeomContainer;
class TowerInfoContainer;
class TTree;

// Per-daughter track-to-calorimeter matching (EMCal live; IHCal/OHCal branches are NaN for now) for the nTuple
class KFParticle_caloInfo
{
 public:
  // Binds the configuration only; it is not read until the branches are booked
  explicit KFParticle_caloInfo(const KFParticle_Config &config)
    : m_config(config)
  {
  }

  void initializeBranches(TTree *m_tree, int daughter_id, const std::string &daughter_number);
  void fill(PHCompositeNode *topNode, const KFParticle &daughter, int daughter_id, bool &isTrackEMCalmatch, const KFParticle &vertex_in);
  void clear(int nTracks);

 private:
  void Get5x5CellInfo(
    RawClusterDefs::keytype key_in,
    float &E3x3,
    float &E5x5,
    std::vector<unsigned int> &v_emcal_5x5Cell_Phi,
    std::vector<unsigned int> &v_emcal_5x5Cell_Eta,
    std::vector<float> &v_emcal_5x5Cell_E
  );

  // Wrap an azimuthal difference into [-pi, pi]
  static float PiRange(float deltaPhi) { return std::remainder(deltaPhi, 2 * M_PI); }

  static constexpr int max_tracks{KFParticle_utilities::max_tracks};

  const KFParticle_Config &m_config;

  // Only used by the commented-out IHCal/OHCal matching
  float m_ihcal_radius_user{117};
  float m_ohcal_radius_user{177.423};
  float m_ihcal_e_low_cut{0.01};
  float m_ohcal_e_low_cut{0.01};

  // Vectors for calo info
  float detector_emcal_deltaphi[max_tracks]{std::numeric_limits<float>::quiet_NaN()};
  float detector_emcal_deltaeta[max_tracks]{std::numeric_limits<float>::quiet_NaN()};
  float detector_emcal_deltaz[max_tracks]{std::numeric_limits<float>::quiet_NaN()};
  float detector_emcal_energy_3x3[max_tracks]{std::numeric_limits<float>::quiet_NaN()};
  float detector_emcal_energy_5x5[max_tracks]{std::numeric_limits<float>::quiet_NaN()};
  float detector_emcal_cluster_energy[max_tracks]{std::numeric_limits<float>::quiet_NaN()};
  float detector_emcal_cluster_prob[max_tracks]{std::numeric_limits<float>::quiet_NaN()};

  float detector_ihcal_deltaphi[max_tracks]{std::numeric_limits<float>::quiet_NaN()};
  float detector_ihcal_deltaeta[max_tracks]{std::numeric_limits<float>::quiet_NaN()};
  float detector_ihcal_energy_3x3[max_tracks]{std::numeric_limits<float>::quiet_NaN()};
  float detector_ihcal_energy_5x5[max_tracks]{std::numeric_limits<float>::quiet_NaN()};
  float detector_ihcal_cluster_energy[max_tracks]{std::numeric_limits<float>::quiet_NaN()};
  float detector_ohcal_deltaphi[max_tracks]{std::numeric_limits<float>::quiet_NaN()};
  float detector_ohcal_deltaeta[max_tracks]{std::numeric_limits<float>::quiet_NaN()};
  float detector_ohcal_energy_3x3[max_tracks]{std::numeric_limits<float>::quiet_NaN()};
  float detector_ohcal_energy_5x5[max_tracks]{std::numeric_limits<float>::quiet_NaN()};
  float detector_ohcal_cluster_energy[max_tracks]{std::numeric_limits<float>::quiet_NaN()};
  // Detailed Calo Info
  float detector_emcal_projection_z[max_tracks]{std::numeric_limits<float>::quiet_NaN()};
  float detector_emcal_projection_eta[max_tracks]{std::numeric_limits<float>::quiet_NaN()};
  float detector_emcal_projection_phi[max_tracks]{std::numeric_limits<float>::quiet_NaN()};
  float detector_emcal_cluster_z[max_tracks]{std::numeric_limits<float>::quiet_NaN()};
  float detector_emcal_cluster_eta[max_tracks]{std::numeric_limits<float>::quiet_NaN()};
  float detector_emcal_cluster_phi[max_tracks]{std::numeric_limits<float>::quiet_NaN()};
  float detector_emcal_chi2[max_tracks]{std::numeric_limits<float>::quiet_NaN()};
  unsigned int detector_emcal_ntowers[max_tracks]{std::numeric_limits<unsigned int>::quiet_NaN()};
  // 5x5 grid
  std::vector<unsigned int> detector_emcal_5x5Cell_Phi[max_tracks];
  std::vector<unsigned int> detector_emcal_5x5Cell_Eta[max_tracks];
  std::vector<float> detector_emcal_5x5Cell_E[max_tracks];

  // Containers for calo info
  RawTowerGeomContainer *EMCalGeo{nullptr};
  RawClusterContainer *clustersEM{nullptr};
  TowerInfoContainer *_towersEM{nullptr};
  RawTowerGeomContainer *IHCalGeo{nullptr};
  RawClusterContainer *clustersIH{nullptr};
  RawTowerContainer *_towersIH{nullptr};
  RawTowerGeomContainer *OHCalGeo{nullptr};
  RawClusterContainer *clustersOH{nullptr};
  RawTowerContainer *_towersOH{nullptr};
};

#endif  // KFPARTICLESPHENIX_KFPARTICLECALOINFO_H
