#ifndef KFPARTICLESPHENIX_KFPARTICLEDETECTORINFO_H
#define KFPARTICLESPHENIX_KFPARTICLEDETECTORINFO_H

#include "KFParticle_Config.h"
#include "KFParticle_utilities.h"

#include <string>
#include <vector>

class KFParticle;
class PHCompositeNode;
class TTree;

// Per-daughter tracking-detector information (hits, states, residuals, sensor IDs) for the nTuple
class KFParticle_detectorInfo
{
 public:
  // Binds the configuration only; it is not read until the branches are booked
  explicit KFParticle_detectorInfo(const KFParticle_Config &config)
    : m_config(config)
  {
  }

  void initializeBranches(TTree *m_tree, int daughter_id, const std::string &daughter_number);
  void fill(PHCompositeNode *topNode, const KFParticle &daughter, int daughter_id);
  void clear(int nTracks);

 private:
  void initializeSubDetectorBranches(TTree *m_tree, const std::string &detectorName, int daughter_id, const std::string &daughter_number);

  static constexpr int max_tracks{KFParticle_utilities::max_tracks};

  const KFParticle_Config &m_config;

  unsigned int detector_nHits_MVTX[max_tracks]{0};
  unsigned int detector_nHits_INTT[max_tracks]{0};
  unsigned int detector_nHits_TPC[max_tracks]{0};
  unsigned int detector_nHits_TPOT[max_tracks]{0};
  unsigned int detector_nStates_MVTX[max_tracks]{0};
  unsigned int detector_nStates_INTT[max_tracks]{0};
  unsigned int detector_nStates_TPC[max_tracks]{0};
  unsigned int detector_nStates_TPOT[max_tracks]{0};
  std::vector<float> residual_x[max_tracks];
  std::vector<float> residual_y[max_tracks];
  std::vector<float> residual_z[max_tracks];
  std::vector<int> detector_layer[max_tracks];
  std::vector<int> mvtx_staveID[max_tracks];
  std::vector<int> mvtx_chipID[max_tracks];
  std::vector<int> intt_ladderZID[max_tracks];
  std::vector<int> intt_ladderPhiID[max_tracks];
  std::vector<int> tpc_sectorID[max_tracks];
  std::vector<int> tpc_side[max_tracks];
};

#endif  // KFPARTICLESPHENIX_KFPARTICLEDETECTORINFO_H
