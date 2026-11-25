#ifndef KFPARTICLESPHENIX_KFPARTICLEDETECTORINFO_H
#define KFPARTICLESPHENIX_KFPARTICLEDETECTORINFO_H

#include "KFParticle_Config.h"
#include "KFParticle_util.h"

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

  static constexpr int maxTracks{KFParticle_util::maxTracks};

  const KFParticle_Config &m_config;

  unsigned int detector_nHits_MVTX[maxTracks]{0};
  unsigned int detector_nHits_INTT[maxTracks]{0};
  unsigned int detector_nHits_TPC[maxTracks]{0};
  unsigned int detector_nHits_TPOT[maxTracks]{0};
  unsigned int detector_nStates_MVTX[maxTracks]{0};
  unsigned int detector_nStates_INTT[maxTracks]{0};
  unsigned int detector_nStates_TPC[maxTracks]{0};
  unsigned int detector_nStates_TPOT[maxTracks]{0};
  std::vector<float> residual_x[maxTracks];
  std::vector<float> residual_y[maxTracks];
  std::vector<float> residual_z[maxTracks];
  std::vector<int> detector_layer[maxTracks];
  std::vector<int> mvtx_staveID[maxTracks];
  std::vector<int> mvtx_chipID[maxTracks];
  std::vector<int> intt_ladderZID[maxTracks];
  std::vector<int> intt_ladderPhiID[maxTracks];
  std::vector<int> tpc_sectorID[maxTracks];
  std::vector<int> tpc_side[maxTracks];
};

#endif  // KFPARTICLESPHENIX_KFPARTICLEDETECTORINFO_H
