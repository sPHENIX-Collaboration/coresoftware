#ifndef KFPARTICLESPHENIX_KFPARTICLETRUTHINFO_H
#define KFPARTICLESPHENIX_KFPARTICLETRUTHINFO_H

#include "KFParticle_Config.h"
#include "KFParticle_utilities.h"

#include <limits>
#include <string>
#include <vector>

class GlobalVertex;
class GlobalVertexMap;
class KFParticle;
class MbdVertexMap;
class PHCompositeNode;
class PHG4Particle;
class PHG4TruthInfoContainer;
class PHG4VtxPoint;
class PHHepMCGenEvent;
class PHHepMCGenEventMap;
class SvtxEvalStack;
class SvtxTrack;
class SvtxTrackEval;
class SvtxTruthEval;
class SvtxVertexEval;
class SvtxVertexMap;
class TTree;

namespace HepMC
{
  class GenParticle;
}

// Per-daughter truth matching (true kinematics, production vertex, true DCA, decay history) for the nTuple
class KFParticle_truthInfo
{
 public:
  // Binds the configuration only; it is not read until the branches are filled
  explicit KFParticle_truthInfo(const KFParticle_Config &config)
    : m_config(config)
  {
  }

  void initializeBranches(TTree *m_tree, int daughter_id, const std::string &daughter_number, bool m_constrain_to_vertex_truthMatch);
  void fill(PHCompositeNode *topNode, const KFParticle &daughter, int daughter_id, const KFParticle &kfvertex, bool m_constrain_to_vertex_truthMatch);
  int getHepMCInfo(PHCompositeNode *topNode, const KFParticle &daughter, int daughter_id);
  void clear(int nTracks);

 private:
  static GlobalVertex *getVertex(unsigned int vertex_id, GlobalVertexMap *vertexmap);
  PHG4Particle *getTruthTrack(SvtxTrack *thisTrack, PHCompositeNode *topNode);
  void fillGeant4Branch(PHG4Particle *particle, int daughter_id);
  void fillHepMCBranch(HepMC::GenParticle *particle, int daughter_id);

  static constexpr int max_tracks{KFParticle_utilities::max_tracks};

  const KFParticle_Config &m_config;

  SvtxEvalStack *m_svtx_evalstack{nullptr};
  SvtxTrackEval *trackeval{nullptr};
  SvtxTruthEval *trutheval{nullptr};
  SvtxVertexEval *vertexeval{nullptr};

  PHG4Particle *g4particle{nullptr};
  PHG4VtxPoint *g4vertex_point{nullptr};

  MbdVertexMap *dst_mbdvertexmap{nullptr};
  SvtxVertexMap *dst_vertexmap{nullptr};

  PHG4TruthInfoContainer *m_truthinfo{nullptr};
  PHHepMCGenEventMap *m_geneventmap{nullptr};
  PHHepMCGenEvent *m_genevt{nullptr};

  float m_true_daughter_vertex_x[max_tracks]{std::numeric_limits<float>::quiet_NaN()};
  float m_true_daughter_vertex_y[max_tracks]{std::numeric_limits<float>::quiet_NaN()};
  float m_true_daughter_vertex_z[max_tracks]{std::numeric_limits<float>::quiet_NaN()};
  float m_true_daughter_ip[max_tracks]{std::numeric_limits<float>::quiet_NaN()};
  float m_true_daughter_ip_xy[max_tracks]{std::numeric_limits<float>::quiet_NaN()};
  float m_true_daughter_px[max_tracks]{std::numeric_limits<float>::quiet_NaN()};
  float m_true_daughter_py[max_tracks]{std::numeric_limits<float>::quiet_NaN()};
  float m_true_daughter_pz[max_tracks]{std::numeric_limits<float>::quiet_NaN()};
  float m_true_daughter_p[max_tracks]{std::numeric_limits<float>::quiet_NaN()};
  float m_true_daughter_pt[max_tracks]{std::numeric_limits<float>::quiet_NaN()};
  int m_true_daughter_id[max_tracks]{std::numeric_limits<int>::quiet_NaN()};
  float m_true_daughter_pv_x[max_tracks]{std::numeric_limits<float>::quiet_NaN()};
  float m_true_daughter_pv_y[max_tracks]{std::numeric_limits<float>::quiet_NaN()};
  float m_true_daughter_pv_z[max_tracks]{std::numeric_limits<float>::quiet_NaN()};

  std::vector<int> m_true_daughter_track_history_PDG_ID[max_tracks];
  std::vector<float> m_true_daughter_track_history_PDG_mass[max_tracks];
  std::vector<float> m_true_daughter_track_history_px[max_tracks];
  std::vector<float> m_true_daughter_track_history_py[max_tracks];
  std::vector<float> m_true_daughter_track_history_pz[max_tracks];
  std::vector<float> m_true_daughter_track_history_pE[max_tracks];
  std::vector<float> m_true_daughter_track_history_pT[max_tracks];
};

#endif  // KFPARTICLESPHENIX_KFPARTICLETRUTHINFO_H
