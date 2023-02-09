#ifndef QA_QAG4SimulationTracking_H
#define QA_QAG4SimulationTracking_H

#include <g4eval/SvtxEvalStack.h>

#include <trackbase/TrkrDefs.h>  // for cluskey

#include <fun4all/SubsysReco.h>

#include <memory>
#include <set>
#include <string>
#include <utility>

class PHCompositeNode;
class SvtxTrackMap;
class PHG4Hit;
class PHG4HitContainer;
class PHG4TruthInfoContainer;
class TrkrClusterContainer;
class TrkrClusterHitAssoc;
class TrkrHitTruthAssoc;
class SvtxVertexMap;

/// \class QAG4SimulationTracking
class QAG4SimulationTracking : public SubsysReco
{
 public:
  QAG4SimulationTracking(const std::string &name = "QAG4SimulationTracking");
  virtual ~QAG4SimulationTracking() = default;

  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);

  // common prefix for QA histograms
  std::string get_histo_prefix();

  //! If added, only process truth particle associated with the selected list of EmbeddingIDs
  //! Call multiple times to add multiple EmbeddingIDs
  //! For EmbeddingID<0, all negative embedding IDs are accepted for pile up events.
  void addEmbeddingID(int embeddingID);

  //! range of the truth track eta to be analyzed
  void setEtaRange(double low, double high)
  {
    m_etaRange.first = low;
    m_etaRange.second = high;
  }

  //! only count unique truth<->reco track pair in tracking efficiency
  void setUniqueTrackingMatch(bool b)
  {
    m_uniqueTrackingMatch = b;
  }

  void set_embed_id_cut(const int id) { m_embed_id_cut = id; }

 private:
  /// load nodes
  int load_nodes(PHCompositeNode *);

  void get_dca(SvtxTrack *track, float &dca3dxy, float &dca3dz,
               float &dca3dxysigma, float &dca3dzsigma);
  // get geant hits associated to a cluster
  using G4HitSet = std::set<PHG4Hit *>;
  G4HitSet find_g4hits(TrkrDefs::cluskey) const;

  std::unique_ptr<SvtxEvalStack> m_svtxEvalStack;
  std::set<int> m_embeddingIDs;

  //! range of the truth track eta to be analyzed
  std::pair<double, double> m_etaRange = {-1, 1};

  //! only count unique truth<->reco track pair in tracking efficiency
  bool m_uniqueTrackingMatch = true;

  //! cut for selecting on foreground
  int m_embed_id_cut = 0;

  PHG4TruthInfoContainer *m_truthContainer = nullptr;
  SvtxTrackMap *m_trackMap = nullptr;
  SvtxVertexMap *m_vertexMap = nullptr;

  TrkrClusterContainer *m_cluster_map = nullptr;
  TrkrClusterHitAssoc *m_cluster_hit_map = nullptr;
  TrkrHitTruthAssoc *m_hit_truth_map = nullptr;

  PHG4HitContainer *m_g4hits_tpc = nullptr;
  PHG4HitContainer *m_g4hits_intt = nullptr;
  PHG4HitContainer *m_g4hits_mvtx = nullptr;
  PHG4HitContainer *m_g4hits_micromegas = nullptr;
};

#endif  // QA_QAG4SimulationTracking_H
