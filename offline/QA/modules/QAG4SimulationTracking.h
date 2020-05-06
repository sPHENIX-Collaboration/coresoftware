#ifndef QA_QAG4SimulationTracking_H
#define QA_QAG4SimulationTracking_H

#include <fun4all/SubsysReco.h>

#include <memory>
#include <set>
#include <string>
#include <utility>

class PHCompositeNode;
class PHG4TruthInfoContainer;
class PHG4Particle;
class CaloEvalStack;
class SvtxEvalStack;
class SvtxTrack;

/// \class QAG4SimulationTracking
class QAG4SimulationTracking : public SubsysReco
{
 public:
  QAG4SimulationTracking(const std::string &name = "QAG4SimulationTracking");
  virtual ~QAG4SimulationTracking() {}

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

 private:

  std::unique_ptr<SvtxEvalStack> _svtxEvalStack;
  std::set<int> m_embeddingIDs;

  //! range of the truth track eta to be analyzed
  std::pair<double, double> m_etaRange = {-1, 1};

  //! only count unique truth<->reco track pair in tracking efficiency
  bool m_uniqueTrackingMatch = true;

  PHG4TruthInfoContainer *_truthContainer = nullptr;
};

#endif  // QA_QAG4SimulationTracking_H
