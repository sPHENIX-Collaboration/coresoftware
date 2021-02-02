#ifndef QA_QAG4SimulationUpsilon_H
#define QA_QAG4SimulationUpsilon_H

#include <fun4all/SubsysReco.h>

#include <memory>
#include <set>
#include <string>
#include <utility>

class PHCompositeNode;
class PHG4TruthInfoContainer;
class SvtxEvalStack;

/// \class QAG4SimulationUpsilon
class QAG4SimulationUpsilon : public SubsysReco
{
 public:
  QAG4SimulationUpsilon(const std::string &name = "QAG4SimulationUpsilon");
  virtual ~QAG4SimulationUpsilon() {}

  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);

  // common prefix for QA histograms
  std::string get_histo_prefix();

  //! If added, only process truth particle associated with the selected list of EmbeddingIDs
  //! Call multiple times to add multiple EmbeddingIDs
  //! For EmbeddingID<0, all negative embedding IDs are accepted for pile up events.
  void addEmbeddingID(int embeddingID);

  void setEtaRange(double low, double high)
  {
    m_etaRange.first = low;
    m_etaRange.second = high;
  }

  void setQuarkoniaPID(const int pid)
  {
    m_quarkoniaPID = pid;
  }

  void setDaughterAbsPID(const int pid)
  {
    m_daughterAbsPID = pid;
  }

 private:
  std::shared_ptr<SvtxEvalStack> _svtxEvalStack;
  std::set<int> m_embeddingIDs;
  std::pair<double, double> m_etaRange;

  PHG4TruthInfoContainer *_truthContainer;

  int m_quarkoniaPID = 553;
  int m_daughterAbsPID = 11;
};

#endif  // QA_QAG4SimulationUpsilon_H
