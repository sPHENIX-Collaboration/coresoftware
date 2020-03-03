#ifndef QA_QAG4SimulationTracking_H
#define QA_QAG4SimulationTracking_H

#include <fun4all/SubsysReco.h>

#include <memory>
#include <string>
#include <set>
#include <utility>

#if !defined(__CINT__) || defined(__CLING__)
#include <cstdint>
#else
#include <stdint.h>
#endif

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
  QAG4SimulationTracking();
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

  void setEtaRange(double low, double high)
  {
    m_etaRange.first = low;
    m_etaRange.second = high;
  }

 private:
#if !defined(__CINT__) || defined(__CLING__)
  //CINT is not c++11 compatible
  std::shared_ptr<SvtxEvalStack> _svtxEvalStack;
  std::set<int> m_embeddingIDs;
#endif
  std::pair<double, double> m_etaRange;

  PHG4TruthInfoContainer *_truthContainer;
};

#endif  // QA_QAG4SimulationTracking_H
