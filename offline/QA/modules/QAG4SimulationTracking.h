#ifndef QA_QAG4SimulationTracking_H
#define QA_QAG4SimulationTracking_H

#include <fun4all/SubsysReco.h>

#include <memory>
#include <string>

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

  virtual ~QAG4SimulationTracking(){}

  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);


  //! common prefix for QA histograms
  std::string
  get_histo_prefix();

 private:

#if !defined(__CINT__) || defined(__CLING__)
  //CINT is not c++11 compatible
  std::shared_ptr<SvtxEvalStack> m_svtxEvalStack;
#endif

  PHG4TruthInfoContainer *m_truthContainer;


};

#endif  // QA_QAG4SimulationTracking_H
