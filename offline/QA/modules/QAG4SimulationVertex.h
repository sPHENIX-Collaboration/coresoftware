// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef QAG4SIMULATIONVERTEX_H
#define QAG4SIMULATIONVERTEX_H

#include <fun4all/SubsysReco.h>

#include <g4eval/SvtxEvalStack.h>

#include <memory>
#include <set>
#include <string>
#include <utility>

class PHCompositeNode;

class QAG4SimulationVertex : public SubsysReco
{
 public:
  QAG4SimulationVertex(const std::string &name = "QAG4SimulationVertex");

  virtual ~QAG4SimulationVertex() = default;

  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);

  std::string get_histo_prefix();

  void addEmbeddingID(int embeddingID);

 private:
  int load_nodes(PHCompositeNode *);

  PHG4TruthInfoContainer *m_truthContainer = nullptr;

  unsigned int _nlayers_maps = 3;

  std::unique_ptr<SvtxEvalStack> m_svtxEvalStack;

  std::set<int> m_embeddingIDs;
};

#endif  // QAG4SIMULATIONVERTEX_H
