#ifndef CALOVTXRECO_CALOVTXRECO_H
#define CALOVTXRECO_CALOVTXRECO_H

#include "CaloVtxAlgo.h"
#include <memory>
#include <vector>
#include <fun4all/SubsysReco.h>

class CaloVertexMap;
class PHCompositeNode;

class CaloVtxReco : public SubsysReco
{
 public:
  enum ALGOTYPE
    {
      JETSKEW = 0,
      CALOZ = 1,
      JET_MLP = 2
    };
  
  CaloVtxReco(const std::string &name = "CaloVtxReco");

  virtual ~CaloVtxReco() = default;

  int createNodes(PHCompositeNode *topNode);

  int InitRun(PHCompositeNode *topNode) override;
  
  int process_event(PHCompositeNode *topNode) override;
  
  void registerAlgo(CaloVtxAlgo* algo) { m_algos.push_back(std::move(algo)); }
  
 private:
  std::vector<CaloVtxAlgo*> m_algos{};
  
  CaloVertexMap *m_calovtxmap{nullptr};
};

#endif  // CALOVTXRECO_CALOVTXRECO_H
