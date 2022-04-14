#ifndef SVTXTRUTHRECOTABLEEVAL_H
#define SVTXTRUTHRECOTABLEEVAL_H

#include <fun4all/SubsysReco.h>

#include <trackbase_historic/PHG4ParticleSvtxMap.h>
#include <trackbase_historic/SvtxPHG4ParticleMap.h>

#include <string>

class PHCompositeNode;
class PHG4TruthInfoContainer;
class SvtxEvalStack;

class SvtxTruthRecoTableEval : public SubsysReco
{
 public:

  SvtxTruthRecoTableEval(const std::string &name = "SvtxTruthRecoTableEval");

  virtual ~SvtxTruthRecoTableEval();

  int Init(PHCompositeNode*) override;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int ResetEvent(PHCompositeNode *topNode) override;
 
  int End(PHCompositeNode *topNode) override;


 private:
  int createNodes(PHCompositeNode* topNode);

  void fillTruthMap(PHCompositeNode *topNode);
  void fillRecoMap(PHCompositeNode *topNode);

  bool m_scanForPrimaries = false;

  PHG4ParticleSvtxMap *m_truthMap = nullptr;
  SvtxPHG4ParticleMap *m_recoMap = nullptr;
  SvtxEvalStack *m_svtxevalstack = nullptr;
};

#endif // SVTXTRUTHRECOTABLEEVAL_H
