#ifndef SVTXTRUTHRECOTABLEEVAL_H
#define SVTXTRUTHRECOTABLEEVAL_H

#include <fun4all/SubsysReco.h>

#include <trackbase_historic/PHG4ParticleSvtxMap.h>
#include <trackbase_historic/SvtxPHG4ParticleMap.h>

#include <string>
#include <memory>

class PHCompositeNode;
class PHG4TruthInfoContainer;
class SvtxEvalStack;

class SvtxTruthRecoTableEval : public SubsysReco
{
 public:
  SvtxTruthRecoTableEval(const std::string &name = "SvtxTruthRecoTableEval");

  virtual ~SvtxTruthRecoTableEval();

  int Init(PHCompositeNode *) override;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int ResetEvent(PHCompositeNode *topNode) override;

  int End(PHCompositeNode *topNode) override;

  //! minimal momentum of particle to trigger an analysis and entry into the truth map.
  //! default to 50MeV, bending before reaching the active TPC
  void setMinMomentumTruthMap(const double v) { m_minMomentumTruthMap = v; }

 private:
  int createNodes(PHCompositeNode *topNode);

  void fillTruthMap(PHCompositeNode *topNode);
  void fillRecoMap(PHCompositeNode *topNode);

  bool m_scanForPrimaries = false;

  //! minimal momentum of particle to trigger an analysis and entry into the truth map.
  //! default to 50MeV, bending before reaching the active TPC
  double m_minMomentumTruthMap = 50e-3;

  PHG4ParticleSvtxMap *m_truthMap = nullptr;
  SvtxPHG4ParticleMap *m_recoMap = nullptr;
  std::unique_ptr<SvtxEvalStack> m_svtxevalstack;
};

#endif  // SVTXTRUTHRECOTABLEEVAL_H
