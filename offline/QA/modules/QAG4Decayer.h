#include <fun4all/SubsysReco.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>

#include <g4eval/SvtxEvalStack.h>

#include <trackbase_historic/ActsTransformations.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxVertex.h>
#include <trackbase_historic/SvtxVertexMap.h>
#include <trackbase_historic/TrackSeed.h>

#include <TBranch.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>

#include <string>

class PHCompositeNode;
class PHG4TruthInfoContainer;
class PHG4Particle;
class PHG4VtxPoint;

class QAG4Decayer : public SubsysReco
{
 public:
  //		QAG4Decayer();
  QAG4Decayer(const std::string &name = "QAG4Decayer");

  virtual ~QAG4Decayer();

  int Init(PHCompositeNode *topNode) override;
  //		int Init() override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;
  std::vector<int> Channel(int pdgid, std::vector<int> Daughter);

 private:
  PHG4TruthInfoContainer *m_truth_info;
  //  PHG4Particle *m_g4particle;
  SvtxEvalStack *_svtxevalstack = NULL;

  TH1F *QAPx;
  TH1F *QAPy;
  TH1F *QAPz;
  TH1F *QAE;
  TH1F *TotalVtxStat;
  TH1F *BadVtxStat;
  TH1F *BadVtxPercent;
  TH1F *HFHadronStat;

  TH2F *DecayBR;

  TH1F *SVtoPVDistance;
  TH1F *ProperLifeTime;
  TH1F *QACosTheta;
  TH1F *OtherHis;
  TFile *fout;
  TTree *QATree;
  int EvtID;
  float LifeTime;
  TH1D *MassHis;
  TH1D *NPartHis;
  TH1D *BR1DHis;
  int NParticles;
  float D0Px;
  float D0Py;
  float D0Pz;
  float D0E;
  std::vector<int> PVDaughtersPDGID;

  bool m_write_nTuple;
  bool m_write_QAHists;
  bool m_SaveFiles;
};
