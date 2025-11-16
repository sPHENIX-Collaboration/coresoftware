#include <fun4all/SubsysReco.h>

#include <string>
#include <vector>

class PHCompositeNode;
class PHG4TruthInfoContainer;
class PHG4Particle;
class PHG4VtxPoint;
class TFile;
class TTree;

class QAG4Decayer : public SubsysReco
{
 public:
  //		QAG4Decayer();
  QAG4Decayer(const std::string &name = "QAG4Decayer");

  ~QAG4Decayer() override = default;

  int Init(PHCompositeNode *topNode) override;
  //		int Init() override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;
  std::vector<int> Channel(int pdgid, std::vector<int> Daughter);

 private:
  PHG4TruthInfoContainer *m_truth_info{nullptr};

  TFile *fout{nullptr};
  TTree *QATree{nullptr};
  int EvtID{-1};
  float LifeTime{-1};
  int NParticles{-1};
  std::vector<int> PVDaughtersPDGID;

  bool m_write_nTuple;
  bool m_SaveFiles;
};
