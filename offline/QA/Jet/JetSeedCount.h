// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef JETSEEDCOUNT_H
#define JETSEEDCOUNT_H

#include <qautils/QAHistManagerDef.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/SubsysReco.h>

#include <fastjet/PseudoJet.hh>

#include <limits>
#include <string>
#include <vector>

class PHCompositeNode;

class JetSeedCount : public SubsysReco
{
 public:
  JetSeedCount(const std::string &recojetname = "AntiKt_Tower_r04",
               const std::string &truthjetname = "AntiKt_Truth_r04",
               const std::string &outputfilename = "myjetanalysis.root");

  ~JetSeedCount() override = default;

  void
  setEtaRange(double low, double high)
  {
    m_etaRange.first = low;
    m_etaRange.second = high;
  }
  void
  setPtRange(double low, double high)
  {
    m_ptRange.first = low;
    m_ptRange.second = high;
  }

  int Init(PHCompositeNode *topNode) override;

  int InitRun(PHCompositeNode *topNode) override;

  int process_event(PHCompositeNode *topNode) override;

  /// Called at the end of all processing.
  int End(PHCompositeNode *topNode) override;

  void Print(const std::string &what = "ALL") const override;

 private:
  Fun4AllHistoManager *m_manager{nullptr};

  int m_event{0};
  int m_seed_sub{std::numeric_limits<int>::max()};
  int m_seed_raw{std::numeric_limits<int>::max()};
  double z_vtx{std::numeric_limits<double>::quiet_NaN()};

  std::string m_recoJetName;
  std::string m_truthJetName;
  std::string m_outputFileName;

  std::pair<double, double> m_etaRange;
  std::pair<double, double> m_ptRange;

  std::vector<double> m_RawEta;
  std::vector<double> m_RawPhi;
  std::vector<double> m_SubEta;
  std::vector<double> m_SubPhi;
  std::vector<int> m_centrality;
  std::vector<int> m_centrality_diff;

  std::vector<int> m_raw_counts;
  std::vector<int> m_sub_counts;
  std::vector<double> m_rawpt;
  std::vector<double> m_subpt;
  std::vector<double> m_rawpt_all;
  std::vector<double> m_subpt_all;
  std::vector<double> m_rawenergy;
  std::vector<double> m_subenergy;
  std::vector<double> m_rawcent;
  std::vector<double> m_subcent;
};

#endif  // JETSEEDCOUNT_H
