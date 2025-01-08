// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef JETSEEDCOUNT_H
#define JETSEEDCOUNT_H

#include "JetQADefs.h"

#include <qautils/QAHistManagerDef.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/SubsysReco.h>

#include <fastjet/PseudoJet.hh>

#include <limits>
#include <string>
#include <vector>


class PHCompositeNode;
class TH1;
class TH2;

class JetSeedCount : public SubsysReco
{
 public:
  JetSeedCount(const std::string &moduleName = "JetSeedCount",
               const std::string &recojetname = "AntiKt_Tower_r04",
               const std::string &rawSeedName = "AntiKt_TowerInfo_HIRecoSeedsRaw_r02",
               const std::string &subSeedName = "AntiKt_TowerInfo_HIRecoSeedsSub_r02",
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
  void
  setWriteToOutputFile(bool write)
  {
    m_writeToOutputFile = write;
  }
  void
  setHistTag(const std::string &tag)
  {
    m_histTag = tag;
  }
  void
  setTrgToSelect(const uint32_t trig = JetQADefs::GL1::MBDNSJet1)
  {
    m_doTrgSelect = true;
    m_trgToSelect = trig;
  }
  void
  setPPMode(const bool pp)
  {
    m_inPPMode = pp;
  }

  int Init(PHCompositeNode *topNode) override;

  int InitRun(PHCompositeNode *topNode) override;

  int process_event(PHCompositeNode *topNode) override;

  /// Called at the end of all processing.
  int End(PHCompositeNode *topNode) override;

  void Print(const std::string &what = "ALL") const override;

 private:
  Fun4AllHistoManager *m_manager{nullptr};

  bool m_writeToOutputFile{false};
  bool m_inPPMode{false};

  int m_event{0};

  std::string m_moduleName;
  std::string m_recoJetName;
  std::string m_rawSeedName;
  std::string m_subSeedName;
  std::string m_truthJetName;
  std::string m_outputFileName;
  std::string m_histTag;

  std::pair<double, double> m_etaRange;
  std::pair<double, double> m_ptRange;

  // trigger selection
  bool m_doTrgSelect {false};
  uint32_t m_trgToSelect {0};

  TH1* m_hRawSeedCount{nullptr};
  TH1* m_hRawPt{nullptr};
  TH1* m_hRawPt_All{nullptr};
  TH2* m_hRawEtaVsPhi{nullptr};
  TH1* m_hSubSeedCount{nullptr};
  TH1* m_hSubPt{nullptr};
  TH1* m_hSubPt_All{nullptr};
  TH2* m_hSubEtaVsPhi{nullptr};
  TH2* m_hRawSeedEnergyVsCent{nullptr};
  TH2* m_hSubSeedEnergyVsCent{nullptr};
  TH1* m_hCentMbd{nullptr};
  TH2* m_hRawSeedVsCent{nullptr};
  TH2* m_hSubSeedVsCent{nullptr};

};

#endif  // JETSEEDCOUNT_H
