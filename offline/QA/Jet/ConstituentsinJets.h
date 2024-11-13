// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef QA_JET_CONSTITUENTSINJETS_H
#define QA_JET_CONSTITUENTSINJETS_H

#include <fun4all/SubsysReco.h>
#include <jetbase/Jet.h>

#include <string>
#include <utility>  // std::pair, std::make_pair
#include <vector>

#include "JetQADefs.h"

class Fun4AllHistoManager;
class PHCompositeNode;
class TH1;
class TH2;

class ConstituentsinJets : public SubsysReco
{
 public:
  ConstituentsinJets(
      const std::string &moduleName = "ConstituentsInJets",
      const std::string &recojetname = "AntiKt_Tower_r04",
      const std::string &towBkgdName = "TowerInfoBackground_Sub2",
      const std::string &histTag = "AllTrig_AntiKt_Tower_R04");
  ~ConstituentsinJets() override{};

  void setRecoJetNodeName(const std::string &name)
  {  // set the name of the node containing the reco jets
    m_recoJetName = name;
  }
  void setHistTag(const std::string &tag)
  {  // set the tag to be applied to histogram names
    m_histTag = tag;
  }

  void setEtaRange(double low, double high)
  {  // set the eta range for the reco jets
    m_etaRange.first = low;
    m_etaRange.second = high;
  }
  void setPtRange(double low, double high)
  {  // set the pt range for the reco jets
    m_ptRange.first = low;
    m_ptRange.second = high;
  }
  void setTrgToSelect(const uint32_t trig = JetQADefs::GL1::MBDNSJet1)
  {  // specifies a trigger to require
    m_doTrgSelect = true;
    m_trgToSelect = trig;
  }

  // standard Fun4All functions
  int Init(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

 private:
  //! Module name, input node strings, and histogram tags
  std::string m_moduleName{"ConstituentsinJets"};
  std::string m_recoJetName{"AntiKt_Tower_r04"};
  std::string m_towBkgdName{"TowerInfoBackground_Sub2"};
  std::string m_histTag{"AllTrig_AntiKt_Tower_R04"};
  // std::string m_outputFileName{ "ConstituentsinJets.root"};

  //! Trigger selection
  bool m_doTrgSelect{false};
  uint32_t m_trgToSelect{JetQADefs::GL1::MBDNSJet1};

  // ! Kinematic cuts and reco jet node name
  std::pair<double, double> m_etaRange = std::make_pair(-1.1, 1.1);
  std::pair<double, double> m_ptRange = std::make_pair(1.0, 1000.0);

  // Jet N constituents
  Fun4AllHistoManager *m_manager{nullptr};

  TH1 *h1_ConstituentsinJets_total{nullptr};
  // TH1 * h1_ConstituentsinJets_CaloTowers{nullptr};
  TH1 *h1_ConstituentsinJets_IHCAL{nullptr};
  TH1 *h1_ConstituentsinJets_OHCAL{nullptr};
  TH1 *h1_ConstituentsinJets_CEMC{nullptr};
  TH2 *h2_ConstituentsinJets_vs_caloLayer{nullptr};

  // Jet E fraction
  TH1 *h1_jetFracE_IHCAL{nullptr};
  TH1 *h1_jetFracE_OHCAL{nullptr};
  TH1 *h1_jetFracE_CEMC{nullptr};
  TH2 *h2_jetFracE_vs_caloLayer{nullptr};
};

#endif  // ConstituentsinJets_H
