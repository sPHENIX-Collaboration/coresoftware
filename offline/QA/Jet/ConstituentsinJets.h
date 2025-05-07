// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef QA_JET_CONSTITUENTSINJETS_H
#define QA_JET_CONSTITUENTSINJETS_H

#include "JetQADefs.h"

#include <fun4all/SubsysReco.h>
#include <jetbase/Jet.h>

#include <string>
#include <utility>  // std::pair, std::make_pair
#include <vector>

class Fun4AllHistoManager;
class PHCompositeNode;
class TH1;
class TH2;
class TriggerAnalyzer;

class ConstituentsinJets : public SubsysReco
{
 public:
  ConstituentsinJets(
      const std::string &moduleName = "ConstituentsInJets",
      const std::string &recojetname = "AntiKt_Tower_r04",
      const std::string &towBkgdName = "TowerInfoBackground_Sub2",
      const std::string &histTag = "AllTrig_AntiKt_Tower_R04",
      const std::string &towCEMCName = "TOWERINFO_CALIB_CEMC_RETOWER",
      const std::string &towIHCALName = "TOWERINFO_CALIB_HCALIN",      
      const std::string &towOHCALName = "TOWERINFO_CALIB_HCALOUT");
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
  void setPPMode(const bool pp)
  {
    m_inPPMode = pp;
  }

  void setTowBkgdNodeName(const std::string &name)
  { // set the name of the node containing the subtracted background towers
    m_towBkgdName = name;
  }

  void setTowNodeNameCEMC(const std::string &name)
  {//set the name of the node containing raw towers from EMCAL
    m_towCEMCName = name;
  }

  void setTowNodeNameIHCAL(const std::string &name)
  {
    m_towIHCALName = name;  
  }

  void setTowNodeNameOHCAL(const std::string &name)
  {
    m_towOHCALName = name;
  }


  // standard Fun4All functions
  int Init(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

 private:
  //! Module name, input node strings, and histogram tags
  std::string m_moduleName;
  std::string m_recoJetName;
  std::string m_towBkgdName;
  std::string m_histTag;
  std::string m_towCEMCName;
  std::string m_towIHCALName;
  std::string m_towOHCALName;
  // std::string m_outputFileName{ "ConstituentsinJets.root"};

  //! Trigger selection
  bool m_doTrgSelect{false};
  bool m_inPPMode{false};
  uint32_t m_trgToSelect{JetQADefs::GL1::MBDNSJet1};

  // ! Kinematic cuts and reco jet node name
  std::pair<double, double> m_etaRange{-1.1, 1.1};
  std::pair<double, double> m_ptRange{1.0, 1000.0};

  // Jet N constituents
  Fun4AllHistoManager *m_manager{nullptr};

  TriggerAnalyzer* m_analyzer{nullptr};

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
