// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef UEVSETACENTRALITY_H
#define UEVSETACENTRALITY_H

// jet qa
#include "jetqa/JetQADefs.h"

#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllHistoManager.h>

// c++ utilities
#include <cstdlib>
#include <map>
#include <string>
#include <vector>

class PHCompositeNode;
class TH1F;
class TH2F;
class TriggerAnalyzer;

class UEvsEtaCentrality : public SubsysReco
{
 public:

  // set options
  struct Config
  {
    ///! turn debug messages on/off
    bool debug {true};

    ///! turn trigger selection on/off
    bool doTrgSelect {false};

    ///! trigger to select
    uint32_t trgToSelect {JetQADefs::GL1::MBDNSJet1};

    ///! histogram tag
    std::string histTag {""};

    ///! module name
    std::string moduleName {"UEvsEtaCentrality"};
  };

  UEvsEtaCentrality(const std::string &moduleName = "UEvsEtaCentrality");
  UEvsEtaCentrality(const Config& config);
  ~UEvsEtaCentrality() override;

  // setters
  void SetConfig(const Config& config) {m_config = config;}

  // getters
  const Config& GetConfig() {return m_config;}

  // f4a methods

  /** Called during initialization.
      Typically this is where you can book histograms, and e.g.
      register them to Fun4AllServer (so they can be output to file
      using Fun4AllServer::dumpHistos() method).
   */
  int Init(PHCompositeNode* /*topNode*/) override;

  /** Called for first event when run number is known.
      Typically this is where you may want to fetch data from
      database, because you know the run number. A place
      to book histograms which have to know the run number.
   */
  int InitRun(PHCompositeNode* /*topNode*/) override;

  /** Called for each event.
      This is where you do the real work.
   */
  int process_event(PHCompositeNode* topNode) override;
  
  /// Clean up internals after each event.
  int ResetEvent(PHCompositeNode* /*topNode*/) override;

  /// Called at the end of each run.
  int EndRun(const int runnumber) override;

  /// Called at the end of all processing.
  int End(PHCompositeNode* /*topNode*/) override;

  /// Reset
  int Reset(PHCompositeNode * /*topNode*/) override;

  void Print(const std::string &what = "ALL") const override;

 private:

  //std::string m_outputFileName;
  
  void InitHistManager();

  ///! module configuration
  Config m_config;

  std::string m_moduleName;

  ///! histogram manager
  Fun4AllHistoManager* m_manager {nullptr};

  TriggerAnalyzer* m_analyzer {nullptr};

  TH2F *hv2_cent = nullptr;
  TH2F *hPsi2_cent = nullptr;
  
  TH2F *hUEiHcalEta_Cent0_20 = nullptr;
  TH2F *hUEoHcalEta_Cent0_20 = nullptr;
  TH2F *hUEemcalEta_Cent0_20 = nullptr;
  
  TH2F *hUEiHcalEta_Cent20_50 = nullptr;
  TH2F *hUEoHcalEta_Cent20_50 = nullptr;
  TH2F *hUEemcalEta_Cent20_50 = nullptr;
  
  TH2F *hUEiHcalEta_Cent50_100 = nullptr;
  TH2F *hUEoHcalEta_Cent50_100 = nullptr;
  TH2F *hUEemcalEta_Cent50_100 = nullptr;

  TH2F *hUEiHcalEta = nullptr;
  TH2F *hUEoHcalEta = nullptr;
  TH2F *hUEemcalEta = nullptr;
};

#endif // UEVSETACENTRALITY_H
