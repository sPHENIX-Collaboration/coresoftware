// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef EMCLUSTERKINEMATICS_H
#define EMCLUSTERKINEMATICS_H

#include "JetQADefs.h"

#include <fun4all/SubsysReco.h>

#include <string>

class Fun4AllHistoManager;
class PHCompositeNode;
class TH1;
class TH2;
class TriggerAnalyzer;

class EMClusterKinematics : public SubsysReco
{
 public:


  EMClusterKinematics(const std::string &modulename = "EMClusterKinematics", const std::string &inputnode = "CLUSTERINFO_CEMC", const std::string &histtag = "");
  ~EMClusterKinematics() override;

  /** Called during initialization.
      Typically this is where you can book histograms, and e.g.
      register them to Fun4AllServer (so they can be output to file
      using Fun4AllServer::dumpHistos() method).
   */
  int Init(PHCompositeNode *topNode) override;

  /** Called for first event when run number is known.
      Typically this is where you may want to fetch data from
      database, because you know the run number. A place
      to book histograms which have to know the run number.
   */
  int InitRun(PHCompositeNode *topNode) override;

  /** Called for each event.
      This is where you do the real work.
   */
  int process_event(PHCompositeNode *topNode) override;

  /// Clean up internals after each event.
  int ResetEvent(PHCompositeNode *topNode) override;

  /// Called at the end of each run.
  int EndRun(const int runnumber) override;

  /// Called at the end of all processing.
  int End(PHCompositeNode *topNode) override;

  /// Reset
  int Reset(PHCompositeNode * /*topNode*/) override;

  void Print(const std::string &what = "ALL") const override;

  /// specifies a trigger to select
  void SetTrgToSelect(const uint32_t trig = JetQADefs::GL1::MBDNSPhoton1)
  {
    m_doTrgSelect = true;
    m_trgToSelect = trig;
  }

  /// set histogram tag
  void SetHistTag(const std::string& tag)
  {
    m_histtag = tag;
  }

  /// turn on/off optional histograms
  void SetDoOptHist(const bool opt)
  {
    m_doOptHist = opt;
  }

 private:
  // std::string outfilename;
  // TFile *outfile;
 // hist manager

 TriggerAnalyzer* m_analyzer {nullptr};
 Fun4AllHistoManager* m_manager {nullptr};
 std::string m_modulename;
 std::string m_inputnode;
 std::string m_histtag;
 uint32_t m_trgToSelect;
 bool m_doTrgSelect;
 bool m_doOptHist;

 ///Output histograms
 TH1 *h_emcal_cluster_chi2 {nullptr};
 TH1 *h_emcal_cluster_energy {nullptr};
 TH2 *h_emcal_cluster_eta_phi {nullptr};
 TH1 *h_emcal_cluster_eta {nullptr};    // <-- Declare eta histogram
 TH1 *h_emcal_cluster_phi {nullptr};    // <-- Declare phi histogram
 TH2 *h_emcal_cluster_energy_eta {nullptr};
 TH2 *h_emcal_cluster_chi2_eta {nullptr};
 TH1 *h_emcal_cluster_eta_with_cuts{nullptr};
 TH1 *h_emcal_cluster_phi_with_cuts{nullptr};
 TH2 *h_emcal_cluster_eta_phi_with_cuts{nullptr};
 TH2 *h_emcal_cluster_energy_phi {nullptr};
 TH2 *h_emcal_cluster_chi2_phi {nullptr};
 TH1 *h_emcal_cluster_eta_with_energy_cut{nullptr};
 TH2 *h_emcal_cluster_chi2_energy{nullptr};
};

#endif // EMCLUSTERKINEMATICS_H
