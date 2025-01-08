
// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef QA_JET_STRUCTUREINJETS_H
#define QA_JET_STRUCTUREINJETS_H

#include "JetQADefs.h"

#include <fun4all/SubsysReco.h>

#include <string>

class Fun4AllHistoManager;
class PHCompositeNode;
class TH2;
class TH3;
class TriggerAnalyzer;

class StructureinJets : public SubsysReco
{
 public:
  StructureinJets(const std::string &moduleName = "StructureInJets",
                  const std::string &recojetname = "AntiKt_Tower_r04",
                  const std::string &histTag = "AllTrig_AntiKt_Tower_r04",
                  const std::string &outputfilename = "tracksinjets.root");

  ~StructureinJets() override;

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

  bool isAA() const { return isAAFlag; }
  void isAA(bool b) { isAAFlag = b; }

  bool writeToOutputFile() const { return writeToOutputFileFlag; }
  void writeToOutputFile(bool b) { writeToOutputFileFlag = b; }

  // set trigger to require
  void setTrgToSelect(const uint32_t trig = JetQADefs::GL1::MBDNSJet1)
  {
    m_doTrgSelect = true;
    m_trgToSelect = trig;
  }

 private:
  std::string m_moduleName;
  std::string m_recoJetName;
  std::string m_histTag;
  float m_trk_pt_cut{2};
  float m_jetRadius{0.4};
  bool isAAFlag{false};
  bool writeToOutputFileFlag{false};
  bool m_doTrgSelect{false};
  uint32_t m_trgToSelect{JetQADefs::GL1::MBDNSJet1};
  std::string m_outputFileName;
  TH3 *m_h_track_vs_calo_pt{nullptr};
  TH2 *m_h_track_pt{nullptr};
  Fun4AllHistoManager *m_manager{nullptr};
  TriggerAnalyzer *m_analyzer{nullptr};
};

#endif  // QA_JET_STRUCTUREINJETS_H
