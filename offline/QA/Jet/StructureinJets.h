
// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef QA_JET_STRUCTUREINJETS_H
#define QA_JET_STRUCTUREINJETS_H

#include "JetQADefs.h"

#include <fun4all/SubsysReco.h>

#include <string>
#include <utility>

class Fun4AllHistoManager;
class PHCompositeNode;
class TH1;
class TH2;
class TH3;
class TriggerAnalyzer;

class StructureinJets : public SubsysReco
{
 public:
  StructureinJets(const std::string& moduleName = "StructureInJets",
                  const std::string& recojetname = "AntiKt_Tower_r04_Sub1",
                  const std::string& trkNodeName = "SvtxTrackMap",
                  const std::string& histTag = "AllTrig_AntiKt_Tower_r04",
                  const std::string& outputfilename = "tracksinjets.root");

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

  bool isAA() const { return m_isAAFlag; }
  void isAA(bool b) { m_isAAFlag = b; }

  bool writeToOutputFile() const { return m_writeToOutputFileFlag; }
  void writeToOutputFile(bool b) { m_writeToOutputFileFlag = b; }

  /// set the name of the node containing the reco jets
  void setRecoJetNodeName(const std::string &name)
  {
    m_recoJetName = name;
  }

  /// set the name of the node containg the tracks
  void setTrkNodeName(const std::string &name)
  {
    m_trkNodeName = name;
  }

  /// set the name of the output file
  void setOutputFileName(const std::string &name)
  {
    m_outputFileName = name;
  }

  /// set the tag to be applied to the histogram names
  void setHistTag(const std::string &tag)
  {
    m_histTag = tag;
  }

  /// set trigger to require
  void setTrgToSelect(const uint32_t trig = JetQADefs::GL1::MBDNSJet1)
  {
    m_doTrgSelect = true;
    m_trgToSelect = trig;
  }

  /// set minimum track pt
  void setTrkPtCut(const float cut)
  {
    m_trk_pt_cut = cut;
  }

  /// set max track quality
  void setTrkQualityCut(const float cut)
  {
    m_trk_qual_cut = cut;
  }

  /// set min no. of silicon clusters
  void setTrkNSilCut(const int cut)
  {
    m_trk_nsil_cut = cut;
  }

  /// set min no. of tpc clusters
  void setTrkNTPCCut(const int cut)
  {
    m_trk_ntpc_cut = cut;
  }

  /// set jet radius
  void setJetRadius(const float radius)
  {
    m_jetRadius = radius;
  }

  /// set jet pt range
  void setJetPtRange(const double low, const double high)
  {
    m_ptJetRange.first = low;
    m_ptJetRange.second = high;
  }

  /// set jet eta range
  void setJetEtaRange(const double low, const double high)
  {
    m_etaJetRange.first = low;
    m_etaJetRange.second = high;
  }

 private:

  /// module name, input node strings, histogram tags, and output file name
  std::string m_moduleName;
  std::string m_recoJetName;
  std::string m_trkNodeName;
  std::string m_histTag;
  std::string m_outputFileName;

  /// track cuts
  float m_trk_pt_cut{0.1};
  float m_trk_qual_cut{6.0};
  int m_trk_nsil_cut{4};
  int m_trk_ntpc_cut{25};
  float m_jetRadius{0.4};

  /// jet kinematic cuts
  std::pair<double, double> m_etaJetRange{-1.1, 1.1};
  std::pair<double, double> m_ptJetRange{1.0, 1000.0};

  /// flags
  bool m_isAAFlag{false};
  bool m_writeToOutputFileFlag{false};
  bool m_doTrgSelect{false};

  /// trigger to select
  uint32_t m_trgToSelect{JetQADefs::GL1::MBDNSJet1};

  /// output histograms
  TH3 *m_hSumTrkVsJetPtVsCent{nullptr};
  TH2 *m_hSumTrkVsJetPt{nullptr};
  TH1 *m_hSumTrkOverJetPt{nullptr};
  TH1 *m_hSumTrkPt{nullptr};

  /// fun4all members
  Fun4AllHistoManager *m_manager{nullptr};
  TriggerAnalyzer *m_analyzer{nullptr};
};

#endif  // QA_JET_STRUCTUREINJETS_H
