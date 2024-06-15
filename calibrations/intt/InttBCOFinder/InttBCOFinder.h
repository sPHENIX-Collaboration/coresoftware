#ifndef INTTBCOFINDER_INTTBCOFINDER_H
#define INTTBCOFINDER_INTTBCOFINDER_H

#include <fun4all/SubsysReco.h>

#include <string>

class CDBTTree;
class PHCompositeNode;
class TFile;
class TH2;

class InttBCOFinder : public SubsysReco
{
 public:
  InttBCOFinder(const std::string &name = "InttBCOFinder", const std::string &fname = "outputfile.root", const std::string &fname2 = "cdbfile.root", int nevent = 10000);

  virtual ~InttBCOFinder();

  int Init(PHCompositeNode *);

  int InitRun(PHCompositeNode *);

  /// SubsysReco event processing method
  int process_event(PHCompositeNode *);

  /// SubsysReco end processing method
  int End(PHCompositeNode *);

  void FindBCOPeak();
  void ADCCut(const bool flag) { IsADCcutON_ = flag; }
  void WriteCDBTTree(const bool flag) { WriteCDBTTree_ = flag; }
  void WriteQAFile(const bool flag) { WriteQAFile_ = flag; }

 private:
  TFile *outFile_{nullptr};
  CDBTTree *cdbttree_{nullptr};
  TH2 *h2_bco_ladder_[8]{};      // histogram for BCO alignment check half ladder by half ladder
  TH2 *h2_bco_ladder_cut_[8]{};  // histogram after BCO cuto

  int nevents_{0};
  int ievent_{0};

  bool IsADCcutON_{false};
  bool WriteCDBTTree_{false};
  bool WriteQAFile_{false};

  std::string m_InttRawNodeName{"INTTRAWHIT"};
  std::string outfname_;
  std::string cdbname_;
};
#endif
