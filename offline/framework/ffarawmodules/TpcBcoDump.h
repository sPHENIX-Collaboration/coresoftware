// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FFARAWMODULES_TPCBCODUMP_H
#define FFARAWMODULES_TPCBCODUMP_H

#include <fun4all/SubsysReco.h>

#include <fstream>
#include <map>
#include <set>
#include <string>

class Fun4AllInputManager;
class PHCompositeNode;
class TFile;
class TTree;

class TpcBcoDump : public SubsysReco
{
 public:
  TpcBcoDump(const std::string &name = "TpcBcoDump");

  ~TpcBcoDump() override {}

  int InitRun(PHCompositeNode *topNode) override;

  int process_event(PHCompositeNode *topNode) override;

  int End(PHCompositeNode *topNode) override;
  //  int ResetEvent(PHCompositeNode *topNode) override;
  void OutFileName(const std::string &name) { outfilename = name; }

 private:
  TFile *outfile{nullptr};
  TTree *ttree = nullptr;
  std::map<int, uint64_t> lastbco;
  std::string outfilename;

  int m_id{0};
  int m_evt{0};
  uint64_t m_bco{0};
  int64_t m_bcodiff{0};
};

#endif  // FFARAWMODULES_TPCBCODUMP_H
