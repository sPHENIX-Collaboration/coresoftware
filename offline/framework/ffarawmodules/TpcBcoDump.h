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
class TNtuple;

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
  TFile *outTfile{nullptr};
  TNtuple *ntup{nullptr};
  std::map<int, uint64_t> lastbco;
  std::string outfilename;
};

#endif  // FFARAWMODULES_TPCBCODUMP_H
