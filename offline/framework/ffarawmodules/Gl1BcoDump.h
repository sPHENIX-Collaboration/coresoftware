// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FFARAWMODULES_GL1BCODUMP_H
#define FFARAWMODULES_GL1BCODUMP_H

#include <fun4all/SubsysReco.h>

#include <fstream>
#include <map>
#include <set>
#include <string>

class Fun4AllInputManager;
class PHCompositeNode;
class TFile;
class TTree;

class Gl1BcoDump : public SubsysReco
{
 public:
  Gl1BcoDump(const std::string &name = "Gl1BcoDump");

  ~Gl1BcoDump() override {}

  int InitRun(PHCompositeNode *topNode) override;

  int process_event(PHCompositeNode *topNode) override;

  int End(PHCompositeNode *topNode) override;
  //  int ResetEvent(PHCompositeNode *topNode) override;
  void OutFileName(const std::string &name) { outfilename = name; }

 private:
  TFile *outfile{nullptr};
  TTree *ttree{nullptr};
  uint64_t lastbco{0};
  std::string outfilename;
  int m_id{0};
  int m_evt{0};
  uint64_t m_bco{0};
  int64_t m_bcodiff{0};
};

#endif  // FFARAWMODULES_GL1BCODUMP_H
