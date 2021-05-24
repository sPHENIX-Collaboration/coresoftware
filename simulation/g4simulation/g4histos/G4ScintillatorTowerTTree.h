#ifndef G4HISTOS_G4SCINTILLATORTOWERTTREE_H
#define G4HISTOS_G4SCINTILLATORTOWERTTREE_H

#include <fun4all/SubsysReco.h>

#include <string>

class Fun4AllHistoManager;
class PHCompositeNode;
class TH1;

class G4ScintillatorTowerTTree : public SubsysReco
{
 public:
  G4ScintillatorTowerTTree(const std::string &name = "SCINTILLATORTOWERTTREE");
  ~G4ScintillatorTowerTTree() override {}

  //! full initialization
  int Init(PHCompositeNode *) override;

  //! event processing method
  int process_event(PHCompositeNode *) override;

  int End(PHCompositeNode *) override;

  void Detector(const std::string &det);

  void SaveScintillatorTowers(const int i = 1) { savetowers = i; }

  void HistoFileName(const std::string &name) { _histofilename = name; }

 protected:
  std::string _detector;
  std::string _outnodename;
  std::string _towernodename;
  std::string _histofilename;
  int savetowers;
  int evtno;
  Fun4AllHistoManager *hm;
  TH1 *etot_hist;
};

#endif
