#ifndef G4HISTOS_G4SCINTILLATORSLATTTREE_H
#define G4HISTOS_G4SCINTILLATORSLATTTREE_H

#include <fun4all/SubsysReco.h>

#include <string>

class Fun4AllHistoManager;
class PHCompositeNode;
class TH1;

class G4ScintillatorSlatTTree : public SubsysReco
{
 public:
  G4ScintillatorSlatTTree(const std::string &name = "SCINTILLATORSLATTTREE");
  ~G4ScintillatorSlatTTree() override = default;

  //! full initialization
  int Init(PHCompositeNode *) override;

  //! event processing method
  int process_event(PHCompositeNode *) override;

  int End(PHCompositeNode *) override;

  void Detector(const std::string &det);

  void SaveScintillatorSlats(const int i = 1) { saveslats = i; }

  void HistoFileName(const std::string &name) { _histofilename = name; }

 private:
  std::string _detector;
  std::string _outnodename;
  std::string _slatnodename;
  std::string _histofilename;
  int saveslats{1};
  int evtno{0};
  Fun4AllHistoManager *hm{nullptr};
  TH1 *etot_hist{nullptr};
};

#endif
