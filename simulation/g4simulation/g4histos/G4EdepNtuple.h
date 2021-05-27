#ifndef G4HISTOS_G4EDEPNTUPLE_H
#define G4HISTOS_G4EDEPNTUPLE_H

#include <fun4all/SubsysReco.h>

#include <map>
#include <set>
#include <string>

// Forward declerations
class Fun4AllHistoManager;
class PHCompositeNode;
class TFile;
class TNtuple;

class G4EdepNtuple : public SubsysReco
{
 public:
  //! constructor
  G4EdepNtuple(const std::string &name = "G4EdepNtuple", const std::string &filename = "G4EdepNtuple.root");

  //! destructor
  ~G4EdepNtuple() override;

  //! full initialization
  int Init(PHCompositeNode *) override;

  //! event processing method
  int process_event(PHCompositeNode *) override;

  //! end of run method
  int End(PHCompositeNode *) override;

  void AddNode(const std::string &name, const int detid = 0);

 protected:
  int nblocks;
  Fun4AllHistoManager *hm;
  std::string _filename;
  std::set<std::string> _node_postfix;
  std::map<std::string, int> _detid;
  TNtuple *ntup;
  TFile *outfile;
};

#endif
