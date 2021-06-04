#ifndef G4HISTOS_G4SNGLNTUPLE_H
#define G4HISTOS_G4SNGLNTUPLE_H

#include <fun4all/SubsysReco.h>

#include <map>
#include <set>
#include <string>
#include <vector>

// Forward declerations
class Fun4AllHistoManager;
class PHCompositeNode;
class TFile;
class TH1;
class TNtuple;

class G4SnglNtuple : public SubsysReco
{
 public:
  //! constructor
  G4SnglNtuple(const std::string &name = "G4SnglNtuple", const std::string &filename = "G4SnglNtuple.root");

  //! destructor
  ~G4SnglNtuple() override;

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
  std::vector<TH1 *> nhits;
  std::vector<TH1 *> eloss;
  //  std::vector<TH2 *> nhit_edep;
  std::string _filename;
  std::set<std::string> _node_postfix;
  std::map<std::string, int> _detid;
  TNtuple *ntup;
  TNtuple *ntup_e;
  TFile *outfile;
};

#endif
