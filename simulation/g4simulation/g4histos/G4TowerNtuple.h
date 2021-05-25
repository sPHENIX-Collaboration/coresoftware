#ifndef G4HISTOS_G4TOWERNTUPLE_H
#define G4HISTOS_G4TOWERNTUPLE_H

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

class G4TowerNtuple : public SubsysReco
{
 public:
  //! constructor
  G4TowerNtuple(const std::string &name = "G4TowerNtuple", const std::string &filename = "G4TowerNtuple.root");

  //! destructor
  ~G4TowerNtuple() override;

  //! full initialization
  int Init(PHCompositeNode *) override;

  //! event processing method
  int process_event(PHCompositeNode *) override;

  //! end of run method
  int End(PHCompositeNode *) override;

  void AddNode(const std::string &name, const std::string &twrtype, const int detid);

 protected:
  int nblocks;
  Fun4AllHistoManager *hm;
  std::vector<TH1 *> nhits;
  std::vector<TH1 *> eloss;
  //  std::vector<TH2 *> nhit_edep;
  std::string _filename;
  std::set<std::string> _node_postfix;
  std::map<std::string, std::string> _tower_type;
  std::map<std::string, int> _detid;
  TNtuple *ntup;
  TFile *outfile;
};

#endif
