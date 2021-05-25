#ifndef G4HISTOS_G4SNGLTREE_H
#define G4HISTOS_G4SNGLTREE_H

#include "G4EvtTree.h"

#include <fun4all/SubsysReco.h>

#include <map>
#include <set>
#include <string>

// Forward declarations
class PHCompositeNode;
class PHG4HitContainer;
class TFile;
class TTree;

class G4SnglTree : public SubsysReco
{
 public:
  //! constructor
  G4SnglTree(const std::string &name = "G4SnglTree", const std::string &filename = "G4SnglTree.root");

  //! destructor
  ~G4SnglTree() override {}

  //! full initialization
  int Init(PHCompositeNode *) override;

  //! event processing method
  int process_event(PHCompositeNode *) override;

  //! hit processing method
  int process_hit(PHG4HitContainer *hits, const std::string &dName, int detid, int &nhits);

  //! end of run method
  int End(PHCompositeNode *) override;

  void AddNode(const std::string &name, const int detid = 0);

 protected:
  int nblocks;
  //  std::vector<TH2 *> nhit_edep;
  std::string _filename;
  std::set<std::string> _node_postfix;
  std::map<std::string, int> _detid;

  TTree *g4tree;
  G4EvtTree mG4EvtTree;
  TFile *outfile;
};

#endif
