#ifndef G4HISTOS_G4CELLNTUPLE_H
#define G4HISTOS_G4CELLNTUPLE_H

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
class TH2;
class TNtuple;

class G4CellNtuple : public SubsysReco
{
 public:
  //! constructor
  G4CellNtuple(const std::string &name = "G4CellNtuple", const std::string &filename = "G4CellNtuple.root");

  //! destructor
  virtual ~G4CellNtuple();

  //! full initialization
  int Init(PHCompositeNode *);

  //! event processing method
  int process_event(PHCompositeNode *);

  //! end of run method
  int End(PHCompositeNode *);

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
  TFile *outfile;
};

#endif
