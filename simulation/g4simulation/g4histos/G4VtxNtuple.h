#ifndef G4HISTOS_G4VTXNTUPLE_H
#define G4HISTOS_G4VTXNTUPLE_H

#include <fun4all/SubsysReco.h>

#include <string>

// Forward declerations
class Fun4AllHistoManager;
class PHCompositeNode;
class TNtuple;

class G4VtxNtuple : public SubsysReco
{
 public:
  //! constructor
  G4VtxNtuple(const std::string &name = "G4VtxNtuple", const std::string &filename = "G4VtxNtuple.root");

  //! destructor
  ~G4VtxNtuple() override;

  //! full initialization
  int Init(PHCompositeNode *) override;

  //! event processing method
  int process_event(PHCompositeNode *) override;

  //! end of run method
  int End(PHCompositeNode *) override;

 protected:
  std::string m_FileName;
  Fun4AllHistoManager *hm = nullptr;
  TNtuple *ntup = nullptr;
};

#endif
