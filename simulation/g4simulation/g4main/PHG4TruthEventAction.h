// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4TRUTHEVENTACTION_H
#define G4MAIN_PHG4TRUTHEVENTACTION_H

#include "PHG4EventAction.h"

#include <map>
#include <set>

class G4Event;
class PHG4HitContainer;
class PHG4TruthInfoContainer;
class PHCompositeNode;

class PHG4TruthEventAction : public PHG4EventAction
{
 public:
  //! constructor
  PHG4TruthEventAction();

  //! destuctor
  ~PHG4TruthEventAction() override {}

  void BeginOfEventAction(const G4Event*) override;

  void EndOfEventAction(const G4Event*) override;

  int ResetEvent(PHCompositeNode*) override;

  //! get relevant nodes from top node passed as argument
  void SetInterfacePointers(PHCompositeNode*) override;

  //! add id into track list
  void AddTrackidToWritelist(const int trackid);

 private:
  void SearchNode(PHCompositeNode* topNode);
  void PruneShowers();
  void ProcessShowers();

  //! set of track ids to be written out
  std::set<int> m_WriteSet;

  //! pointer to truth information container
  PHG4TruthInfoContainer* m_TruthInfoContainer;

  int m_LowerKeyPrevExist;
  int m_UpperKeyPrevExist;

  std::map<int, PHG4HitContainer*> m_HitContainerMap;
};

#endif
