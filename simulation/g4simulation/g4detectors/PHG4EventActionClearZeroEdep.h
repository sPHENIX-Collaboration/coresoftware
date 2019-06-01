// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4EVENTACTIONCLEARZEROEDEP_H
#define G4DETECTORS_PHG4EVENTACTIONCLEARZEROEDEP_H

#include <g4main/PHG4EventAction.h>

#include <string>
#include <set>

class G4Event;
class PHCompositeNode;

class PHG4EventActionClearZeroEdep: public PHG4EventAction
{

public:

  //! constructor
  PHG4EventActionClearZeroEdep( PHCompositeNode *topNode, const std::string &name);

  void AddNode(const std::string &name);

  //! destuctor
  virtual ~PHG4EventActionClearZeroEdep()
  {}

  void EndOfEventAction(const G4Event*);
  
 private:
  
  std::set<std::string> nodename_set;
  PHCompositeNode *topNode;
  
};


#endif
