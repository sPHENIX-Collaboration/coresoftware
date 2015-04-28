#ifndef PHG4CylinderEventAction_h
#define PHG4CylinderEventAction_h

#include <g4main/PHG4EventAction.h>

#include <string>
#include <set>

class G4Event;
class PHCompositeNode;

class PHG4CylinderEventAction: public PHG4EventAction
{

public:

  //! constructor
  PHG4CylinderEventAction( PHCompositeNode *topNode, const std::string &name);

  void AddNode(const std::string &name);

  //! destuctor
  virtual ~PHG4CylinderEventAction()
  {}

  void EndOfEventAction(const G4Event*);
  
 private:
  
  std::set<std::string> nodename_set;
  PHCompositeNode *topNode;
  
};


#endif
