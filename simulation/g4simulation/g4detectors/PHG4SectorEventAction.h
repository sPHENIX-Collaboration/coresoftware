#ifndef PHG4SectorEventAction_h
#define PHG4SectorEventAction_h

#include <g4main/PHG4EventAction.h>
#include <string>

class G4Event;
class PHCompositeNode;

class PHG4SectorEventAction: public PHG4EventAction
{

public:

  //! constructor
  PHG4SectorEventAction( PHCompositeNode *topNode, const std::string &name);

  //! destuctor
  virtual ~PHG4SectorEventAction()
  {}

  void EndOfEventAction(const G4Event*);
  
 private:
  
  std::string nodename;
  PHCompositeNode *topNode;
  
};


#endif
