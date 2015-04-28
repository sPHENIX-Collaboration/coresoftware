#ifndef PHG4ConeEventAction_h
#define PHG4ConeEventAction_h

#include <g4main/PHG4EventAction.h>
#include <string>

class G4Event;
class PHCompositeNode;

class PHG4ConeEventAction: public PHG4EventAction
{

public:

  //! constructor
  PHG4ConeEventAction( PHCompositeNode *topNode, const std::string &name);

  //! destuctor
  virtual ~PHG4ConeEventAction()
  {}

  void EndOfEventAction(const G4Event*);
  
 private:
  
  std::string nodename;
  PHCompositeNode *topNode;
  
};


#endif
