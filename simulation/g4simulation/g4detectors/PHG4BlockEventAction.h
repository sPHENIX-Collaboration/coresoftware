#ifndef PHG4BlockEventAction_h
#define PHG4BlockEventAction_h

#include <g4main/PHG4EventAction.h>
#include <string>

class G4Event;
class PHCompositeNode;

class PHG4BlockEventAction: public PHG4EventAction
{

public:

  //! constructor
  PHG4BlockEventAction( PHCompositeNode *topNode, const std::string &name);

  //! destuctor
  virtual ~PHG4BlockEventAction()
  {}

  void EndOfEventAction(const G4Event*);
  
 private:
  
  std::string nodename;
  PHCompositeNode *topNode;
  
};


#endif
