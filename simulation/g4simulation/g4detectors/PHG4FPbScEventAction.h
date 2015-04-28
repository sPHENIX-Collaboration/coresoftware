#ifndef PHG4FPbScEventAction_h
#define PHG4FPbScEventAction_h

#include <g4main/PHG4EventAction.h>
#include <string>

class G4Event;
class PHCompositeNode;

class PHG4FPbScEventAction: public PHG4EventAction
{

public:

  //! constructor
  PHG4FPbScEventAction( PHCompositeNode *topNode, const std::string &name);

  //! destuctor
  virtual ~PHG4FPbScEventAction()
  {}

  void EndOfEventAction(const G4Event*);
  
 private:
  
  std::string nodename;
  PHCompositeNode *topNode;
  
};


#endif
