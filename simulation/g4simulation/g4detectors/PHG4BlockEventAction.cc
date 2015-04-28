#include "PHG4BlockEventAction.h"
#include <g4main/PHG4HitContainer.h>

#include <Geant4/G4Event.hh>

#include <fun4all/getClass.h>

using namespace std;

//___________________________________________________
PHG4BlockEventAction::PHG4BlockEventAction( PHCompositeNode *node, const string &name ):
  nodename(name),
  topNode(node)
  {}


//___________________________________________________
void PHG4BlockEventAction::EndOfEventAction(const G4Event* evt)
{
  PHG4HitContainer* block_hits =  findNode::getClass<PHG4HitContainer>( topNode ,nodename.c_str());
  if (block_hits)
    {
      block_hits->RemoveZeroEDep();
    }
}
