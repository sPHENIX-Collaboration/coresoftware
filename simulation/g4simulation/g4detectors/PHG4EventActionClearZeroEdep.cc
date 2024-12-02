#include "PHG4EventActionClearZeroEdep.h"
#include <g4main/PHG4HitContainer.h>

#include <phool/getClass.h>

//___________________________________________________
PHG4EventActionClearZeroEdep::PHG4EventActionClearZeroEdep(PHCompositeNode *node, const std::string &name)
  : topNode(node)
{
  AddNode(name);
}

void PHG4EventActionClearZeroEdep::AddNode(const std::string &name)
{
  nodename_set.insert(name);
}

//___________________________________________________
void PHG4EventActionClearZeroEdep::EndOfEventAction(const G4Event * /*evt*/)
{
  for (auto &node : nodename_set)
  {
    PHG4HitContainer *generic_hits = findNode::getClass<PHG4HitContainer>(topNode, node);
    if (generic_hits)
    {
      generic_hits->RemoveZeroEDep();
    }
  }
}
