//  Implementation of class PHNodeReset
//  Author: Matthias Messer

#include "PHNodeReset.h"

#include "PHDataNode.h"
#include "PHIODataNode.h"
#include "PHObject.h"

#include <iostream>

using namespace std;

void PHNodeReset::perform(PHNode* node)
{
  if (node->getResetFlag() != True) return;
  if (verbosity > 0)
  {
    cout << "PHNodeReset: Resetting " << node->getName() << endl;
  }
  if (node->getType() == "PHDataNode" || node->getType() == "PHIODataNode")
  {
    if (node->getObjectType() == "PHObject")
    {
      (static_cast<PHDataNode<PHObject>*>(node))->getData()->Reset();
    }
  }
}
