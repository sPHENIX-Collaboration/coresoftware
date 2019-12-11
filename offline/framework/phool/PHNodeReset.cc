//  Implementation of class PHNodeReset
//  Author: Matthias Messer

#include "PHNodeReset.h"

#include "PHDataNode.h"
#include "PHNode.h"
#include "PHObject.h"

#include <iostream>
#include <string>

using namespace std;

void PHNodeReset::perform(PHNode* node)
{
  if (node->getResetFlag() != true) return;
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
