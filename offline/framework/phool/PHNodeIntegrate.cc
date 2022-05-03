#include "PHNodeIntegrate.h"

#include "PHCompositeNode.h"
#include "PHDataNode.h"
#include "PHIODataNode.h"
#include "PHNode.h"
#include "PHObject.h"
#include "getClass.h"

#include <iostream>
#include <string>

void PHNodeIntegrate::perform(PHNode *node)
{
  if (verbosity > 0)
  {
    std::cout << "PHNodeIntegrate: Integrating " << node->getName() << std::endl;
  }
  if (node->getType() == "PHDataNode" || node->getType() == "PHIODataNode")
  {
    if (node->getObjectType() == "PHObject")
    {
      PHObject *obj = static_cast<PHDataNode<PHObject> *>(node)->getData();
      if (obj->Integrate())
      {
        PHObject *sumobj = findNode::getClass<PHObject>(runsumnode, node->getName());
        if (sumobj)
        {
          sumobj->Integrate(obj);
          // now copy the summed object to the runwise node
          // since the runwise nodes were copied from the input file
          // we are guaranteed that it exists (and we handle only
          // objects which come from the current input file
          PHObject *runobj = findNode::getClass<PHObject>(runnode, node->getName());
          runobj->CopyFrom(sumobj);
        }
        else
        {
          // since this object was also copied to the node tree we only need
          // to store it in case a second file gets opened where this one then
          // serves as the object which contains the sum
          sumobj = obj->CloneMe();
          PHIODataNode<PHObject> *sumobjnode = new PHIODataNode<PHObject>(sumobj, node->getName(), "PHObject");
          runsumnode->addNode(sumobjnode);
        }
      }
    }
  }
}
