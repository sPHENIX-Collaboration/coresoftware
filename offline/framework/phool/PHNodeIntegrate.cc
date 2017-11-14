//  Implementation of class PHNodeIntegrate
//  Author: Matthias Messer

#include "PHNodeIntegrate.h"

#include "getClass.h"
#include "PHCompositeNode.h"
#include "PHDataNode.h"
#include "PHIODataNode.h"
#include "PHObject.h"

#include <iostream>

using namespace std;

void PHNodeIntegrate::perform(PHNode* node)
{
  if (verbosity > 0)
  {
    cout << "PHNodeIntegrate: Integrateting " << node->getName() << endl;
  }
  if (node->getType() == "PHDataNode" || node->getType() == "PHIODataNode")
  {
    if (node->getObjectType() == "PHObject")
    {
      PHObject *obj = static_cast<PHDataNode<PHObject>*>(node)->getData();
      if (obj->Integrate(nullptr) > 0)
      {
	cout << "finding " << node->getName() << " on " <<  runsumnode->getName() << endl;
	PHObject *sumobj = findNode::getClass<PHObject>( runsumnode,node->getName());
	if (sumobj)
	{
	  cout << "found sumobj on " <<  runsumnode->getName() << endl;
	  sumobj->identify();
	  sumobj->Integrate(obj);
	  cout << "after integration" << endl;
	  sumobj->identify();
// now finally copy the summed object to the runwise node
// since the runwise nodes were copied from the input file we are guaranteed
// that it exists (and we handle only objects which come from the current
// input file
	PHObject *runobj = findNode::getClass<PHObject>( runnode,node->getName());
	runobj->CopyContent(sumobj);
	}
	else
	{
          sumobj =  obj->clone();
 PHIODataNode<PHObject> *intnode = new  PHIODataNode<PHObject>(sumobj,node->getName(),"PHObject");
 runsumnode->addNode(intnode);
	  cout << "needed to clone " << runsumnode->getName() << endl;
	  sumobj->identify();
	}
	//     int iret = (static_cast<PHDataNode<PHObject>*>(node))->getData()->Integrate(obj);

PHObject *runobj = findNode::getClass<PHObject>( runnode,node->getName());
cout << "final run object:" << endl;
runobj->identify();
      }
    }
  }
}
