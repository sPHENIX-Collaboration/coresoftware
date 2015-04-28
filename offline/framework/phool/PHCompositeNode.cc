//-----------------------------------------------------------------------------
//
//  The PHOOL's Software
//  Copyright (C) PHENIX collaboration, 1999
//
//  Implementation of class PHCompositeNode
//
//-----------------------------------------------------------------------------
#include "PHCompositeNode.h" 
#include "PHPointerListIterator.h"

#include <iostream>

using namespace std;

PHCompositeNode::PHCompositeNode() : PHNode("NULL")
{
}

PHCompositeNode::PHCompositeNode(const PHString& name) : PHNode(name)
{
  type = "PHCompositeNode";
  objecttype = "PHCompositeNode";
  deleteMe = 0;
}

PHCompositeNode::~PHCompositeNode() 
{
  // we need to mark this node to be to deleted
  // The recursive taking out of deleted subnodes via
  // forgetMe interferes with the way the PHPointerList::clearAndDestroy()
  // works but it has to be executed in case the PHCompositeNode is
  // a parent and supposed to stay. Then the deleted node has to take itself
  // out of the node list
  deleteMe = 1;
  subNodes.clearAndDestroy();
}

PHBoolean
PHCompositeNode::addNode(PHNode* newNode)
{
  //
  // Check all existing subNodes for name-conflict.
  //
  PHPointerListIterator<PHNode> nodeIter(subNodes);
  PHNode* thisNode;
  while ((thisNode = nodeIter())) 
    {
      if (thisNode->getName() == newNode->getName()) 
	{
	  cout << PHWHERE << "Node " << newNode->getName() 
               << " already exists" << endl;
	  return False;
	}
    }
  //
  // No conflict, so we can append the new node.
  //
  newNode->setParent(this);
  return (subNodes.append(newNode));
}

void
PHCompositeNode::prune()
{
  PHPointerListIterator<PHNode> nodeIter(subNodes);
  PHNode* thisNode;
  while ((thisNode = nodeIter())) 
    {
      if (!thisNode->isPersistent()) 
	{
	  subNodes.removeAt(nodeIter.pos());
	  --nodeIter;
	  delete thisNode;
	}
      else
	{
	  thisNode->prune();
	}
    }
}

void
PHCompositeNode::forgetMe(PHNode* child)
{
// if this PHCompositeNode is supposed to be deleted,
// do not remove the child from the list,
// otherwise the clearanddestroy() bookkeeping gets
// confused and deletes only every other node
  if (deleteMe) 
    {
      return;
    }
  PHPointerListIterator<PHNode> nodeIter(subNodes);
  PHNode* thisNode;
  while (child && (thisNode = nodeIter())) 
    {
      if (thisNode == child) 
	{
	  subNodes.removeAt(nodeIter.pos());
	  child = 0;
	}
    }   
}

PHBoolean
PHCompositeNode::write(PHIOManager * IOManager, const PHString& path)
{
   PHString newPath;
   if (path != "")
     {
       newPath = path + "/" + name;
     }
   else
     {
       newPath = name;
     }

   PHPointerListIterator<PHNode> nodeIter(subNodes);
   PHNode* thisNode;
   PHBoolean success = True;
   while ((thisNode = nodeIter())) 
     {
       if (!(thisNode->write(IOManager, newPath)))
	 {
	   success = False;
	 }
     }
   return success;
}

void
PHCompositeNode::print(const PHString& path)
{
  PHString newPath = "   " + path;
  cout << path << name << " (" << type << ")/" << endl;
  PHPointerListIterator<PHNode> nodeIter(subNodes);
  PHNode* thisNode;
  while ((thisNode = nodeIter())) 
    {
      thisNode->print(newPath);
    }
}
