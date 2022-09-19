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
#include "phool.h"
#include "phooldefs.h"

#include <iostream>

PHCompositeNode::PHCompositeNode(const std::string& n)
  : PHNode(n, "PHCompositeNode")
{
 type = "PHCompositeNode";
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

bool PHCompositeNode::addNode(PHNode* newNode)
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
      std::cout << PHWHERE << "Node " << newNode->getName()
           << " already exists" << std::endl;
      return false;
    }
  }
  //
  // No conflict, so we can append the new node.
  //
  newNode->setParent(this);
  return (subNodes.append(newNode));
}

void PHCompositeNode::prune()
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

void PHCompositeNode::forgetMe(PHNode* child)
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
      child = nullptr;
    }
  }
}

bool PHCompositeNode::write(PHIOManager* IOManager, const std::string& path)
{
  std::string newPath = name;
  if (!path.empty())
  {
    newPath = path + phooldefs::branchpathdelim + name;
  }
  PHPointerListIterator<PHNode> nodeIter(subNodes);
  PHNode* thisNode;
  bool success = true;
  while ((thisNode = nodeIter()))
  {
    if (!(thisNode->write(IOManager, newPath)))
    {
      success = false;
    }
  }
  return success;
}

void PHCompositeNode::print(const std::string& path)
{
  std::string newPath = "   " + path;
  std::cout << path << name << " (" << type << ")/" << std::endl;
  PHPointerListIterator<PHNode> nodeIter(subNodes);
  PHNode* thisNode;
  while ((thisNode = nodeIter()))
  {
    thisNode->print(newPath);
  }
}
