//-----------------------------------------------------------------------------
//
//  The PHOOL's Software
//  Copyright (C) PHENIX collaboration, 1999
//
//  Implementation of class PHNodeIterator
//
//  Author: Matthias Messer
//-----------------------------------------------------------------------------
#include "PHNodeIterator.h"

#include "PHCompositeNode.h"
#include "PHNode.h"
#include "PHNodeOperation.h"
#include "PHPointerListIterator.h"
#include "phooldefs.h"

#include <boost/algorithm/string.hpp>

#include <vector>

PHNodeIterator::PHNodeIterator(PHCompositeNode* node)
  : currentNode(node)
{
}

PHPointerList<PHNode>&
PHNodeIterator::ls()
{
  PHPointerListIterator<PHNode> iter(currentNode->subNodes);
  subNodeList.clear();
  PHNode* thisNode;
  while ((thisNode = iter()))
  {
    subNodeList.append(thisNode);
  }
  return subNodeList;
}

void PHNodeIterator::print()
{
  currentNode->print();
}

// NOLINTNEXTLINE(misc-no-recursion)
PHNode* PHNodeIterator::findFirst(const std::string& requiredType, const std::string& requiredName)
{
  PHPointerListIterator<PHNode> iter(currentNode->subNodes);
  PHNode* thisNode;
  while ((thisNode = iter()))
  {
    if (thisNode->getType() == requiredType && thisNode->getName() == requiredName)
    {
      return thisNode;
    }

    if (thisNode->getType() == "PHCompositeNode")
    {
      PHNodeIterator nodeIter(dynamic_cast<PHCompositeNode*>(thisNode));
      PHNode* nodeFoundInSubTree = nodeIter.findFirst(requiredType, requiredName);
      if (nodeFoundInSubTree)
      {
        return nodeFoundInSubTree;
      }
    }
  }
  return nullptr;
}

// NOLINTNEXTLINE(misc-no-recursion)
PHNode* PHNodeIterator::findFirst(const std::string& requiredName)
{
  PHPointerListIterator<PHNode> iter(currentNode->subNodes);
  PHNode* thisNode;
  while ((thisNode = iter()))
  {
    if (thisNode->getName() == requiredName)
    {
      return thisNode;
    }

    if (thisNode->getType() == "PHCompositeNode")
    {
      PHNodeIterator nodeIter(dynamic_cast<PHCompositeNode*>(thisNode));
      PHNode* nodeFoundInSubTree = nodeIter.findFirst(requiredName);
      if (nodeFoundInSubTree)
      {
        return nodeFoundInSubTree;
      }
    }
  }
  return nullptr;
}

bool PHNodeIterator::cd(const std::string& pathString)
{
  bool success = true;
  if (pathString.empty())
  {
    while (currentNode->getParent())
    {
      currentNode = dynamic_cast<PHCompositeNode*>(currentNode->getParent());
    }
  }
  else
  {
    std::vector<std::string> splitpath;
    boost::split(splitpath, pathString, boost::is_any_of(phooldefs::nodetreepathdelim));
    bool pathFound;
    PHNode* subNode;
    for (const auto& iter : splitpath)
    {
      if (iter == "..")
      {
        if (currentNode->getParent())
        {
          currentNode = dynamic_cast<PHCompositeNode*>(currentNode->getParent());
        }
        else
        {
          success = false;
        }
      }
      else
      {
        PHPointerListIterator<PHNode> subNodeIter(currentNode->subNodes);
        pathFound = false;
        while ((subNode = subNodeIter()))
        {
          if (subNode->getType() == "PHCompositeNode" && subNode->getName() == iter)
          {
            currentNode = dynamic_cast<PHCompositeNode*>(subNode);
            pathFound = true;
          }
        }
        if (!pathFound)
        {
          success = false;
        }
      }
    }
  }
  return success;
}

bool PHNodeIterator::addNode(PHNode* newNode)
{
  return currentNode->addNode(newNode);
}

// NOLINTNEXTLINE(misc-no-recursion)
void PHNodeIterator::forEach(PHNodeOperation& operation)
{
  operation(currentNode);
  PHPointerListIterator<PHNode> iter(currentNode->subNodes);
  PHNode* thisNode;
  while ((thisNode = iter()))
  {
    if (thisNode->getType() == "PHCompositeNode")
    {
      PHNodeIterator subNodeIter(dynamic_cast<PHCompositeNode*>(thisNode));
      subNodeIter.forEach(operation);
    }
    else
    {
      operation(thisNode);
    }
  }
}

void PHNodeIterator::for_each(PHNodeOperation& operation)
{
  forEach(operation);
}
