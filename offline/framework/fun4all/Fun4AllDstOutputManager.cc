#include "Fun4AllDstOutputManager.h"

#include "Fun4AllServer.h"

#include <phool/PHNode.h>
#include <phool/PHNodeIOManager.h>
#include <phool/PHNodeIterator.h>
#include <phool/phool.h>  // for PHWHERE, PHReadOnly, PHRunTree

#include <TSystem.h>

#include <cstdlib>
#include <iostream>
#include <string>

Fun4AllDstOutputManager::Fun4AllDstOutputManager(const std::string &myname, const std::string &fname)
  : Fun4AllOutputManager(myname, fname)
{
  dstOut = new PHNodeIOManager(fname, PHWrite);
  if (!dstOut->isFunctional())
  {
    delete dstOut;
    std::cout << PHWHERE << " Could not open " << fname
              << " exiting now" << std::endl;
    gSystem->Exit(1);
    exit(1);  // cppcheck does not know gSystem->Exit(1)
  }
  dstOut->SetCompressionLevel(3);
  return;
}

Fun4AllDstOutputManager::~Fun4AllDstOutputManager()
{
  delete dstOut;
  return;
}

int Fun4AllDstOutputManager::AddNode(const std::string &nodename)
{
  savenodes.insert(nodename);
  return 0;
}

int Fun4AllDstOutputManager::AddRunNode(const std::string &nodename)
{
  saverunnodes.insert(nodename);
  return 0;
}

int Fun4AllDstOutputManager::StripNode(const std::string &nodename)
{
  stripnodes.insert(nodename);
  return 0;
}

int Fun4AllDstOutputManager::StripRunNode(const std::string &nodename)
{
  striprunnodes.insert(nodename);
  return 0;
}

int Fun4AllDstOutputManager::outfileopen(const std::string &fname)
{
  delete dstOut;
  dstOut = new PHNodeIOManager(fname, PHWrite);
  if (!dstOut->isFunctional())
  {
    delete dstOut;
    dstOut = nullptr;
    std::cout << PHWHERE << " Could not open " << fname << std::endl;
    return -1;
  }

  dstOut->SetCompressionLevel(3);
  return 0;
}

void Fun4AllDstOutputManager::Print(const std::string &what) const
{
  if (what == "ALL" || what == "WRITENODES")
  {
    std::cout << Name() << " writes " << OutFileName() << std::endl;
    if (savenodes.empty())
    {
      if (stripnodes.empty())
      {
        std::cout << Name() << ": All Nodes will be written out" << std::endl;
      }
      else
      {
        for (auto &nodename : stripnodes)
        {
          std::cout << Name() << ": Node " << nodename << " will be stripped" << std::endl;
        }
      }
    }
    else
    {
      for (auto &nodename : savenodes)
      {
        std::cout << Name() << ": Node " << nodename << " is written out" << std::endl;
      }
    }
  }
  // base class print method
  Fun4AllOutputManager::Print(what);

  return;
}

// All nodes are set to transient by the framework
// here we first change the nodes we want to write out
// to persistent and then call the write method
// of the io manager
// afterwards the nodes we just wrote out are changed back
// to transient
// if we want to strip nodes (only meaningful if we take the default
// that everything is written out), those nodes are declared transient
int Fun4AllDstOutputManager::Write(PHCompositeNode *startNode)
{
  PHNodeIterator nodeiter(startNode);
  if (savenodes.empty())
  {
    Fun4AllServer *se = Fun4AllServer::instance();
    se->MakeNodesPersistent(startNode);
    if (!stripnodes.empty())
    {
      for (auto &nodename : stripnodes)
      {
        PHNode *ChosenNode = nodeiter.findFirst("PHIODataNode", nodename);
        if (ChosenNode)
        {
          ChosenNode->makeTransient();
        }
        else
        {
          if (Verbosity() > 0)
          {
            std::cout << PHWHERE << Name() << ": Node " << nodename
                      << " does not exist" << std::endl;
          }
        }
      }
    }
  }
  else
  {
    for (auto &nodename : savenodes)
    {
      PHNode *ChosenNode = nodeiter.findFirst("PHIODataNode", nodename);
      if (ChosenNode)
      {
        ChosenNode->makePersistent();
      }
      else
      {
        if (Verbosity() > 0)
        {
          std::cout << PHWHERE << Name() << ": Node " << nodename
                    << " does not exist" << std::endl;
        }
      }
    }
  }
  dstOut->write(startNode);
  // to save some cpu cycles we only make it globally transient if
  // all nodes have been written (savenodes set is empty)
  // else we only make the nodes transient which we have written (all
  // others are transient by construction)
  if (savenodes.empty())
  {
    Fun4AllServer *se = Fun4AllServer::instance();
    se->MakeNodesTransient(startNode);
  }
  else
  {
    for (auto &nodename : savenodes)
    {
      PHNode *ChosenNode = nodeiter.findFirst("PHIODataNode", nodename);
      if (ChosenNode)
      {
        ChosenNode->makeTransient();
      }
    }
  }
  return 0;
}

int Fun4AllDstOutputManager::WriteNode(PHCompositeNode *thisNode)
{
  delete dstOut;
  if (!m_SaveRunNodeFlag)
  {
    dstOut = nullptr;
    return 0;
  }
  dstOut = new PHNodeIOManager(OutFileName(), PHUpdate, PHRunTree);
  Fun4AllServer *se = Fun4AllServer::instance();
  PHNodeIterator nodeiter(thisNode);
  if (saverunnodes.empty())
  {
    se->MakeNodesPersistent(thisNode);
    if (!striprunnodes.empty())
    {
      for (auto &nodename : striprunnodes)
      {
        PHNode *ChosenNode = nodeiter.findFirst("PHIODataNode", nodename);
        if (ChosenNode)
        {
          ChosenNode->makeTransient();
        }
        else
        {
          if (Verbosity() > 0)
          {
            std::cout << PHWHERE << Name() << ": Node " << nodename
                      << " does not exist" << std::endl;
          }
        }
      }
    }
  }
  else
  {
    for (auto &nodename : saverunnodes)
    {
      PHNode *ChosenNode = nodeiter.findFirst("PHIODataNode", nodename);
      if (ChosenNode)
      {
        ChosenNode->makePersistent();
      }
      else
      {
        if (Verbosity() > 0)
        {
          std::cout << PHWHERE << Name() << ": Node " << nodename
                    << " does not exist" << std::endl;
        }
      }
    }
  }
  dstOut->write(thisNode);
  se->MakeNodesTransient(thisNode);
  delete dstOut;
  dstOut = nullptr;
  return 0;
}
