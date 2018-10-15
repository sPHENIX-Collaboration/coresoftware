#include "Fun4AllDstOutputManager.h"

#include "Fun4AllServer.h"

#include <phool/PHNode.h>
#include <phool/PHNodeIOManager.h>
#include <phool/PHNodeIterator.h>

#include <boost/foreach.hpp>

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

Fun4AllDstOutputManager::Fun4AllDstOutputManager(const string &myname, const string &fname)
  : Fun4AllOutputManager(myname)
{
  outfilename = fname;
  dstOut = new PHNodeIOManager(fname, PHWrite);
  if (!dstOut->isFunctional())
  {
    delete dstOut;
    cout << PHWHERE << " Could not open " << fname
         << " exiting now" << endl;
    exit(1);
  }
  dstOut->SetCompressionLevel(3);
  return;
}

Fun4AllDstOutputManager::~Fun4AllDstOutputManager()
{
  delete dstOut;
  return;
}

int Fun4AllDstOutputManager::AddNode(const string &nodename)
{
  savenodes.insert(nodename);
  return 0;
}

int Fun4AllDstOutputManager::StripNode(const string &nodename)
{
  stripnodes.insert(nodename);
  return 0;
}

int Fun4AllDstOutputManager::StripRunNode(const string &nodename)
{
  striprunnodes.insert(nodename);
  return 0;
}

int Fun4AllDstOutputManager::outfileopen(const string &fname)
{
  dstOut = new PHNodeIOManager(fname, PHWrite);
  if (!dstOut->isFunctional())
  {
    delete dstOut;
    dstOut = 0;
    cout << PHWHERE << " Could not open " << fname << endl;
    return -1;
  }

  dstOut->SetCompressionLevel(3);
  return 0;
}

void Fun4AllDstOutputManager::Print(const string &what) const
{
  if (what == "ALL" || what == "WRITENODES")
  {
//    vector<string>::const_iterator iter;
    cout << Name() << " writes " << outfilename << endl;
    if (savenodes.empty())
    {
      if (stripnodes.empty())
      {
        cout << Name() << ": All Nodes will be written out" << endl;
      }
      else
      {
        BOOST_FOREACH (string nodename, stripnodes)
        {
          cout << Name() << ": Node " << nodename << " will be stripped" << endl;
        }
      }
    }
    else
    {
      BOOST_FOREACH (string nodename, savenodes)
      {
        cout << Name() << ": Node " << nodename << " is written out" << endl;
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
//  vector<string>::iterator iter;
  PHNode *ChosenNode = 0;
  if (savenodes.empty())
  {
    Fun4AllServer *se = Fun4AllServer::instance();
    se->MakeNodesPersistent(startNode);
    if (!stripnodes.empty())
    {
      BOOST_FOREACH (string nodename, stripnodes)
      {
        ChosenNode = nodeiter.findFirst("PHIODataNode", nodename);
        if (ChosenNode)
        {
          ChosenNode->makeTransient();
        }
        else
        {
          if (Verbosity() > 0)
          {
            cout << PHWHERE << Name() << ": Node " << nodename
                 << " does not exist" << endl;
          }
        }
      }
    }
  }
  else
  {
    BOOST_FOREACH (string nodename, savenodes)
    {
      ChosenNode = nodeiter.findFirst("PHIODataNode", nodename);
      if (ChosenNode)
      {
        ChosenNode->makePersistent();
      }
      else
      {
        if (Verbosity() > 0)
        {
          cout << PHWHERE << Name() << ": Node " << nodename
               << " does not exist" << endl;
        }
      }
    }
  }
  dstOut->write(startNode);
  if (savenodes.empty())
  {
    Fun4AllServer *se = Fun4AllServer::instance();
    se->MakeNodesTransient(startNode);
  }
  else
  {
    BOOST_FOREACH (string nodename, savenodes)
    {
      ChosenNode = nodeiter.findFirst("PHIODataNode", nodename);
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
  dstOut = new PHNodeIOManager(outfilename, PHUpdate, PHRunTree);
  Fun4AllServer *se = Fun4AllServer::instance();
  se->MakeNodesPersistent(thisNode);
  if (!striprunnodes.empty())
  {
    PHNodeIterator nodeiter(thisNode);
    BOOST_FOREACH (string nodename, striprunnodes)
    {
      PHNode *ChosenNode = nodeiter.findFirst("PHIODataNode", nodename);
      if (ChosenNode)
      {
        ChosenNode->makeTransient();
      }
    }
  }
  dstOut->write(thisNode);
  se->MakeNodesPersistent(thisNode);
  delete dstOut;
  dstOut = nullptr;
  return 0;
}
