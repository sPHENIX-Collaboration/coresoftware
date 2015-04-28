#include "Fun4AllDstOutputManager.h"
#include "Fun4AllServer.h"

#include <phool/PHNode.h>
#include <phool/PHNodeIOManager.h>
#include <phool/PHNodeIterator.h>

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

Fun4AllDstOutputManager::Fun4AllDstOutputManager(const string &myname, const string &fname): 
 Fun4AllOutputManager( myname )
{
  outfilename = fname;
  dstOut = new PHNodeIOManager(fname.c_str(), PHWrite);
  if (!dstOut->isFunctional())
    {
      delete dstOut;
      cout << PHWHERE << " Could not open " << fname 
	   << " exiting now" << endl;
      exit(1);
    }
  dstOut->SetCompressionLevel(3);
  return ;
}

Fun4AllDstOutputManager::~Fun4AllDstOutputManager()
{
  delete dstOut;
  return ;
}

int
Fun4AllDstOutputManager::AddNode(const string &nodename)
{
  string newnode = nodename;
  vector<string>::const_iterator iter;
  for (iter = savenodes.begin(); iter != savenodes.end(); ++iter)
    {
      if ( *iter == newnode)
        {
          cout << "Node " << newnode << " allready in list" << endl;
          return -1;
        }
    }
  savenodes.push_back(newnode);
  return 0;
}

int
Fun4AllDstOutputManager::StripNode(const string &nodename)
{
  string newnode = nodename;
  vector<string>::const_iterator iter;
  for (iter = stripnodes.begin(); iter != stripnodes.end(); ++iter)
    {
      if ( *iter == newnode)
        {
          cout << "Node " << newnode << " allready in list" << endl;
          return -1;
        }
    }
  stripnodes.push_back(newnode);
  return 0;
}

int
Fun4AllDstOutputManager::outfileopen(const string &fname)
{
  dstOut = new PHNodeIOManager(fname.c_str(), PHWrite);
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

int
Fun4AllDstOutputManager::RemoveNode(const string &nodename)
{
  string node = nodename;
  vector<string>::iterator iter;
  for (iter = savenodes.begin(); iter != savenodes.end(); ++iter)
  if ( *iter == node) {
    savenodes.erase(iter);
    cout << "Removing " << node << " from list" << endl;
    return 0;
  }
  
  cout << "Could not find " << node << " in list" << endl;
  
  return -1;

}

void
Fun4AllDstOutputManager::Print(const string &what) const
{
  if (what == "ALL" || what == "WRITENODES")
    {
      vector<string>::const_iterator iter;
      cout << ThisName << " writes " << outfilename << endl;
      if (savenodes.empty())
        {
          if (stripnodes.empty())
            {
              cout << ThisName << ": All Nodes will be written out" << endl;
            }
          else
            {
              for (iter = stripnodes.begin(); iter != stripnodes.end(); ++iter)
                {
                  cout << ThisName << ": Node " << *iter << " will be stripped" << endl;
                }
            }
        }
      else
        {
          for (iter = savenodes.begin(); iter != savenodes.end(); ++iter)
            {
              cout << ThisName << ": Node " << *iter << " is written out" << endl;
            }
        }
    }
  // base class print method
  Fun4AllOutputManager::Print( what );

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
int 
Fun4AllDstOutputManager::Write(PHCompositeNode *startNode)
{
  PHNodeIterator nodeiter(startNode);
  vector<string>::iterator iter;
  PHNode *ChosenNode = 0;
  if (savenodes.empty())
    {
      Fun4AllServer *se = Fun4AllServer::instance();
      se->MakeNodesPersistent(startNode);
      if (! stripnodes.empty())
	{
	  for (iter = stripnodes.begin(); iter != stripnodes.end(); ++iter)
	    {
	      ChosenNode = nodeiter.findFirst("PHIODataNode", iter->c_str());
	      if (ChosenNode)
		{
		  ChosenNode->makeTransient();
		}
	      else
		{
		  if (verbosity > 0)
		    {
		      cout << PHWHERE << ThisName << ": Node " << *iter
			   << " does not exist" << endl;
		    }
		}

	    }
	}
    }
  else
    {
      for (iter = savenodes.begin(); iter != savenodes.end(); ++iter)
        {
          ChosenNode = nodeiter.findFirst("PHIODataNode", iter->c_str());
          if (ChosenNode)
            {
              ChosenNode->makePersistent();
            }
	  else
	    {
	      if (verbosity > 0)
		{
		  cout << PHWHERE << ThisName << ": Node " << *iter 
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
      for (iter = savenodes.begin(); iter != savenodes.end(); ++iter)
        {
          ChosenNode = nodeiter.findFirst("PHIODataNode", iter->c_str());
          if (ChosenNode)
            {
              ChosenNode->makeTransient();
            }
        }
    }
  return 0;
}

int
Fun4AllDstOutputManager::WriteNode(PHCompositeNode *thisNode)
{
  delete dstOut;

  dstOut = new PHNodeIOManager(outfilename.c_str(), PHUpdate, PHRunTree);
  dstOut->write(thisNode);
  delete dstOut;
  dstOut = 0;
  return 0;
}

