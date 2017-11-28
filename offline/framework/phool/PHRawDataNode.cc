//  Implementation of class PHRawDataNode
//  Author: Matthias Messer

#include "PHRawDataNode.h"
#include "PHRawOManager.h"

PHRawDataNode::PHRawDataNode()
  : length(0)
  , ID(0)
  , wordLength(0)
  , hitFormat(0)
{
}

PHRawDataNode::~PHRawDataNode()
{
  // set the data poitner to 0 so
  // the dtor of the PHDataNode parent class doesn't try
  // to delete it
  setData(nullptr);
}

PHRawDataNode::PHRawDataNode(PHDWORD* d, const string& n,
                             const int l, const int i, const int w, const int h)
  : PHDataNode<PHDWORD>(d, n)
  , length(l)
  , ID(i)
  , wordLength(w)
  , hitFormat(h)
{
}

bool PHRawDataNode::write(PHIOManager* IOManager, const string&)
{
  PHRawOManager* rawOManager = dynamic_cast<PHRawOManager*>(IOManager);
  bool bret = false;
  if (rawOManager)
  {
    bret = rawOManager->write(this);
  }

  return bret;
}
