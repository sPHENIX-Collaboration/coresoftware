//  Implementation of class PHRawDataNode
//  Author: Matthias Messer

#include "PHRawDataNode.h"

#include "PHRawOManager.h"

#include <phool/PHDataNode.h>   // for PHDataNode
#include <phool/PHIOManager.h>

#include <Event/phenixTypes.h>

PHRawDataNode::PHRawDataNode(PHDWORD* d, const std::string& n,
                             const int l, const int i, const int w, const int h)
  : PHDataNode<PHDWORD>(d, n)
  , length(l)
  , ID(i)
  , wordLength(w)
  , hitFormat(h)
{
}

PHRawDataNode::~PHRawDataNode()
{
  // set the data poitner to 0 so
  // the dtor of the PHDataNode parent class doesn't try
  // to delete it
  setData(nullptr);
}

bool PHRawDataNode::write(PHIOManager* IOManager, const std::string&)
{
  PHRawOManager* rawOManager = dynamic_cast<PHRawOManager*>(IOManager);
  bool bret = false;
  if (rawOManager)
  {
    bret = rawOManager->write(this);
  }

  return bret;
}
