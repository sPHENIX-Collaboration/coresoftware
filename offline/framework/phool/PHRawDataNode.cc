//  Implementation of class PHRawDataNode
//  Author: Matthias Messer

#include "PHRawDataNode.h"
#include "PHRawOManager.h"

PHRawDataNode::PHRawDataNode() 
{
   length     = 0;
   ID         = 0;
   wordLength = 0;
   hitFormat  = 0;
}

PHRawDataNode::~PHRawDataNode()
{
  // set the data poitner to 0 so
  // the dtor of the PHDataNode parent class doesn't try
  // to delete it
  setData(0);
}

PHRawDataNode::PHRawDataNode(PHDWORD * d, const PHString& n, 
			     const int l, const int i, const int w, const int h) 
  : PHDataNode<PHDWORD>(d, n)
{
   length     = l;
   ID         = i;
   wordLength = w;
   hitFormat  = h;
}

PHBoolean 
PHRawDataNode::write(PHIOManager * IOManager, const PHString&)
{
  PHRawOManager* rawOManager = dynamic_cast<PHRawOManager*>(IOManager);

  if (rawOManager) 
    {
      rawOManager->write(this);
    }
  
  return True;
}
