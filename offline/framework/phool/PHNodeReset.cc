//  Implementation of class PHNodeReset
//  Author: Matthias Messer

#include "PHNodeReset.h" 
#include "PHDataNode.h"
#include "PHIODataNode.h" 
#include "PHObject.h"

void
PHNodeReset::perform(PHNode* node)
{
  if ( node->getResetFlag() != True ) return;

   if (node->getType() == "PHDataNode")
     {  
       if (node->getObjectType() == "PHObject")
	 {
	   (static_cast<PHDataNode<PHObject>*>(node))->getData()->Reset();
         }
     }
   else if (node->getType() == "PHIODataNode")
     {
       if (node->getObjectType() == "PHObject")
	 {
	   (static_cast<PHIODataNode<PHObject>*>(node))->getData()->Reset();
	 }
     }
}
