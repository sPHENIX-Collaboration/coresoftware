#ifndef __PHCOMPOSITENODE_H__
#define __PHCOMPOSITENODE_H__

//  Declaration of class PHCompositeNode
//  Purpose: a node which can hold other nodes

#include "phool.h"
#include "PHNode.h"
#include "PHPointerList.h"

class PHIOManager;
class PHNodeIterator;

class PHCompositeNode : public PHNode { 

   friend class PHNodeIterator;
   
public: 
   PHCompositeNode(const PHString &); 
   virtual ~PHCompositeNode(); 

   //
   // The user is only allowed to add new nodes, not to delete existing ones.
   //
   PHBoolean addNode(PHNode*);

   //
   // This recursively calls the prune function of all the subnodes.
   // If a subnode is found to be marked as transient (non persistent)
   // the entire sub-tree is deleted.
   //
   virtual void prune();

   //
   // I/O functions
   //
   void print(const PHString& = "");
   virtual PHBoolean write(PHIOManager *, const PHString& = "");

protected:
   virtual void forgetMe(PHNode*);
   PHPointerList<PHNode> subNodes;
   int deleteMe;

private:
   PHCompositeNode();
}; 

#endif /* __PHCOMPOSITENODE_H__ */
