#ifndef __PHNODEITERATOR_H__
#define __PHNODEITERATOR_H__

//  Declaration of class PHNodeIterator
//  Purpose: iterator to navigate a node tree
//  Author: Matthias Messer

#include "phool.h"
#include "PHCompositeNode.h"
#include "PHPointerList.h"
#include "PHString.h"

class PHNodeOperation;

class PHNodeIterator { 

public: 
   PHNodeIterator(PHCompositeNode*); 
   virtual ~PHNodeIterator() {} 
   PHNodeIterator(); 
   
public:
   void                   print();
   PHPointerList<PHNode>& ls();
   PHNode*                findFirst(const PHString&, const PHString&);
   PHNode*                findFirst(const PHString&);
   PHBoolean              cd(const PHString = "");
   PHBoolean              addNode(PHNode*);
   void                   forEach(PHNodeOperation&);
   void                   for_each(PHNodeOperation&);
   
protected: 
   PHCompositeNode *currentNode;
   PHPointerList<PHNode> subNodeList;

}; 

#endif /* __PHNODEITERATOR_H__ */
