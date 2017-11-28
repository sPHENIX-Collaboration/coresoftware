#ifndef PHNODEITERATOR_H__
#define PHNODEITERATOR_H__

//  Declaration of class PHNodeIterator
//  Purpose: iterator to navigate a node tree
//  Author: Matthias Messer

#include "PHCompositeNode.h"
#include "PHPointerList.h"

class PHNodeOperation;

class PHNodeIterator
{
 public:
  PHNodeIterator(PHCompositeNode*);
  virtual ~PHNodeIterator() {}
  PHNodeIterator();

 public:
  void print();
  PHPointerList<PHNode>& ls();
  PHNode* findFirst(const std::string&, const std::string&);
  PHNode* findFirst(const std::string&);
  bool cd(const std::string& pathString = "");
  bool addNode(PHNode*);
  void forEach(PHNodeOperation&);
  void for_each(PHNodeOperation&);
  PHCompositeNode* get_currentNode() const { return currentNode; }
 protected:
  PHCompositeNode* currentNode;
  PHPointerList<PHNode> subNodeList;
};

#endif /* __PHNODEITERATOR_H__ */
