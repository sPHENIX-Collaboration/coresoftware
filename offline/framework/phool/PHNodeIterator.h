#ifndef PHOOL_PHNODEITERATOR_H
#define PHOOL_PHNODEITERATOR_H

//  Declaration of class PHNodeIterator
//  Purpose: iterator to navigate a node tree
//  Author: Matthias Messer

//#include "PHCompositeNode.h"
#include "PHPointerList.h"

#include <string> 

class PHCompositeNode;
class PHNode;
class PHNodeOperation;

class PHNodeIterator
{
 public:
  explicit PHNodeIterator(PHCompositeNode*);
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

#endif
