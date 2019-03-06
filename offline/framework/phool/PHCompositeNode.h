#ifndef PHOOL_PHCOMPOSITENODE_H
#define PHOOL_PHCOMPOSITENODE_H

//  Declaration of class PHCompositeNode
//  Purpose: a node which can hold other nodes

#include "PHNode.h"
#include "PHPointerList.h"

class PHIOManager;
class PHNodeIterator;

class PHCompositeNode : public PHNode
{
  friend class PHNodeIterator;

 public:
  explicit PHCompositeNode(const std::string &);
  virtual ~PHCompositeNode();

  //
  // The user is only allowed to add new nodes, not to delete existing ones.
  //
  bool addNode(PHNode *);

  //
  // This recursively calls the prune function of all the subnodes.
  // If a subnode is found to be marked as transient (non persistent)
  // the entire sub-tree is deleted.
  //
  virtual void prune();

  //
  // I/O functions
  //
  void print(const std::string & = "");
  virtual bool write(PHIOManager *, const std::string & = "");

 protected:
  virtual void forgetMe(PHNode *);
  PHPointerList<PHNode> subNodes;
  int deleteMe;

 private:
  PHCompositeNode() = delete;
};

#endif
