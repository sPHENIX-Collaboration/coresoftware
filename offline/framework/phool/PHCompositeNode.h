#ifndef PHOOL_PHCOMPOSITENODE_H
#define PHOOL_PHCOMPOSITENODE_H

//  Declaration of class PHCompositeNode
//  Purpose: a node which can hold other nodes

#include "PHNode.h"
#include "PHPointerList.h"

#include <string>

class PHIOManager;

class PHCompositeNode : public PHNode
{
  friend class PHNodeIterator;

 public:
  explicit PHCompositeNode(const std::string &);
  ~PHCompositeNode() override;

  //
  // The user is only allowed to add new nodes, not to delete existing ones.
  //
  bool addNode(PHNode *);

  //
  // This recursively calls the prune function of all the subnodes.
  // If a subnode is found to be marked as transient (non persistent)
  // the entire sub-tree is deleted.
  //
  void prune() override;

  //
  // I/O functions
  //
  void print(const std::string & = "") override;
  bool write(PHIOManager *, const std::string & = "") override;

 protected:
  void forgetMe(PHNode *) override;
  PHPointerList<PHNode> subNodes;
  int deleteMe = 0;

 private:
  PHCompositeNode() = delete;
};

#endif
