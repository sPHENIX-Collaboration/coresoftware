#ifndef PHNodeIntegrate_h
#define PHNodeIntegrate_h

//  Declaration of class PHNodeIntegrate
//  Purpose: strategy which calls Integrate() on a PHNode if it is
//  a PHObject

#include "PHNodeOperation.h"

class PHNode;
class PHCompositeNode;

class PHNodeIntegrate : public PHNodeOperation
{
 public:
  PHNodeIntegrate():runnode(nullptr),runsumnode(nullptr) {}
  virtual ~PHNodeIntegrate() {}
  void RunNode(PHCompositeNode *node)
  {
    runnode = node;
  }
  void RunSumNode(PHCompositeNode *node)
  {
    runsumnode = node;
  }

 protected:
  virtual void perform(PHNode *);

 private:
  PHCompositeNode *runnode;
  PHCompositeNode *runsumnode;
};

#endif /* PHNodeIntegrate_h */
