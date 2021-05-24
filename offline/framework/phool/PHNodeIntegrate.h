#ifndef PHOOL_PHNODEINTEGRATE_H
#define PHOOL_PHNODEINTEGRATE_H

//  Declaration of class PHNodeIntegrate
//  Purpose: strategy which calls Integrate() on a PHNode if it is
//  a PHObject

#include "PHNodeOperation.h"

class PHNode;
class PHCompositeNode;

class PHNodeIntegrate : public PHNodeOperation
{
 public:
  PHNodeIntegrate()
    : runnode(nullptr)
    , runsumnode(nullptr)
  {
  }
  ~PHNodeIntegrate() override {}
  void RunNode(PHCompositeNode *node)
  {
    runnode = node;
  }
  void RunSumNode(PHCompositeNode *node)
  {
    runsumnode = node;
  }

 protected:
  void perform(PHNode *) override;

 private:
  PHCompositeNode *runnode;
  PHCompositeNode *runsumnode;
};

#endif
