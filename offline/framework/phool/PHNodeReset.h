#ifndef PHNodeReset_h
#define PHNodeReset_h

//  Declaration of class PHNodeReset
//  Purpose: strategy which calls reset() on a PHNode
//  Author: Matthias Messer

#include "PHNodeOperation.h"

class PHNode;

class PHNodeReset : public PHNodeOperation
{
 public:
  PHNodeReset() {}
  virtual ~PHNodeReset() {}
 protected:
  virtual void perform(PHNode*);
};

#endif /* PHNodeReset_h */
