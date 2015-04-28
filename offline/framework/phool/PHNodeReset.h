#ifndef PHNODERESET_H__
#define PHNODERESET_H__

//  Declaration of class PHNodeReset
//  Purpose: strategy which calls reset() on a PHNode 
//  Author: Matthias Messer

#include "PHNodeOperation.h"

class PHNode;

class PHNodeReset : public PHNodeOperation 
{ 
public: 
  PHNodeReset(){} 
   virtual ~PHNodeReset(){} 

protected: 
   virtual void perform(PHNode*);

}; 

#endif /* PHNODERESET_H__ */
