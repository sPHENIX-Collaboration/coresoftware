#ifndef TRACKRECO_PHACTSVERTEXFINDER_H
#define TRACKRECO_PHACTSVERTEXFINDER_H

#include "PHInitVertexing.h"



#include <trackbase/TrkrDefs.h>

class PHActsVertexFinder: public PHInitVertexing 
{
  
 public:
  PHActsVertexFinder(const std::string &name);
  virtual ~PHActsVertexFinder() {}
  
 protected:
  int Setup(PHCompositeNode *topNode) override;
  int Process(PHCompositeNode *topNode) override;

 private:
  
  
};

#endif // TRACKRECO_PHACTSVERTEXFINDER_H
