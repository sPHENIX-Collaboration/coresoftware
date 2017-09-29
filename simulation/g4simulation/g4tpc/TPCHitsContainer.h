// Organizes TPC HITS per block
// Author: Carlos Perez

#ifndef __TPCHITSCONTAINER_H__
#define __TPCHITSCONTAINER_H__

#include <vector>
#include <phool/PHObject.h>
#include "TPCHit.h"

class TPCHitsContainer : public PHObject {
 public:
  TPCHitsContainer();
  virtual ~TPCHitsContainer();
  void Add(TPCHit *hit);
  TPCHit* GetHit(int blk, unsigned int at);
  unsigned int GetNHits(Int_t blk) {return fNHits[blk];}
  void Reset();

 protected:
  unsigned int fNHits[2];
  std::vector<TPCHit*> fHits[2];

  ClassDef(TPCHitsContainer,1);
};

#endif
