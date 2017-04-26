// TPC HITSCONTAINER class
// Stores TPC HITS per block
// features avoidance of ctor dtor calls
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
  TPCHit* GetHit(Int_t blk, UInt_t at);
  UInt_t GetNHits(Int_t blk) {return fNHits[blk];}
  void Reset();

 protected:
  UInt_t fNHits[2];
  std::vector<TPCHit*> fHits[2];
};

#endif
