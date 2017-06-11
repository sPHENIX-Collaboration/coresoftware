// TPC HITSCONTAINER class
// Stores TPC HITS per block
// features avoidance of ctor dtor calls
// Author: Carlos Perez

#include <vector>
#include "TPCHitsContainer.h"

//=====
TPCHitsContainer::TPCHitsContainer()
{
  Reset();
}
//=====
TPCHitsContainer::~TPCHitsContainer()
{
  for(unsigned int i=0; i!=fHits[0].size(); ++i) delete fHits[0].at(i);
  for(unsigned int i=0; i!=fHits[1].size(); ++i) delete fHits[1].at(i);
}
//=====
void TPCHitsContainer::Reset()
{
  fNHits[0] = 0;
  fNHits[1] = 0;
}
//=====
void TPCHitsContainer::Add(TPCHit *hit) {
  Int_t blk = hit->GetZ()<0? 0: 1;
  //-- appending hit into right container
  TPCHit *mine;
  if(fNHits[blk]<fHits[blk].size()) {
    mine = fHits[blk].at(fNHits[blk]);
  } else {
    mine = new TPCHit();
    fHits[blk].push_back( mine );
  }
  fNHits[blk]++;
  //-- coying info
  mine->CopyFrom( hit );
  // notice: hit is not owned by container!
}
//=====
TPCHit* TPCHitsContainer::GetHit(Int_t blk, UInt_t at) {
  if(fNHits[blk]<at) return NULL;
  return fHits[blk].at(at);
}
