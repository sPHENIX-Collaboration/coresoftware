// TPC HIT class
// Stores  one hit in TPC fiducial volume
// Author: Carlos Perez
#include "TPCHit.h"

//=====
TPCHit::TPCHit() : fDEnergy(0.0) {
}
//=====
TPCHit::~TPCHit() {}
//=====
void TPCHit::CopyFrom(TPCHit *th) {
  vHit::CopyFrom(th);
  fDEnergy = th->fDEnergy;
}
