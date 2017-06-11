// TPC CLOUD class
// Stores  one cloud in TPC fiducial volume
// Author: Carlos Perez
#include <TPCbase/TPCHit.h>
#include "TPCCloud.h"

//=====
TPCCloud::TPCCloud() : 
  fMST(0.0),
  fMSL0(0.0),
  fMSL1(0.0)
{
}
//=====
TPCCloud::~TPCCloud()
{
}
//=====
void TPCCloud::CopyFrom(TPCHit *th)
{
  TPCHit::CopyFrom(th);
  fMST = 0.0;
  fMSL0 = 0.0;
  fMSL1 = 0.0;
}
