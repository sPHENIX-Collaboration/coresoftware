// Base CLUSTER class
// Stores 3D point related reconstructed cluster for central arm
// Origin: Carlos Perez

#include "vCluster.h"

//==========
vCluster::vCluster()
{
  fSgn = 0.0;
  for(int i=0; i!=3; ++i) fX[i] = 0.0;
  for(int i=0; i!=6; ++i) fCov[i] = 0.0;
  for(int i=0; i!=3; ++i) fSize[i] = 0;
}
//==========
vCluster::~vCluster()
{
}
//==========
void vCluster::CopyFrom(vCluster *vc)
{
  fSgn = vc->fSgn;
  for(int i=0; i!=3; ++i) fX[i] = vc->fX[i];
  for(int i=0; i!=6; ++i) fCov[i] = vc->fCov[i];
  for(int i=0; i!=3; ++i) fSize[i] = vc->fSize[i];
}
