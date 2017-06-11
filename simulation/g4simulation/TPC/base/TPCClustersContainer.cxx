// TPC CLUSTERSCONTAINER class
// Stores TPC CLUSTERS per layer and sector
// notice that it stores a clusters in vector and not in map
// no need of index access since tracking will inspect per layer and sector only
// i.e. stores compressed info
// features avoidance of ctor dtor calls
// Author: Carlos Perez

#include "TPCClustersContainer.h"

//=====
TPCClustersContainer::TPCClustersContainer() :
{
  for(int i=0; i!=48; ++i)
    for(int j=0; j!=12; ++j)
      for(int k=0; k!=7; ++k)
	fNClusters[i][j][k] = 0;
}
//=====
TPCClustersContainer::~TPCClustersContainer() :
  for(int i=0; i!=48; ++i)
    for(int j=0; j!=12; ++j)
      for(int k=0; k!=7; ++k)
	for(int l=0; l!=fClusters[i][j][k].size(); ++l)
	  delete fClusters[i][j][k].at(l);
}
//=====
TPCClustersContainer::Reset() :
  for(int i=0; i!=48; ++i)
    for(int j=0; j!=12; ++j)
      for(int k=0; k!=7; ++k)
	fNClusters[i][j][k] = 0;
}
//=====
void Add(Int_t lyr, Int_t sec, Int_t blk, TPCCluster *clu) {
  //-- appending cluster into right container
  TPCCluster *mine;
  if(fNClusters[lyr][sec][blk]<fClusters[lyr][sec][blk].size()) {
    mine = fClusters[lyr][sec][blk].at(fNClusters[lyr][sec][blk]);
  } else {
    mine = new TPCCluster();
    fClusters[lyr][sec][blk].push_back( mine );
  }
  fNClusters[lyr][sec][blk]++;
  //-- coying info
  mine->CopyFrom( clu );
  // notice: clu is not owned by container!
}
