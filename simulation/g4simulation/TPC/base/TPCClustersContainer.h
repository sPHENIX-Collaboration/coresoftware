// TPC CLUSTERSCONTAINER class
// Stores TPC CLUSTERS per layer and sector
// notice that it stores a clusters in vector and not in map
// no need of index access since tracking will inspect per layer and sector only
// i.e. stores compressed info
// features avoidance of ctor dtor calls
// Author: Carlos Perez

#ifndef __TPCCLUSTERSCONTAINER_H__
#define __TPCCLUSTERSCONTAINER_H__

#include <vector.h>
#include <phool/PHObject.h>
#include "TPCCluster.h"

class TPCClustersContainer : public PHObject {
 public:
  TPCClustersContainer();
  virtual ~TPCClustersContainer();
  void Add(Int_t lyr, Int_t sec, Int_t blk, TPCCluster *clu);
  std::vector<TPCCluster*> GetClusters(Int_t lyr, Int_t sec, Int_t blk) {return fClusters[lyr][sec][blk];}
  Int_t GetNClusters(Int_t lyr, Int_t sec, Int_t blk) {return fNClusters[lyr][sec][blk];}
  void Reset();

 protected:
  UShort_t fNClusters[48][12][7];
  std::vector<vCluster*> fClusters[7][48][7];
  ClassDef(TPCClustersContainer,1);
};

#endif
