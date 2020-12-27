// Tell emacs that this is a C++ source
//  -*- C++ -*-.

/*!
 *  \file		  TpcClusterCleaner.h
 *  \brief		 Class for removing bad TPC clusters
 *  \author	 Tony Frawley <afrawley@fsu.edu>
 */

#ifndef TPCCLUSTERCLEANER_H
#define TPCCLUSTERCLEANER_H

#include <fun4all/SubsysReco.h>

#include <string>
#include <vector>

class PHCompositeNode;
class TrkrCluster;
class TrkrClusterContainer;

class TpcClusterCleaner : public SubsysReco
{
 public:

  TpcClusterCleaner(const std::string &name = "TpcClusterCleaner");

  virtual ~TpcClusterCleaner();

  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);


 private:

  int GetNodes(PHCompositeNode* topNode);

  TrkrClusterContainer *_cluster_map;


};

#endif // TPCCLUSTERCLEANER_H
