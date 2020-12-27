#include "TpcClusterCleaner.h"   

/// Tracking includes

#include <trackbase/TrkrCluster.h>            // for TrkrCluster
#include <trackbase/TrkrDefs.h>               // for cluskey, getLayer, TrkrId
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterHitAssoc.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/phool.h>

#include <cmath>                              // for sqrt, fabs, atan2, cos
#include <iostream>                           // for operator<<, basic_ostream
#include <map>                                // for map
#include <set>                                // for _Rb_tree_const_iterator
#include <utility>                            // for pair, make_pair

//____________________________________________________________________________..
TpcClusterCleaner::TpcClusterCleaner(const std::string &name)
  : SubsysReco(name)
{

}

//____________________________________________________________________________..
TpcClusterCleaner::~TpcClusterCleaner()
{

}

//____________________________________________________________________________..
int TpcClusterCleaner::InitRun(PHCompositeNode *topNode)
{
  int ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  return ret;
}

//____________________________________________________________________________..
int TpcClusterCleaner::process_event(PHCompositeNode *topNode)
{

  // loop over all TPC clusters

  TrkrClusterContainer::ConstRange clusRange = _cluster_map->getClusters();
  TrkrClusterContainer::ConstIterator clusiter;
  
  for (clusiter = clusRange.first; 
       clusiter != clusRange.second; ++clusiter)
    {
      TrkrDefs::cluskey cluskey = clusiter->first;
      TrkrCluster *cluster = clusiter->second;
      
      unsigned int trkrId = TrkrDefs::getTrkrId(cluskey);
      unsigned int layer = TrkrDefs::getLayer(cluskey);
      
      if(trkrId != TrkrDefs::tpcId) continue;  // we want only TPC clusters
 
      //if (Verbosity() >= 1)
	{
	  std::cout << " cluster : " << cluskey << " layer " << layer
		    << " position x,y,z " << cluster->getX() << "  " << cluster->getY() << "  " << cluster->getZ()
		    << std::endl;
	    std::cout << "       errors: r-phi " << cluster->getRPhiError() << " Z " << cluster->getZError() 
		    << " ADC " << cluster->getAdc()
		    << " phi size " << cluster->getPhiSize() << " Z size " << cluster->getZSize()
		    << std::endl;
	}



      // We have a TPC cluster
      // Look for reasons to discard it

      // Energy too large to be from a primary particle


      // size too large to be from a primary particle


      // errors too small





  

    }
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcClusterCleaner::End(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int  TpcClusterCleaner::GetNodes(PHCompositeNode* topNode)
{
  _cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!_cluster_map)
  {
    std::cout << PHWHERE << " ERROR: Can't find node TRKR_CLUSTER" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

