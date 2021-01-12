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
  std::set<TrkrDefs::cluskey>  discard_set;

  unsigned int count_discards = 0;

  // loop over all TPC clusters
  if(Verbosity() > 0) 
    std::cout << std::endl << "original size of cluster map: " << _cluster_map->size() << std::endl;  

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
      
      if (Verbosity() >= 1)
	{
	  std::cout << " cluster : " << cluskey << " layer " << layer
		    << " position x,y,z " << cluster->getX() << "  " << cluster->getY() << "  " << cluster->getZ()
		      << " ADC " << cluster->getAdc()
		    << std::endl;
	  std::cout << "       errors: r-phi " << cluster->getRPhiError() << " Z " << cluster->getZError() 
		    << " phi size " << cluster->getPhiSize() << " Z size " << cluster->getZSize()
		    << std::endl;
	}
      
      bool discard_cluster = false;
      
      // We have a TPC cluster, look for reasons to discard it
	
      // errors too small
      // associated with very large ADC values
      if(cluster->getRPhiError() < _rphi_error_low_cut)
	discard_cluster = true;
      
      // errors too large
      // associated with very small ADC values
      if(cluster->getRPhiError() > _rphi_error_high_cut)
	discard_cluster = true;
      
      if(discard_cluster)
	{
	  count_discards++;
	  // mark it for removal
	  discard_set.insert(cluskey);
	  if(Verbosity() > 0) 
	    std::cout << " will remove cluster " << cluskey << " with ephi " << cluster->getRPhiError() << " adc " << cluster->getAdc() 
		    << " phisize " << cluster->getPhiSize() << " Z size " << cluster->getZSize() << std::endl;
	}
    }
  
  for(auto iter = discard_set.begin(); iter != discard_set.end(); ++iter)
    {
      // remove bad clusters from the node tree map
      _cluster_map->removeCluster(*iter);
    }

  if(Verbosity() > 0)
    std::cout << "Clusters discarded this event: " << count_discards << std::endl;

  /*
  // check the map on the node tree 
  std::cout << " Readback modified cluster map - size is: " << _cluster_map->size() << std::endl;  
  clusRange = _cluster_map->getClusters();
  
  for (auto clusiter = clusRange.first; 
       clusiter != clusRange.second; ++clusiter)
    {
      TrkrDefs::cluskey cluskey = clusiter->first;
      TrkrCluster *cluster = clusiter->second;
      
      unsigned int trkrId = TrkrDefs::getTrkrId(cluskey);
      unsigned int layer = TrkrDefs::getLayer(cluskey);
      
      if(trkrId != TrkrDefs::tpcId) continue;  // we want only TPC clusters
      
      if (Verbosity() >= 1)
	{
	  std::cout << " cluster : " << cluskey << " layer " << layer
		    << " position x,y,z " << cluster->getX() << "  " << cluster->getY() << "  " << cluster->getZ()
		      << " ADC " << cluster->getAdc()
		    << std::endl;
	  //std::cout << "       errors: r-phi " << cluster->getRPhiError() << " Z " << cluster->getZError() 
	  //	      << " phi size " << cluster->getPhiSize() << " Z size " << cluster->getZSize()
	  //	      << std::endl;
	}
    }
  */
  
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

