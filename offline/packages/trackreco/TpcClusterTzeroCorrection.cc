#include "TpcClusterTzeroCorrection.h"

/// Tracking includes

#include <trackbase/TrkrCluster.h>            // for TrkrCluster
#include <trackbase/TrkrClusterContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/phool.h>

#include <iostream>                           // for operator<<, basic_ostream
#include <map>                                // for map
#include <set>                                // for _Rb_tree_const_iterator
#include <utility>                            // for pair, make_pair
#include <cassert>                          


//____________________________________________________________________________..
TpcClusterTzeroCorrection::TpcClusterTzeroCorrection(const std::string &name)
  : SubsysReco(name)
  , PHParameterInterface(name)
{ InitializeParameters(); }


//____________________________________________________________________________..
int TpcClusterTzeroCorrection::InitRun(PHCompositeNode* /*unused*/)
{
  UpdateParametersWithMacro();
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int TpcClusterTzeroCorrection::process_event(PHCompositeNode *topNode )
{
  // load nodes
  const auto res = load_nodes(topNode);
  if( res != Fun4AllReturnCodes::EVENT_OK ) { return res; }

  process_clusters();
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int TpcClusterTzeroCorrection::End(PHCompositeNode* /*topNode*/ )
{ return Fun4AllReturnCodes::EVENT_OK; }

//_____________________________________________________________________
void TpcClusterTzeroCorrection::SetDefaultParameters()
{
  // Data on gasses @20 C and 760 Torr from the following source:
  // http://www.slac.stanford.edu/pubs/icfa/summer98/paper3/paper3.pdf
  // diffusion and drift velocity for 400kV for NeCF4 50/50 from calculations:
  // http://skipper.physics.sunysb.edu/~prakhar/tpc/HTML_Gases/split.html
  return;
}

//_____________________________________________________________________
int TpcClusterTzeroCorrection::load_nodes( PHCompositeNode* topNode )
{
  m_cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, m_clusterContainerName);
  assert(m_cluster_map);

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
void TpcClusterTzeroCorrection::process_clusters()
{
  if( !m_cluster_map ) { return; }

 for (const auto& hitsetkey : m_cluster_map->getHitSetKeys(TrkrDefs::TrkrId::tpcId))
  {
    auto range = m_cluster_map->getClusters(hitsetkey);
    for (auto clusIter = range.first; clusIter != range.second; ++clusIter)
      {
	TrkrDefs::cluskey cluster_key = clusIter->first;
	
	// consider TPC clusters only
	if( TrkrDefs::getTrkrId(cluster_key) != TrkrDefs::tpcId ) { continue;}

	TrkrCluster* cluster = clusIter->second;	
	float initialLocalY = cluster->getLocalY();
	cluster->setLocalY( cluster->getLocalY() - m_tzero);

	if( Verbosity() )
	  { 
	    std::cout << "TpcClusterTzeroCorrection::process_cluster: " << cluster_key
	    << " initial localY: " << initialLocalY << " m_tzero " << m_tzero 
	    << " corrected localY " << cluster->getLocalY() << std::endl;
	  }
      }
  }
}
