/**
 * @file mvtx/MvtxClusterPruner.cc
 * @author Hugo Pereira Da Costa <hugo.pereira-da-costa@lanl.gov>
 * @date May 2025
 * @brief Implementation of MvtxClusterPruner
 */

#include "MvtxClusterPruner.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>

#include <trackbase/MvtxDefs.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterHitAssoc.h>

#include <set>
#include <algorithm>

namespace
{
  //! range adaptor to be able to use range-based for loop
  template<class T> class range_adaptor
  {
    public:
    range_adaptor( const T& range ):m_range(range){}
    const typename T::first_type& begin() {return m_range.first;}
    const typename T::second_type& end() {return m_range.second;}
    private:
    T m_range;
  };
}

//_____________________________________________________________________________
MvtxClusterPruner::MvtxClusterPruner(const std::string &name)
  : SubsysReco(name)
{
}

//_____________________________________________________________________________
int MvtxClusterPruner::InitRun(PHCompositeNode * /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________________
int MvtxClusterPruner::process_event(PHCompositeNode *topNode)
{
  // load relevant nodes
  auto trkrclusters = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if( !trkrclusters )
  {
    std::cout << "MvtxClusterPruner::process_event - TRKR_CLUSTER not found. Doing nothing" << std::endl;
    return Fun4AllReturnCodes::EVENT_OK;
  }

  auto clusterhitassoc = findNode::getClass<TrkrClusterHitAssoc>(topNode, "TRKR_CLUSTERHITASSOC");
  if( !clusterhitassoc )
  {
    std::cout << "MvtxClusterPruner::process_event - TRKR_CLUSTERHITASSOC not found. Doing nothing" << std::endl;
    return Fun4AllReturnCodes::EVENT_OK;
  }

  // loop over MVTX hitset keys
  const auto hitsetkeys = trkrclusters->getHitSetKeys(TrkrDefs::mvtxId);
  for( const auto& hitsetkey:hitsetkeys )
  {
    // get layer, stave, chip and current strobe
    const auto layer = TrkrDefs::getLayer(hitsetkey);
    const auto stave = MvtxDefs::getStaveId(hitsetkey);
    const auto chip = MvtxDefs::getChipId(hitsetkey);
    const auto current_strobe = MvtxDefs::getStrobeId(hitsetkey);

    // get clusters for this hitsetkey
    const auto cluster_range1= trkrclusters->getClusters(hitsetkey);

    // get clusters for the next strobe
    int next_strobe = current_strobe+1;
    const auto hitsetkey_next_strobe = MvtxDefs::genHitSetKey(layer, stave, chip, next_strobe);
    const auto cluster_range2 = trkrclusters->getClusters(hitsetkey_next_strobe);

    // loop over clusters from first range
    for( const auto& [ckey1,cluster1]:range_adaptor(cluster_range1) )
    {

      // get associated hits
      const auto& hit_range1 = clusterhitassoc->getHits(ckey1);
      std::set<TrkrDefs::hitkey> hitkeys1;
      std::transform(hit_range1.first, hit_range1.second, std::inserter(hitkeys1,hitkeys1.end()),
        [](const TrkrClusterHitAssoc::Map::value_type& pair ){ return pair.second; });

      // loop over clusters from second range
      for( const auto& [ckey2,cluster2]:range_adaptor(cluster_range2) )
      {
        // get associated hits
        const auto& hit_range2 = clusterhitassoc->getHits(ckey2);
        std::set<TrkrDefs::hitkey> hitkeys2;
        std::transform(hit_range2.first, hit_range2.second, std::inserter(hitkeys2,hitkeys2.end()),
          [](const TrkrClusterHitAssoc::Map::value_type& pair ){ return pair.second; });

        // make sure first set is larger than second
        const bool swapped = hitkeys2.size() > hitkeys1.size();
        if( swapped ) { std::swap(hitkeys2,hitkeys1); }

        // see if hitkeys2 is a subset of hitkeys1
        if( std::includes(hitkeys1.begin(), hitkeys1.end(), hitkeys2.begin(), hitkeys2.end()) )
        {
          if( swapped )
          {
            // remove first cluster
            trkrclusters->removeCluster(ckey1);
            break;
          } else {
            // remove second cluster
            trkrclusters->removeCluster(ckey2);
          }
        }
      } // second cluster loop
    } // first cluster loop
  } // hitsetkey loop

  return Fun4AllReturnCodes::EVENT_OK;

}

//_____________________________________________________________________________
int MvtxClusterPruner::End(PHCompositeNode * /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}
