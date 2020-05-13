/*!
 *  \file       PHTpcTrackerUtil.cc
 *  \brief      
 *  \author     Dmitry Arkhipkin <arkhipkin@gmail.com>
 */

#include "PHTpcTrackerUtil.h"

#include <trackbase/TrkrClusterContainer.h>  // for TrkrClusterContainer
#include <trackbase/TrkrCluster.h>

#include <phool/PHLog.h>

#include <cstdint>                          // for uint64_t
#include <cstring>
#include <map>                               // for _Rb_tree_const_iterator
#include <utility>                           // for pair

namespace PHTpcTrackerUtil
{
  std::vector<std::vector<double> > convert_clusters_to_hits(TrkrClusterContainer* cluster_map)
  {
    //***** convert clusters to kdhits
    std::vector<std::vector<double> > kdhits;
    if (!cluster_map)
    {
      LOG_WARN("tracking.PHTpcTrackerUtil.convert_clusters_to_hits") << "cluster map is not provided";
      return kdhits;
    }
    TrkrClusterContainer::ConstRange clusrange = cluster_map->getClusters();
    for (TrkrClusterContainer::ConstIterator it = clusrange.first; it != clusrange.second; ++it)
    {
      // TrkrDefs::cluskey cluskey = it->first;
      TrkrCluster* cluster = it->second;
      if ((std::pow((double) cluster->getPosition(0), 2) +
           std::pow((double) cluster->getPosition(1), 2) +
           std::pow((double) cluster->getPosition(2), 2)) > (25.0 * 25.0))
      {
        std::vector<double> kdhit(4);
        kdhit[0] = cluster->getPosition(0);
        kdhit[1] = cluster->getPosition(1);
        kdhit[2] = cluster->getPosition(2);
        uint64_t key = cluster->getClusKey();
        std::memcpy(&kdhit[3], &key, sizeof(key));

        //	HINT: way to get original uint64_t value from double:
        //
        //			LOG_DEBUG("tracking.PHTpcTrackerUtil.convert_clusters_to_hits")
        //				<< "orig: " << cluster->getClusKey() << ", readback: " << (*((int64_t*)&kdhit[3]));

        kdhits.push_back(kdhit);
      }
    }
    return kdhits;
  }

}  // namespace PHTpcTrackerUtil
