#ifndef TRACKRECO_ASSOCINFOCONTAINERV1_H
#define TRACKRECO_ASSOCINFOCONTAINERV1_H

#include "AssocInfoContainer.h"

#include <trackbase/TrkrDefs.h>

#include <iostream>  // for cout, ostream
#include <map>
#include <utility>  // for pair
#include <vector>   // for vector

class AssocInfoContainerv1 : public AssocInfoContainer
{
 public:
  typedef std::multimap<TrkrDefs::cluskey, unsigned int> ClusterTrackMap;

  AssocInfoContainerv1();
  ~AssocInfoContainerv1() override;

  void Reset() override;
  void identify(std::ostream& os = std::cout) const override;

  void SetClusterTrackAssoc(const TrkrDefs::cluskey& cluster_id, const unsigned int& track_id) override
  {
    _map_cluster_id_track_id.insert(ClusterTrackMap::value_type(cluster_id, track_id));
  }

  std::vector<unsigned int> GetTracksFromCluster(const TrkrDefs::cluskey& cluster_id) const override
  {
    std::vector<unsigned int> ret;
    for (auto iter = _map_cluster_id_track_id.lower_bound(cluster_id);
         iter != _map_cluster_id_track_id.upper_bound(cluster_id); ++iter)
    {
      ret.push_back(iter->second);
    }
    return ret;
  }

 private:
  ClusterTrackMap _map_cluster_id_track_id;

  ClassDefOverride(AssocInfoContainerv1, 1)
};

#endif
