#ifndef TRACKRECO_ASSOCINFOCONTAINER_H
#define TRACKRECO_ASSOCINFOCONTAINER_H

#include <phool/PHObject.h>
#include <trackbase/TrkrDefs.h>

#include <map>

class AssocInfoContainer : public PHObject
{
 public:
  typedef std::multimap<TrkrDefs::cluskey, unsigned int> ClusterTrackMap;

  AssocInfoContainer();
  virtual ~AssocInfoContainer();

  void Reset();
  void identify(std::ostream& os = std::cout) const;

  void SetClusterTrackAssoc(const TrkrDefs::cluskey& cluster_id, const unsigned int& track_id)
  {
    _map_cluster_id_track_id.insert(ClusterTrackMap::value_type(cluster_id, track_id));
  }

  std::vector<unsigned int> GetTracksFromCluster(const TrkrDefs::cluskey& cluster_id) const
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

  ClassDef(AssocInfoContainer, 1)
};

#endif
