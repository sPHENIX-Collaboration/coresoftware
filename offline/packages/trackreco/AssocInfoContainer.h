#ifndef _H_AssocInfoContainer_H_
#define _H_AssocInfoContainer_H_

#include <phool/PHObject.h>

#include <map>

class AssocInfoContainer : public PHObject
{
 public:
  typedef std::multimap<unsigned int, unsigned int> ClusterTrackMap;

  AssocInfoContainer();
  virtual ~AssocInfoContainer();

  void Reset();
  void identify(std::ostream& os = std::cout) const;

  void SetClusterTrackAssoc(const unsigned int& cluster_id, const unsigned int& track_id)
  {
    _map_cluster_id_track_id.insert(ClusterTrackMap::value_type(cluster_id, track_id));
  }

  std::vector<unsigned int> GetTracksFromCluster(const unsigned int& cluster_id) const
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

#endif  //_H_AssocInfoContainer_H_
