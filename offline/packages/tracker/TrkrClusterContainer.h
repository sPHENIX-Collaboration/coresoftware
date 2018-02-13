#ifndef __TrkrClusterContainer_H__
#define __TrkrClusterContainer_H__

#include "TrkrCluster.h"
#include "TrkrDefUtil.h"

#include <phool/PHObject.h>

#include <map>
#include <set>

class TrkrClusterContainer : public PHObject
{
 public:
  typedef std::map<TrkrDefs::cluskey, TrkrCluster *> Map;
  typedef Map::iterator Iterator;
  typedef Map::const_iterator ConstIterator;
  typedef std::pair<Iterator, Iterator> Range;
  typedef std::pair<ConstIterator, ConstIterator> ConstRange;

  TrkrClusterContainer();

  virtual ~TrkrClusterContainer() {}
  void Reset();

  void identify(std::ostream &os = std::cout) const;

  ConstIterator AddCluster(TrkrCluster *newClus);
  ConstIterator AddClusterSpecifyKey(const TrkrDefs::cluskey key, TrkrCluster *newClus);

  //! preferred removal method, key is currently the clus id
  void RemoveCluster(TrkrDefs::cluskey key)
  {
    clusmap_.erase(key);
  }

  //! inefficent, use key where possible instead
  void RemoveCluster(TrkrCluster *clus)
  {
    Iterator its = clusmap_.begin();
    Iterator last = clusmap_.end();
    for (; its != last;)
    {
      if (its->second == clus)
      {
        clusmap_.erase(its++);
      }
      else
      {
        ++its;
      }
    }
  }

  Iterator FindOrAddCluster(TrkrDefs::cluskey key);

  //! return all Clusters matching a given detid
  ConstRange GetClusters(const TrkrDefs::TRKRID trackerid) const;

  //! return all Clusters matching a given detid and layer
  ConstRange GetClusters(const TrkrDefs::TRKRID trackerid, const char layer) const;

  //! return all clusters
  ConstRange GetClusters(void) const;

  TrkrCluster *FindCluster(TrkrDefs::cluskey key);

  unsigned int size(void) const
  {
    return clusmap_.size();
  }

 protected:
  Map clusmap_;
  ClassDef(TrkrClusterContainer, 1)
};

#endif
