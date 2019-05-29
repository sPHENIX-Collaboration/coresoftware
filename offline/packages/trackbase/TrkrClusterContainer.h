/**
 * @file trackbase/TrkrClusterContainer.h
 * @author D. McGlinchey
 * @date June 2018
 * @brief Cluster container object
 */
#ifndef TRACKBASE_TRKRCLUSTERCONTAINER_H
#define TRACKBASE_TRKRCLUSTERCONTAINER_H

#include "TrkrDefs.h"

#include <phool/PHObject.h>

#include <map>
#include <iostream>          // for cout, ostream
#include <utility>           // for pair

class TrkrCluster;

/**
 * @brief Cluster container object
 *
 * Container for TrkrCluster objects
 */
class TrkrClusterContainer : public PHObject
{
 public:
  typedef std::map<TrkrDefs::cluskey, TrkrCluster *> Map;
  typedef Map::iterator Iterator;
  typedef Map::const_iterator ConstIterator;
  typedef std::pair<Iterator, Iterator> Range;
  typedef std::pair<ConstIterator, ConstIterator> ConstRange;

  TrkrClusterContainer(){}

  virtual ~TrkrClusterContainer() {}
  void Reset();

  void identify(std::ostream &os = std::cout) const;

  ConstIterator addCluster(TrkrCluster *newClus);
  ConstIterator addClusterSpecifyKey(const TrkrDefs::cluskey key, TrkrCluster *newClus);

  //! preferred removal method, key is currently the clus id
  void removeCluster(TrkrDefs::cluskey key)
  {
    m_clusmap.erase(key);
  }

  //! inefficent, use key where possible instead
  void removeCluster(TrkrCluster *clus)
  {
    Iterator its = m_clusmap.begin();
    Iterator last = m_clusmap.end();
    for (; its != last;)
    {
      if (its->second == clus)
      {
        m_clusmap.erase(its++);
      }
      else
      {
        ++its;
      }
    }
  }

  Iterator findOrAddCluster(TrkrDefs::cluskey key);

  //! return all Clusters matching a given detid
  ConstRange getClusters(const TrkrDefs::TrkrId trackerid) const;

  //! return all Clusters matching a given detid and layer
  ConstRange getClusters(const TrkrDefs::TrkrId trackerid, const char layer) const;

  //! return all clusters
  ConstRange getClusters(void) const;

  TrkrCluster *findCluster(TrkrDefs::cluskey key);

  unsigned int size(void) const
  {
    return m_clusmap.size();
  }

 protected:
  Map m_clusmap;
  ClassDef(TrkrClusterContainer, 1)
};

#endif //TRACKBASE_TRKRCLUSTERCONTAINER_H
