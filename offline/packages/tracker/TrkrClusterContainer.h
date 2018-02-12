#ifndef __TrkrClusterContainer_H__
#define __TrkrClusterContainer_H__

#include "TrkrCluster.h"
#include "TrackerDefs.h"

#include <phool/PHObject.h>

#include <map>
#include <set>

class TrkrClusterContainer: public PHObject
{

  public:
  typedef std::map<TrackerDefs::keytype,TrkrCluster *> Map;
  typedef Map::iterator Iterator;
  typedef Map::const_iterator ConstIterator;
  typedef std::pair<Iterator, Iterator> Range;
  typedef std::pair<ConstIterator, ConstIterator> ConstRange;

  TrkrClusterContainer();

  virtual ~TrkrClusterContainer() {}

  void Reset();

  void identify(std::ostream& os = std::cout) const;

  ConstIterator AddCluster(TrkrCluster *newClus);
  ConstIterator AddClusterSpecifyKey(const TrackerDefs::keytype key, TrkrCluster *newClus);
  
  //! preferred removal method, key is currently the clus id
  void RemoveCluster(TrackerDefs::keytype key) {
    clusmap.erase(key);
  }

  //! inefficent, use key where possible instead
  void RemoveCluster(TrkrCluster *clus)
  {
    Iterator its = clusmap.begin();
    Iterator last = clusmap.end();
    for (; its != last;)
      {
	if (its->second == clus)
	  {
	    clusmap.erase(its++);
	  }
	else
	  {
	    ++its;
	  }
      }
  }


  Iterator findOrAddCluster(TrackerDefs::keytype key);

  //! return all Clusters matching a given detid
  ConstRange getClusters(const TrackerDefs::TRACKERID trackerid) const;

  //! return all clusters
  ConstRange getClusters( void ) const;

  TrkrCluster* findCluster(TrackerDefs::keytype key);

  unsigned int size( void ) const
  { return clusmap.size(); }

 protected:
  Map clusmap;
  ClassDef(TrkrClusterContainer,1)
};

#endif
