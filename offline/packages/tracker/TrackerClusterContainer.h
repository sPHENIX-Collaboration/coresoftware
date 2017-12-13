#ifndef TRACKERCLUSTERCONTAINER_H__
#define TRACKERCLUSTERCONTAINER_H__

#include "TrackerCluster.h"
#include "TrackerDefs.h"

#include <phool/PHObject.h>

#include <map>
#include <set>

class TrackerClusterContainer: public PHObject
{

  public:
  typedef std::map<TrackerDefs::keytype,TrackerCluster *> Map;
  typedef Map::iterator Iterator;
  typedef Map::const_iterator ConstIterator;
  typedef std::pair<Iterator, Iterator> Range;
  typedef std::pair<ConstIterator, ConstIterator> ConstRange;

  TrackerClusterContainer();

  virtual ~TrackerClusterContainer() {}

  void Reset();

  void identify(std::ostream& os = std::cout) const;

  ConstIterator AddCluster(TrackerCluster *newClus);
  ConstIterator AddClusterSpecifyKey(const TrackerDefs::keytype key, TrackerCluster *newClus);
  
  //! preferred removal method, key is currently the clus id
  void RemoveCluster(TrackerDefs::keytype key) {
    clusmap.erase(key);
  }

  //! inefficent, use key where possible instead
  void RemoveCluster(TrackerCluster *clus)
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

  TrackerCluster* findCluster(TrackerDefs::keytype key);

  unsigned int size( void ) const
  { return clusmap.size(); }

 protected:
  Map clusmap;
  ClassDef(TrackerClusterContainer,1)
};

#endif
