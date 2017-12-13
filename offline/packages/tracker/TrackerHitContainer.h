#ifndef TRACKERHITCONTAINER_H__
#define TRACKERHITCONTAINER_H__

#include "TrackerHit.h"

#include <phool/PHObject.h>

#include <map>
#include <set>

class TrackerHitContainer: public PHObject
{

  public:
  typedef std::map<TrackerDefs::keytype,TrackerHit *> Map;
  typedef Map::iterator Iterator;
  typedef Map::const_iterator ConstIterator;
  typedef std::pair<Iterator, Iterator> Range;
  typedef std::pair<ConstIterator, ConstIterator> ConstRange;

  TrackerHitContainer();

  virtual ~TrackerHitContainer() {}

  void Reset();

  void identify(std::ostream& os = std::cout) const;

  ConstIterator AddHit(TrackerHit *newHit);
  ConstIterator AddHitSpecifyKey(const TrackerDefs::keytype key, TrackerHit *newHit);
  
  //! preferred removal method, key is currently the hit id
  void RemoveHit(TrackerDefs::keytype key) {
    hitmap.erase(key);
  }

  //! inefficent, use key where possible instead
  void RemoveHit(TrackerHit *hit)
  {
    Iterator its = hitmap.begin();
    Iterator last = hitmap.end();
    for (; its != last;)
      {
	if (its->second == hit)
	  {
	    hitmap.erase(its++);
	  }
	else
	  {
	    ++its;
	  }
      }
  }


  Iterator findOrAddHit(TrackerDefs::keytype key);

  //! return all Hits matching a given detid
  ConstRange getHits(const TrackerDefs::TRACKERID trackerid) const;

  //! return all Hits matching a given detid, layer
  ConstRange getHits(const TrackerDefs::TRACKERID trackerid, const char layer) const;

  //! return all hist
  ConstRange getHits( void ) const;

  TrackerHit* findHit(TrackerDefs::keytype key);

  unsigned int size( void ) const
  { return hitmap.size(); }

 protected:
  Map hitmap;
  ClassDef(TrackerHitContainer,1)
};

#endif
