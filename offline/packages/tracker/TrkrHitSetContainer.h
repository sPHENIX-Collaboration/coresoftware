#ifndef __TrkrHitSetContainer_H__
#define __TrkrHitSetContainer_H__

#include "TrkrHitSet.h"

#include <phool/PHObject.h>

#include <map>
#include <set>

class TrkrHitSetContainer: public PHObject
{

  public:
  typedef std::map<TrkrDefs::hitsetkey,TrkrHitSet *> Map;
  typedef Map::iterator Iterator;
  typedef Map::const_iterator ConstIterator;
  typedef std::pair<Iterator, Iterator> Range;
  typedef std::pair<ConstIterator, ConstIterator> ConstRange;

  TrkrHitSetContainer();

  virtual ~TrkrHitSetContainer() {}

  void Reset();

  void identify(std::ostream& os = std::cout) const;

  ConstIterator AddHitSet(TrkrHitSet *newHit);
  ConstIterator AddHitSetSpecifyKey(const TrkrDefs::hitsetkey key, TrkrHitSet *newHit);
  
  //! preferred removal method, key is currently the hit id
  void RemoveHitSet(TrkrDefs::hitsetkey key) {
    hitmap.erase(key);
  }

  //! inefficent, use key where possible instead
  void RemoveHitSet(TrkrHitSet *hit)
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


  Iterator findOrAddHitSet(TrkrDefs::hitsetkey key);

  //! return all Hits matching a given detid
  ConstRange getHitSets(const TrkrDefs::TRKRID trackerid) const;

  //! return all Hits matching a given detid, layer
  ConstRange getHitSets(const TrkrDefs::TRKRID trackerid, const char layer) const;

  //! return all hist
  ConstRange getHitSets( void ) const;

  TrkrHitSet* findHitSet(TrkrDefs::hitsetkey key);

  unsigned int size( void ) const
  { return hitmap.size(); }

 protected:
  Map hitmap;
  ClassDef(TrkrHitSetContainer,1)
};

#endif
