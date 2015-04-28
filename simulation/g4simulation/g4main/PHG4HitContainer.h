#ifndef __PHG4HITCONTAINER_H__
#define __PHG4HITCONTAINER_H__


#include <phool/PHObject.h>
#include <map>
#include <set>
class PHG4Hit;

class PHG4HitContainer: public PHObject
{

  public:
  typedef std::map<unsigned long long,PHG4Hit *> Map;
  typedef Map::iterator Iterator;
  typedef Map::const_iterator ConstIterator;
  typedef std::pair<Iterator, Iterator> Range;
  typedef std::pair<ConstIterator, ConstIterator> ConstRange;
  typedef std::set<int>::const_iterator LayerIter;

  PHG4HitContainer();

  virtual ~PHG4HitContainer() {}

  void Reset();

  void identify(std::ostream& os = std::cout) const;

  ConstIterator AddHit(const int detid, PHG4Hit *newhit);
  
  Iterator findOrAddHit(unsigned long long key);

  PHG4Hit* findHit( unsigned long long key );

  unsigned long long genkey(const int detid);

  //! return all hits matching a given detid
  ConstRange getHits(const int detid) const;

  //! return all hist
  ConstRange getHits( void ) const;

  unsigned int size( void ) const
  { return hitmap.size(); }
  unsigned int num_layers(void) const
  { return layers.size(); }
  std::pair<LayerIter, LayerIter> getLayers() const
    { return make_pair(layers.begin(), layers.end());}
  void AddLayer(const int ilayer) {layers.insert(ilayer);}
  void RemoveZeroEDep();
  unsigned long long getmaxkey(const int detid);

 protected:
  Map hitmap;
  std::set<int> layers; // layers is not reset since layers must not change event by event

  ClassDef(PHG4HitContainer,1)
};

#endif
