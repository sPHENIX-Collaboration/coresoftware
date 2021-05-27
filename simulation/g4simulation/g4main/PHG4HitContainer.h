// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4HITCONTAINER_H
#define G4MAIN_PHG4HITCONTAINER_H

#include "PHG4HitDefs.h"

#include <phool/PHObject.h>

#include <iostream>
#include <map>
#include <set>
#include <string>
#include <utility>

class PHG4Hit;

class PHG4HitContainer: public PHObject
{

  public:
  typedef std::map<PHG4HitDefs::keytype, PHG4Hit *> Map;
  typedef Map::iterator Iterator;
  typedef Map::const_iterator ConstIterator;
  typedef std::pair<Iterator, Iterator> Range;
  typedef std::pair<ConstIterator, ConstIterator> ConstRange;
  typedef std::set<unsigned int>::const_iterator LayerIter;

  PHG4HitContainer(); //< used only by ROOT for DST readback
  PHG4HitContainer(const std::string &nodename);

  ~PHG4HitContainer() override {}

  void Reset() override;

  void identify(std::ostream& os = std::cout) const override;

  //! container ID should follow definition of PHG4HitDefs::get_volume_id(DST nodename)
  void SetID(int i) {id = i;}
  int GetID() const {return id;}
  
  ConstIterator AddHit(PHG4Hit *newhit);

  ConstIterator AddHit(const unsigned int detid, PHG4Hit *newhit);
  
  Iterator findOrAddHit(PHG4HitDefs::keytype key);

  PHG4Hit* findHit(PHG4HitDefs::keytype key );

  PHG4HitDefs::keytype genkey(const unsigned int detid);

  //! return all hits matching a given detid
  ConstRange getHits(const unsigned int detid) const;

  //! return all hist
  ConstRange getHits( void ) const;

  unsigned int size( void ) const
  { return hitmap.size(); }
  unsigned int num_layers(void) const
  { return layers.size(); }
  std::pair<LayerIter, LayerIter> getLayers() const
     { return make_pair(layers.begin(), layers.end());} 
  void AddLayer(const unsigned int ilayer) {layers.insert(ilayer);}
  void RemoveZeroEDep();
  PHG4HitDefs::keytype getmaxkey(const unsigned int detid);

 protected:

  int id; //< unique identifier from hash of node name. Defined following PHG4HitDefs::get_volume_id
  Map hitmap;
  std::set<unsigned int> layers; // layers is not reset since layers must not change event by event

  ClassDefOverride(PHG4HitContainer,1)
};

#endif
