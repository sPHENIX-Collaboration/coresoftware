// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4BLOCKGEOMCONTAINER_H
#define G4DETECTORS_PHG4BLOCKGEOMCONTAINER_H

#include <phool/PHObject.h>

#include <iostream>  // for cout, ostream
#include <map>
#include <utility>  // for make_pair, pair

class PHG4BlockGeom;

class PHG4BlockGeomContainer : public PHObject
{
 public:
  typedef std::map<int, PHG4BlockGeom *> Map;
  typedef Map::iterator Iterator;
  typedef Map::const_iterator ConstIterator;
  typedef std::pair<Iterator, Iterator> Range;
  typedef std::pair<ConstIterator, ConstIterator> ConstRange;

  PHG4BlockGeomContainer();
  ~PHG4BlockGeomContainer() override;

  // from PHObject
  void identify(std::ostream &os = std::cout) const override;

  int AddLayerGeom(const int i, PHG4BlockGeom *mygeom);
  int AddLayerGeom(PHG4BlockGeom *mygeom);
  PHG4BlockGeom *GetLayerGeom(const int i);
  int get_NLayers() const { return _layergeoms.size(); }
  std::pair<std::map<int, PHG4BlockGeom *>::const_iterator, std::map<int, PHG4BlockGeom *>::const_iterator> get_begin_end() const { return std::make_pair(_layergeoms.begin(), _layergeoms.end()); }

 protected:
  std::map<int, PHG4BlockGeom *> _layergeoms;
  float _magfield;

  ClassDefOverride(PHG4BlockGeomContainer, 1)
};

#endif
