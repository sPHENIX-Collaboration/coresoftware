// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4TPCGEOMCONTAINER_H
#define G4DETECTORS_PHG4TPCGEOMCONTAINER_H

#include <phool/PHObject.h>

#include <iostream>  // for cout, ostream
#include <map>
#include <utility>  // for make_pair, pair

class PHG4TpcGeom;

class PHG4TpcGeomContainer : public PHObject
{
 public:
  typedef std::map<int, PHG4TpcGeom *> Map;
  typedef Map::iterator Iterator;
  typedef Map::const_iterator ConstIterator;
  typedef std::pair<Iterator, Iterator> Range;
  typedef std::pair<ConstIterator, ConstIterator> ConstRange;

  PHG4TpcGeomContainer() {}
  ~PHG4TpcGeomContainer() override;

  // from PHObject
  void identify(std::ostream &os = std::cout) const override;

  int AddLayerCellGeom(const int i, PHG4TpcGeom *mygeom);
  int AddLayerCellGeom(PHG4TpcGeom *mygeom);
  PHG4TpcGeom *GetLayerCellGeom(const int i);
  PHG4TpcGeom *GetFirstLayerCellGeom();
  int get_NLayers() const { return layergeoms.size(); }
  ConstRange get_begin_end() const { return {layergeoms.begin(), layergeoms.end()}; }

 protected:
  Map layergeoms;
  ClassDefOverride(PHG4TpcGeomContainer, 1)
};

#endif
