// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4BLOCKCELLGEOMCONTAINER_H
#define G4DETECTORS_PHG4BLOCKCELLGEOMCONTAINER_H

#include <phool/PHObject.h>

#include <iostream>  // for cout, ostream
#include <map>
#include <utility>  // for make_pair, pair

class PHG4BlockCellGeom;

class PHG4BlockCellGeomContainer : public PHObject
{
 public:
  typedef std::map<int, PHG4BlockCellGeom *> Map;
  typedef Map::iterator Iterator;
  typedef Map::const_iterator ConstIterator;
  typedef std::pair<Iterator, Iterator> Range;
  typedef std::pair<ConstIterator, ConstIterator> ConstRange;

  PHG4BlockCellGeomContainer() {}
  ~PHG4BlockCellGeomContainer() override;

  // from PHObject
  void identify(std::ostream &os = std::cout) const override;

  int AddLayerCellGeom(const int i, PHG4BlockCellGeom *mygeom);
  int AddLayerCellGeom(PHG4BlockCellGeom *mygeom);
  PHG4BlockCellGeom *GetLayerCellGeom(const int i);
  int get_NLayers() const { return layergeoms.size(); }
  std::pair<std::map<int, PHG4BlockCellGeom *>::const_iterator, std::map<int, PHG4BlockCellGeom *>::const_iterator> get_begin_end() const { return std::make_pair(layergeoms.begin(), layergeoms.end()); }

 protected:
  std::map<int, PHG4BlockCellGeom *> layergeoms;
  ClassDefOverride(PHG4BlockCellGeomContainer, 1)
};

#endif
