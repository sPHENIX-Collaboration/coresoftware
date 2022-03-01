// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4CYLINDERCELLGEOMCONTAINER_H
#define G4DETECTORS_PHG4CYLINDERCELLGEOMCONTAINER_H

#include <phool/PHObject.h>

#include <iostream>  // for cout, ostream
#include <map>
#include <utility>  // for make_pair, pair

class PHG4CylinderCellGeom;

class PHG4CylinderCellGeomContainer : public PHObject
{
 public:
  typedef std::map<int, PHG4CylinderCellGeom *> Map;
  typedef Map::iterator Iterator;
  typedef Map::const_iterator ConstIterator;
  typedef std::pair<Iterator, Iterator> Range;
  typedef std::pair<ConstIterator, ConstIterator> ConstRange;

  PHG4CylinderCellGeomContainer() {}
  ~PHG4CylinderCellGeomContainer() override;

  // from PHObject
  void identify(std::ostream &os = std::cout) const override;

  int AddLayerCellGeom(const int i, PHG4CylinderCellGeom *mygeom);
  int AddLayerCellGeom(PHG4CylinderCellGeom *mygeom);
  PHG4CylinderCellGeom *GetLayerCellGeom(const int i);
  PHG4CylinderCellGeom *GetFirstLayerCellGeom();
  int get_NLayers() const { return layergeoms.size(); }
  std::pair<std::map<int, PHG4CylinderCellGeom *>::const_iterator, std::map<int, PHG4CylinderCellGeom *>::const_iterator> get_begin_end() const { return std::make_pair(layergeoms.begin(), layergeoms.end()); }

 protected:
  std::map<int, PHG4CylinderCellGeom *> layergeoms;
  ClassDefOverride(PHG4CylinderCellGeomContainer, 1)
};

#endif
