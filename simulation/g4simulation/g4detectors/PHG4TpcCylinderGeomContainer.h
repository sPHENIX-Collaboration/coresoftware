// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4TPCCYLINDERGEOMCONTAINER_H
#define G4DETECTORS_PHG4TPCCYLINDERGEOMCONTAINER_H

#include <phool/PHObject.h>

#include <iostream>  // for cout, ostream
#include <map>
#include <utility>  // for make_pair, pair

class PHG4TpcCylinderGeom;

class PHG4TpcCylinderGeomContainer : public PHObject
{
 public:
  typedef std::map<int, PHG4TpcCylinderGeom *> Map;
  typedef Map::iterator Iterator;
  typedef Map::const_iterator ConstIterator;
  typedef std::pair<Iterator, Iterator> Range;
  typedef std::pair<ConstIterator, ConstIterator> ConstRange;

  PHG4TpcCylinderGeomContainer() {}
  ~PHG4TpcCylinderGeomContainer() override;

  // from PHObject
  void identify(std::ostream &os = std::cout) const override;

  int AddLayerCellGeom(const int i, PHG4TpcCylinderGeom *mygeom);
  int AddLayerCellGeom(PHG4TpcCylinderGeom *mygeom);
  PHG4TpcCylinderGeom *GetLayerCellGeom(const int i);
  PHG4TpcCylinderGeom *GetFirstLayerCellGeom();
  int get_NLayers() const { return layergeoms.size(); }
  std::pair<std::map<int, PHG4TpcCylinderGeom *>::const_iterator, std::map<int, PHG4TpcCylinderGeom *>::const_iterator> get_begin_end() const { return std::make_pair(layergeoms.begin(), layergeoms.end()); }

 protected:
  std::map<int, PHG4TpcCylinderGeom *> layergeoms;
  ClassDefOverride(PHG4TpcCylinderGeomContainer, 1)
};

#endif
