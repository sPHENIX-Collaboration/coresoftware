#ifndef PHG4CylinderCellGeomContainer_H__
#define PHG4CylinderCellGeomContainer_H__

#include <phool/PHObject.h>

#include <map>

class PHG4CylinderCellGeom;

class PHG4CylinderCellGeomContainer: public PHObject
{
 public:
  typedef std::map<int, PHG4CylinderCellGeom*> Map;
  typedef Map::iterator Iterator;
  typedef Map::const_iterator ConstIterator;
  typedef std::pair<Iterator, Iterator> Range;
  typedef std::pair<ConstIterator, ConstIterator> ConstRange;

  PHG4CylinderCellGeomContainer();
  virtual ~PHG4CylinderCellGeomContainer();

  void identify(std::ostream& os = std::cout) const;

  int AddLayerCellGeom(const int i, PHG4CylinderCellGeom *mygeom);
  int AddLayerCellGeom(PHG4CylinderCellGeom *mygeom);
  PHG4CylinderCellGeom *GetLayerCellGeom(const int i);
  PHG4CylinderCellGeom *GetFirstLayerCellGeom();
  int get_NLayers() const {return layergeoms.size();}
  std::pair<std::map<int,PHG4CylinderCellGeom *>::const_iterator, std::map<int,PHG4CylinderCellGeom *>::const_iterator> get_begin_end() const {return std::make_pair(layergeoms.begin(), layergeoms.end());}

 protected:
  std::map<int,PHG4CylinderCellGeom *> layergeoms ;
  ClassDef(PHG4CylinderCellGeomContainer,1)
};

#endif
