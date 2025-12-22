#include "PHG4BlockGeomContainer.h"

#include "PHG4BlockGeom.h"

PHG4BlockGeomContainer::~PHG4BlockGeomContainer()
{
  while (_layergeoms.begin() != _layergeoms.end())
  {
    delete _layergeoms.begin()->second;
    _layergeoms.erase(_layergeoms.begin());
  }
  return;
}

void PHG4BlockGeomContainer::identify(std::ostream &os) const
{
  os << "mag field: " << _magfield << std::endl;
  os << "number of layers: " << _layergeoms.size() << std::endl;
  for (auto _layergeom : _layergeoms)
  {
    (_layergeom.second)->identify(os);
  }
  return;
}

int PHG4BlockGeomContainer::AddLayerGeom(const int i, PHG4BlockGeom *mygeom)
{
  if (_layergeoms.contains(i))
  {
    std::cout << "layer " << i << " already added to PHBlockGeomContainer" << std::endl;
    return -1;
  }
  mygeom->set_layer(i);
  _layergeoms[i] = mygeom;
  return 0;
}

int PHG4BlockGeomContainer::AddLayerGeom(PHG4BlockGeom *mygeom)
{
  int layer = mygeom->get_layer();
  if (_layergeoms.contains(layer))
  {
    std::cout << "layer " << layer << " already added to PHBlockGeomContainer" << std::endl;
    return -1;
  }
  _layergeoms[layer] = mygeom;
  return 0;
}

PHG4BlockGeom *
PHG4BlockGeomContainer::GetLayerGeom(const int i)
{
  auto const iter = _layergeoms.find(i);
  if (iter != _layergeoms.end())
  {
    return iter->second;
  }
  std::cout << "Could not locate layer " << i << " in PHG4BlockGeomContainer" << std::endl;
  return nullptr;
}
