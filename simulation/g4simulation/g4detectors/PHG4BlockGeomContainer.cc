#include "PHG4BlockGeomContainer.h"

#include "PHG4BlockGeom.h"

#include <cmath>

using namespace std;

PHG4BlockGeomContainer::PHG4BlockGeomContainer()
  : _magfield(NAN)
{
}

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
  os << "mag field: " << _magfield << endl;
  os << "number of layers: " << _layergeoms.size() << endl;
  map<int, PHG4BlockGeom *>::const_iterator iter;
  for (iter = _layergeoms.begin(); iter != _layergeoms.end(); ++iter)
  {
    (iter->second)->identify(os);
  }
  return;
}

int PHG4BlockGeomContainer::AddLayerGeom(const int i, PHG4BlockGeom *mygeom)
{
  if (_layergeoms.find(i) != _layergeoms.end())
  {
    cout << "layer " << i << " already added to PHBlockGeomContainer" << endl;
    return -1;
  }
  mygeom->set_layer(i);
  _layergeoms[i] = mygeom;
  return 0;
}

int PHG4BlockGeomContainer::AddLayerGeom(PHG4BlockGeom *mygeom)
{
  int layer = mygeom->get_layer();
  if (_layergeoms.find(layer) != _layergeoms.end())
  {
    cout << "layer " << layer << " already added to PHBlockGeomContainer" << endl;
    return -1;
  }
  _layergeoms[layer] = mygeom;
  return 0;
}

PHG4BlockGeom *
PHG4BlockGeomContainer::GetLayerGeom(const int i)
{
  map<int, PHG4BlockGeom *>::const_iterator iter = _layergeoms.find(i);
  if (iter != _layergeoms.end())
  {
    return iter->second;
  }
  cout << "Could not locate layer " << i << " in PHG4BlockGeomContainer" << endl;
  return nullptr;
}
