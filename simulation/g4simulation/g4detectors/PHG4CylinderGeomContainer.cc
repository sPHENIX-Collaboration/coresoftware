#include "PHG4CylinderGeomContainer.h"
#include <cmath>
#include "PHG4CylinderGeom.h"

using namespace std;

PHG4CylinderGeomContainer::PHG4CylinderGeomContainer()
  : magfield(NAN)
{
}

PHG4CylinderGeomContainer::~PHG4CylinderGeomContainer()
{
  while (layergeoms.begin() != layergeoms.end())
  {
    delete layergeoms.begin()->second;
    layergeoms.erase(layergeoms.begin());
  }
  return;
}

void PHG4CylinderGeomContainer::identify(std::ostream &os) const
{
  os << "mag field: " << magfield << endl;
  os << "number of layers: " << layergeoms.size() << endl;
  map<int, PHG4CylinderGeom *>::const_iterator iter;
  for (iter = layergeoms.begin(); iter != layergeoms.end(); ++iter)
  {
    (iter->second)->identify(os);
  }

  return;
}

int PHG4CylinderGeomContainer::AddLayerGeom(const int i, PHG4CylinderGeom *mygeom)
{
  if (layergeoms.find(i) != layergeoms.end())
  {
    cout << "layer " << i << " already added to PHCylinderGeomContainer" << endl;
    return -1;
  }
  mygeom->set_layer(i);
  layergeoms[i] = mygeom;
  return 0;
}

int PHG4CylinderGeomContainer::AddLayerGeom(PHG4CylinderGeom *mygeom)
{
  int layer = mygeom->get_layer();
  if (layergeoms.find(layer) != layergeoms.end())
  {
    cout << "layer " << layer << " already added to PHCylinderGeomContainer" << endl;
    return -1;
  }
  layergeoms[layer] = mygeom;
  return 0;
}

PHG4CylinderGeom *
PHG4CylinderGeomContainer::GetLayerGeom(const int i)
{
  map<int, PHG4CylinderGeom *>::const_iterator iter = layergeoms.find(i);
  if (iter != layergeoms.end())
  {
    return iter->second;
  }
  cout << "Could not locate layer " << i << " in PHG4CylinderGeomContainer" << endl;
  return nullptr;
}

PHG4CylinderGeom *
PHG4CylinderGeomContainer::GetFirstLayerGeom()
{
  if (layergeoms.empty())
  {
    return nullptr;
  }
  return layergeoms.begin()->second;
}
