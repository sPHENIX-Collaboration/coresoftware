#include "PHG4CylinderGeomContainer.h"
#include "PHG4CylinderGeom.h"

#include <cmath>

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
  os << "mag field: " << magfield << std::endl;
  os << "number of layers: " << layergeoms.size() << std::endl;
  std::map<int, PHG4CylinderGeom *>::const_iterator iter;
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
    std::cout << "layer " << i << " already added to PHCylinderGeomContainer" << std::endl;
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
    std::cout << "layer " << layer << " already added to PHCylinderGeomContainer" << std::endl;
    return -1;
  }
  layergeoms[layer] = mygeom;
  return 0;
}

PHG4CylinderGeom *
PHG4CylinderGeomContainer::GetLayerGeom(const int i)
{
  std::map<int, PHG4CylinderGeom *>::const_iterator iter = layergeoms.find(i);
  if (iter != layergeoms.end())
  {
    return iter->second;
  }
  std::cout << "Could not locate layer " << i << " in PHG4CylinderGeomContainer" << std::endl;
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
