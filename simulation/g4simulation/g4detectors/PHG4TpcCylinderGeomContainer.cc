#include "PHG4TpcCylinderGeomContainer.h"

#include "PHG4TpcCylinderGeom.h"

PHG4TpcCylinderGeomContainer::~PHG4TpcCylinderGeomContainer()
{
  while (layergeoms.begin() != layergeoms.end())
  {
    delete layergeoms.begin()->second;
    layergeoms.erase(layergeoms.begin());
  }
  return;
}

void PHG4TpcCylinderGeomContainer::identify(std::ostream &os) const
{
  std::map<int, PHG4TpcCylinderGeom *>::const_iterator iter;
  for (iter = layergeoms.begin(); iter != layergeoms.end(); ++iter)
  {
    std::cout << "layer " << iter->first << std::endl;
    (iter->second)->identify(os);
  }
  return;
}

int PHG4TpcCylinderGeomContainer::AddLayerCellGeom(const int i, PHG4TpcCylinderGeom *mygeom)
{
  if (layergeoms.find(i) != layergeoms.end())
  {
    std::cout << "layer " << i << " already added to PHCylinderCellGeomContainer" << std::endl;
    return -1;
  }
  mygeom->set_layer(i);
  layergeoms[i] = mygeom;
  return 0;
}

int PHG4TpcCylinderGeomContainer::AddLayerCellGeom(PHG4TpcCylinderGeom *mygeom)
{
  int layer = mygeom->get_layer();
  if (layergeoms.find(layer) != layergeoms.end())
  {
    std::cout << "layer " << layer << " already added to PHCylinderCellGeomContainer" << std::endl;
    return -1;
  }
  layergeoms[layer] = mygeom;
  return 0;
}

PHG4TpcCylinderGeom *
PHG4TpcCylinderGeomContainer::GetLayerCellGeom(const int i)
{
  std::map<int, PHG4TpcCylinderGeom *>::const_iterator iter = layergeoms.find(i);
  if (iter != layergeoms.end())
  {
    return iter->second;
  }
  std::cout << "Could not locate layer " << i << " in PHG4TpcCylinderGeomContainer" << std::endl;
  return nullptr;
}

PHG4TpcCylinderGeom *
PHG4TpcCylinderGeomContainer::GetFirstLayerCellGeom()
{
  if (layergeoms.empty())
  {
    return nullptr;
  }
  return layergeoms.begin()->second;
}
