#include "PHG4CylinderCellGeomContainer.h"

#include "PHG4CylinderCellGeom.h"

PHG4CylinderCellGeomContainer::~PHG4CylinderCellGeomContainer()
{
  while (layergeoms.begin() != layergeoms.end())
  {
    delete layergeoms.begin()->second;
    layergeoms.erase(layergeoms.begin());
  }
  return;
}

void PHG4CylinderCellGeomContainer::identify(std::ostream &os) const
{
  for (auto layergeom : layergeoms)
  {
    std::cout << "layer " << layergeom.first << std::endl;
    (layergeom.second)->identify(os);
  }
  return;
}

int PHG4CylinderCellGeomContainer::AddLayerCellGeom(const int i, PHG4CylinderCellGeom *mygeom)
{
  if (layergeoms.contains(i))
  {
    std::cout << "layer " << i << " already added to PHCylinderCellGeomContainer" << std::endl;
    return -1;
  }
  mygeom->set_layer(i);
  layergeoms[i] = mygeom;
  return 0;
}

int PHG4CylinderCellGeomContainer::AddLayerCellGeom(PHG4CylinderCellGeom *mygeom)
{
  int layer = mygeom->get_layer();
  if (layergeoms.contains(layer))
  {
    std::cout << "layer " << layer << " already added to PHCylinderCellGeomContainer" << std::endl;
    return -1;
  }
  layergeoms[layer] = mygeom;
  return 0;
}

PHG4CylinderCellGeom *
PHG4CylinderCellGeomContainer::GetLayerCellGeom(const int i)
{
  const auto iter = layergeoms.find(i);
  if (iter != layergeoms.end())
  {
    return iter->second;
  }
  std::cout << "Could not locate layer " << i << " in PHG4CylinderCellGeomContainer" << std::endl;
  return nullptr;
}

PHG4CylinderCellGeom *
PHG4CylinderCellGeomContainer::GetFirstLayerCellGeom()
{
  if (layergeoms.empty())
  {
    return nullptr;
  }
  return layergeoms.begin()->second;
}
