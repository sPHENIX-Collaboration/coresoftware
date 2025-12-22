#include "PHG4BlockCellGeomContainer.h"

#include "PHG4BlockCellGeom.h"

PHG4BlockCellGeomContainer::~PHG4BlockCellGeomContainer()
{
  while (layergeoms.begin() != layergeoms.end())
  {
    delete layergeoms.begin()->second;
    layergeoms.erase(layergeoms.begin());
  }
  return;
}

void PHG4BlockCellGeomContainer::identify(std::ostream &os) const
{
  for (auto layergeom : layergeoms)
  {
    os << "layer " << layergeom.first << std::endl;
    (layergeom.second)->identify(os);
  }
  return;
}

int PHG4BlockCellGeomContainer::AddLayerCellGeom(const int i, PHG4BlockCellGeom *mygeom)
{
  if (layergeoms.contains(i))
  {
    std::cout << "layer " << i << " already added to PHBlockCellGeomContainer" << std::endl;
    return -1;
  }
  mygeom->set_layer(i);
  layergeoms[i] = mygeom;
  return 0;
}

int PHG4BlockCellGeomContainer::AddLayerCellGeom(PHG4BlockCellGeom *mygeom)
{
  int layer = mygeom->get_layer();
  if (layergeoms.contains(layer))
  {
    std::cout << "layer " << layer << " already added to PHBlockCellGeomContainer" << std::endl;
    return -1;
  }
  layergeoms[layer] = mygeom;
  return 0;
}

PHG4BlockCellGeom *
PHG4BlockCellGeomContainer::GetLayerCellGeom(const int i)
{
  const auto iter = layergeoms.find(i);
  if (iter != layergeoms.end())
  {
    return iter->second;
  }
  std::cout << "Could not locate layer " << i << " in PHG4BlockCellGeomContainer" << std::endl;
  return nullptr;
}
