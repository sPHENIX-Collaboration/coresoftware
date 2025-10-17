#include "PHG4TpcGeomContainer.h"

#include "PHG4TpcGeom.h"

PHG4TpcGeomContainer::~PHG4TpcGeomContainer()
{
  while (layergeoms.begin() != layergeoms.end())
  {
    delete layergeoms.begin()->second;
    layergeoms.erase(layergeoms.begin());
  }
  return;
}

void PHG4TpcGeomContainer::identify(std::ostream &os) const
{
  std::map<int, PHG4TpcGeom *>::const_iterator iter;
  for (iter = layergeoms.begin(); iter != layergeoms.end(); ++iter)
  {
    std::cout << "layer " << iter->first << std::endl;
    (iter->second)->identify(os);
  }
  return;
}

int PHG4TpcGeomContainer::AddLayerCellGeom(const int i, PHG4TpcGeom *mygeom)
{
  if (layergeoms.find(i) != layergeoms.end())
  {
    std::cout << "layer " << i << " already added to PHG4TpcGeomContainer" << std::endl;
    return -1;
  }
  mygeom->set_layer(i);
  layergeoms[i] = mygeom;
  return 0;
}

int PHG4TpcGeomContainer::AddLayerCellGeom(PHG4TpcGeom *mygeom)
{
  int layer = mygeom->get_layer();
  if (layergeoms.find(layer) != layergeoms.end())
  {
    std::cout << "layer " << layer << " already added to PHG4TpcGeomContainer" << std::endl;
    return -1;
  }
  layergeoms[layer] = mygeom;
  return 0;
}

PHG4TpcGeom *
PHG4TpcGeomContainer::GetLayerCellGeom(const int i)
{
  std::map<int, PHG4TpcGeom *>::const_iterator iter = layergeoms.find(i);
  if (iter != layergeoms.end())
  {
    return iter->second;
  }
  std::cout << "Could not locate layer " << i << " in PHG4TpcGeomContainer" << std::endl;
  return nullptr;
}

PHG4TpcGeom *
PHG4TpcGeomContainer::GetFirstLayerCellGeom()
{
  if (layergeoms.empty())
  {
    return nullptr;
  }
  return layergeoms.begin()->second;
}
