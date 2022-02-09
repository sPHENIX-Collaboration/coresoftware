#include "PHG4CylinderCellGeomContainer.h"

#include "PHG4CylinderCellGeom.h"

using namespace std;

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
  map<int, PHG4CylinderCellGeom *>::const_iterator iter;
  for (iter = layergeoms.begin(); iter != layergeoms.end(); ++iter)
  {
    cout << "layer " << iter->first << endl;
    (iter->second)->identify(os);
  }
  return;
}

int PHG4CylinderCellGeomContainer::AddLayerCellGeom(const int i, PHG4CylinderCellGeom *mygeom)
{
  if (layergeoms.find(i) != layergeoms.end())
  {
    cout << "layer " << i << " already added to PHCylinderCellGeomContainer" << endl;
    return -1;
  }
  mygeom->set_layer(i);
  layergeoms[i] = mygeom;
  return 0;
}

int PHG4CylinderCellGeomContainer::AddLayerCellGeom(PHG4CylinderCellGeom *mygeom)
{
  int layer = mygeom->get_layer();
  if (layergeoms.find(layer) != layergeoms.end())
  {
    cout << "layer " << layer << " already added to PHCylinderCellGeomContainer" << endl;
    return -1;
  }
  layergeoms[layer] = mygeom;
  return 0;
}

PHG4CylinderCellGeom *
PHG4CylinderCellGeomContainer::GetLayerCellGeom(const int i)
{
  map<int, PHG4CylinderCellGeom *>::const_iterator iter = layergeoms.find(i);
  if (iter != layergeoms.end())
  {
    return iter->second;
  }
  cout << "Could not locate layer " << i << " in PHG4CylinderCellGeomContainer" << endl;
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
