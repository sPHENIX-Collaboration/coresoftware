#include "PHG4TpcCylinderGeomContainer.h"

#include "PHG4TpcCylinderGeom.h"

using namespace std;

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
  map<int, PHG4TpcCylinderGeom *>::const_iterator iter;
  for (iter = layergeoms.begin(); iter != layergeoms.end(); ++iter)
  {
    cout << "layer " << iter->first << endl;
    (iter->second)->identify(os);
  }
  return;
}

int PHG4TpcCylinderGeomContainer::AddLayerCellGeom(const int i, PHG4TpcCylinderGeom *mygeom)
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

int PHG4TpcCylinderGeomContainer::AddLayerCellGeom(PHG4TpcCylinderGeom *mygeom)
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

PHG4TpcCylinderGeom *
PHG4TpcCylinderGeomContainer::GetLayerCellGeom(const int i)
{
  map<int, PHG4TpcCylinderGeom *>::const_iterator iter = layergeoms.find(i);
  if (iter != layergeoms.end())
  {
    return iter->second;
  }
  cout << "Could not locate layer " << i << " in PHG4TpcCylinderGeomContainer" << endl;
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
