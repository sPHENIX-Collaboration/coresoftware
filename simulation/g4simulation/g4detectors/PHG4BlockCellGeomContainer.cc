#include "PHG4BlockCellGeomContainer.h"

#include "PHG4BlockCellGeom.h"

using namespace std;

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
  map<int, PHG4BlockCellGeom *>::const_iterator iter;
  for (iter = layergeoms.begin(); iter != layergeoms.end(); ++iter)
  {
    cout << "layer " << iter->first << endl;
    (iter->second)->identify(os);
  }
  return;
}

int PHG4BlockCellGeomContainer::AddLayerCellGeom(const int i, PHG4BlockCellGeom *mygeom)
{
  if (layergeoms.find(i) != layergeoms.end())
  {
    cout << "layer " << i << " already added to PHBlockCellGeomContainer" << endl;
    return -1;
  }
  mygeom->set_layer(i);
  layergeoms[i] = mygeom;
  return 0;
}

int PHG4BlockCellGeomContainer::AddLayerCellGeom(PHG4BlockCellGeom *mygeom)
{
  int layer = mygeom->get_layer();
  if (layergeoms.find(layer) != layergeoms.end())
  {
    cout << "layer " << layer << " already added to PHBlockCellGeomContainer" << endl;
    return -1;
  }
  layergeoms[layer] = mygeom;
  return 0;
}

PHG4BlockCellGeom *
PHG4BlockCellGeomContainer::GetLayerCellGeom(const int i)
{
  map<int, PHG4BlockCellGeom *>::const_iterator iter = layergeoms.find(i);
  if (iter != layergeoms.end())
  {
    return iter->second;
  }
  cout << "Could not locate layer " << i << " in PHG4BlockCellGeomContainer" << endl;
  return nullptr;
}
