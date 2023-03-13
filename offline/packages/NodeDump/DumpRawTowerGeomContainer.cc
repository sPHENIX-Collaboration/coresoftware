#include "DumpRawTowerGeomContainer.h"

#include <phool/PHIODataNode.h>

#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>

#include <map>
#include <ostream>
#include <string>
#include <utility>

using MyNode_t = PHIODataNode<RawTowerGeomContainer>;

DumpRawTowerGeomContainer::DumpRawTowerGeomContainer(const std::string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpRawTowerGeomContainer::process_Node(PHNode *myNode)
{
  RawTowerGeomContainer *rawtowergeom = nullptr;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    rawtowergeom = thisNode->getData();
  }
  if (rawtowergeom)
  {
    *fout << "Calorimeter ID: " << rawtowergeom->get_calorimeter_id() << std::endl;
    *fout << "size: " << rawtowergeom->size() << std::endl;
    rawtowergeom->identify(*fout);
    RawTowerGeomContainer::ConstRange all_towers = rawtowergeom->get_tower_geometries();
    for (RawTowerGeomContainer::ConstIterator it = all_towers.first;
         it != all_towers.second; ++it)
    {
      it->second->identify(*fout);
    }
  }
  return 0;
}
