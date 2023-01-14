#include "DumpPHG4BlockGeomContainer.h"

#include <phool/PHIODataNode.h>

#include <g4detectors/PHG4BlockGeom.h>
#include <g4detectors/PHG4BlockGeomContainer.h>

#include <map>
#include <ostream>
#include <string>
#include <utility>

using MyNode_t = PHIODataNode<PHG4BlockGeomContainer>;

DumpPHG4BlockGeomContainer::DumpPHG4BlockGeomContainer(const std::string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpPHG4BlockGeomContainer::process_Node(PHNode *myNode)
{
  PHG4BlockGeomContainer *phg4geomcontainer = nullptr;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    phg4geomcontainer = thisNode->getData();
  }
  if (phg4geomcontainer)
  {
    PHG4BlockGeomContainer::ConstIterator hiter;
    PHG4BlockGeomContainer::ConstRange geom_begin_end = phg4geomcontainer->get_begin_end();
    *fout << "num layers: " << phg4geomcontainer->get_NLayers() << std::endl;
    for (hiter = geom_begin_end.first; hiter != geom_begin_end.second; hiter++)
    {
      *fout << "layer: " << hiter->second->get_layer() << std::endl;
      *fout << "size_x: " << hiter->second->get_size_x() << std::endl;
      *fout << "size_y: " << hiter->second->get_size_y() << std::endl;
      *fout << "size_z: " << hiter->second->get_size_z() << std::endl;
      *fout << "center_x: " << hiter->second->get_center_x() << std::endl;
      *fout << "center_y: " << hiter->second->get_center_y() << std::endl;
      *fout << "center_z: " << hiter->second->get_center_z() << std::endl;
      *fout << "z_rot: " << hiter->second->get_z_rot() << std::endl;
      *fout << "width: " << hiter->second->get_width() << std::endl;
      *fout << "length: " << hiter->second->get_length() << std::endl;
      *fout << "thickness: " << hiter->second->get_thickness() << std::endl;
    }
  }
  return 0;
}
