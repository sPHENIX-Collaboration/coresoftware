#include "DumpPHG4BlockGeomContainer.h"

#include <phool/PHIODataNode.h>

#include <g4detectors/PHG4BlockGeom.h>
#include <g4detectors/PHG4BlockGeomContainer.h>

#include <string>

using namespace std;

typedef PHIODataNode<PHG4BlockGeomContainer> MyNode_t;

DumpPHG4BlockGeomContainer::DumpPHG4BlockGeomContainer(const string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpPHG4BlockGeomContainer::process_Node(PHNode *myNode)
{
  PHG4BlockGeomContainer *phg4geomcontainer = NULL;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    phg4geomcontainer = thisNode->getData();
  }
  if (phg4geomcontainer)
  {
    PHG4BlockGeomContainer::ConstIterator hiter;
    PHG4BlockGeomContainer::ConstRange geom_begin_end = phg4geomcontainer->get_begin_end();
    *fout << "num layers: " << phg4geomcontainer->get_NLayers() << endl;
    for (hiter = geom_begin_end.first; hiter != geom_begin_end.second; hiter++)
    {
      *fout << "layer: " << hiter->second->get_layer() << endl;
      *fout << "size_x: " << hiter->second->get_size_x() << endl;
      *fout << "size_y: " << hiter->second->get_size_y() << endl;
      *fout << "size_z: " << hiter->second->get_size_z() << endl;
      *fout << "center_x: " << hiter->second->get_center_x() << endl;
      *fout << "center_y: " << hiter->second->get_center_y() << endl;
      *fout << "center_z: " << hiter->second->get_center_z() << endl;
      *fout << "z_rot: " << hiter->second->get_z_rot() << endl;
      *fout << "width: " << hiter->second->get_width() << endl;
      *fout << "length: " << hiter->second->get_length() << endl;
      *fout << "thickness: " << hiter->second->get_thickness() << endl;
    }
  }
  return 0;
}
