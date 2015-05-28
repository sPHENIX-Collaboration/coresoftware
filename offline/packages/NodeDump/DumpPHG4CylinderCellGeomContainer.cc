#include "DumpPHG4CylinderGeomContainer.h"

#include <phool/PHIODataNode.h>

#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <g4detectors/PHG4CylinderGeom.h>

#include <string>

using namespace std;

typedef PHIODataNode<PHG4CylinderGeomContainer> MyNode_t;

DumpPHG4CylinderGeomContainer::DumpPHG4CylinderGeomContainer(const string &NodeName): DumpObject(NodeName)
{
  return ;
}

int DumpPHG4CylinderGeomContainer::process_Node(PHNode *myNode)
{
  PHG4CylinderGeomContainer *phg4geomcontainer = NULL;
  MyNode_t *thisNode = static_cast <MyNode_t *> (myNode);
  if (thisNode)
    {
      phg4geomcontainer = thisNode->getData();
    }
  if (phg4geomcontainer)
    {
      PHG4CylinderGeomContainer::ConstIterator hiter;
      PHG4CylinderGeomContainer::ConstRange geom_begin_end = phg4geomcontainer->get_begin_end();
      *fout << "num layers: " << phg4geomcontainer->get_NLayers() << endl;
      for (hiter = geom_begin_end.first; hiter != geom_begin_end.second; hiter++)
        {
          *fout << "layer: " << hiter->second->get_layer() << endl;
          *fout << "radius: " << hiter->second->get_radius() << endl;
          *fout << "thickness: " << hiter->second->get_thickness() << endl;
          *fout << "zmin: " << hiter->second->get_zmin() << endl;
          *fout << "zmax: " << hiter->second->get_zmax() << endl;
        }
    }
  return 0;
}

