#include "DumpPHG4CylinderGeomContainer.h"

#include <phool/PHIODataNode.h>

#include <g4detectors/PHG4CylinderGeom.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>

#include <map>
#include <ostream>
#include <string>
#include <utility>

using namespace std;

typedef PHIODataNode<PHG4CylinderGeomContainer> MyNode_t;

DumpPHG4CylinderGeomContainer::DumpPHG4CylinderGeomContainer(const string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpPHG4CylinderGeomContainer::process_Node(PHNode *myNode)
{
  PHG4CylinderGeomContainer *phg4geomcontainer = nullptr;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
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
      *fout << "nscint: " << hiter->second->get_nscint() << endl;
      *fout << "tiltangle: " << hiter->second->get_tiltangle() << endl;
      *fout << "strip_y_spacing: " << hiter->second->get_strip_y_spacing() << endl;
      *fout << "strip_z_spacing: " << hiter->second->get_strip_z_spacing() << endl;
      *fout << "strip_tilt: " << hiter->second->get_strip_tilt() << endl;
      *fout << "N_strip_columns: " << hiter->second->get_N_strip_columns() << endl;
      *fout << "N_strips_per_column: " << hiter->second->get_N_strips_per_column() << endl;
      *fout << "N_sensors_in_layer: " << hiter->second->get_N_sensors_in_layer() << endl;
      *fout << "pixel_z: " << hiter->second->get_pixel_z() << endl;
      *fout << "pixel_x: " << hiter->second->get_pixel_x() << endl;
      *fout << "pixel_thickness: " << hiter->second->get_pixel_thickness() << endl;
    }
  }
  return 0;
}
