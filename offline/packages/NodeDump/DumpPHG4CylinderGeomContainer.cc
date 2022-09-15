#include "DumpPHG4CylinderGeomContainer.h"

#include <phool/PHIODataNode.h>

#include <g4detectors/PHG4CylinderGeom.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>

#include <map>
#include <ostream>
#include <string>
#include <utility>

using MyNode_t = PHIODataNode<PHG4CylinderGeomContainer>;

DumpPHG4CylinderGeomContainer::DumpPHG4CylinderGeomContainer(const std::string &NodeName)
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
    *fout << "num layers: " << phg4geomcontainer->get_NLayers() << std::endl;
    for (hiter = geom_begin_end.first; hiter != geom_begin_end.second; hiter++)
    {
      *fout << "layer: " << hiter->second->get_layer() << std::endl;
      *fout << "radius: " << hiter->second->get_radius() << std::endl;
      *fout << "thickness: " << hiter->second->get_thickness() << std::endl;
      *fout << "zmin: " << hiter->second->get_zmin() << std::endl;
      *fout << "zmax: " << hiter->second->get_zmax() << std::endl;
      *fout << "nscint: " << hiter->second->get_nscint() << std::endl;
      *fout << "tiltangle: " << hiter->second->get_tiltangle() << std::endl;
      *fout << "strip_y_spacing: " << hiter->second->get_strip_y_spacing() << std::endl;
      *fout << "strip_z_spacing: " << hiter->second->get_strip_z_spacing() << std::endl;
      *fout << "strip_tilt: " << hiter->second->get_strip_tilt() << std::endl;
      *fout << "N_strip_columns: " << hiter->second->get_N_strip_columns() << std::endl;
      *fout << "N_strips_per_column: " << hiter->second->get_N_strips_per_column() << std::endl;
      *fout << "N_sensors_in_layer: " << hiter->second->get_N_sensors_in_layer() << std::endl;
      *fout << "pixel_z: " << hiter->second->get_pixel_z() << std::endl;
      *fout << "pixel_x: " << hiter->second->get_pixel_x() << std::endl;
      *fout << "pixel_thickness: " << hiter->second->get_pixel_thickness() << std::endl;
    }
  }
  return 0;
}
