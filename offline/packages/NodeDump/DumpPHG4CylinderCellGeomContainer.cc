#include "DumpPHG4CylinderCellGeomContainer.h"

#include <phool/PHIODataNode.h>

#include <g4detectors/PHG4CylinderCellDefs.h>
#include <g4detectors/PHG4CylinderCellGeom.h>
#include <g4detectors/PHG4CylinderCellGeomContainer.h>

#include <map>
#include <ostream>
#include <string>
#include <utility>

using MyNode_t = PHIODataNode<PHG4CylinderCellGeomContainer>;

DumpPHG4CylinderCellGeomContainer::DumpPHG4CylinderCellGeomContainer(const std::string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpPHG4CylinderCellGeomContainer::process_Node(PHNode *myNode)
{
  PHG4CylinderCellGeomContainer *phg4geomcontainer = nullptr;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    phg4geomcontainer = thisNode->getData();
  }
  if (phg4geomcontainer)
  {
    PHG4CylinderCellGeomContainer::ConstIterator hiter;
    PHG4CylinderCellGeomContainer::ConstRange geom_begin_end = phg4geomcontainer->get_begin_end();
    *fout << "num layers: " << phg4geomcontainer->get_NLayers() << std::endl;
    for (hiter = geom_begin_end.first; hiter != geom_begin_end.second; hiter++)
    {
      *fout << "layer: " << hiter->second->get_layer() << std::endl;
      *fout << "radius: " << hiter->second->get_radius() << std::endl;
      *fout << "thickness: " << hiter->second->get_thickness() << std::endl;
      int binning = hiter->second->get_binning();
      *fout << "binning: " << binning << std::endl;
      switch (binning)
      {
      case PHG4CylinderCellDefs::sizebinning:

        *fout << "zbins: " << hiter->second->get_zbins() << std::endl;
        *fout << "zmin: " << hiter->second->get_zmin() << std::endl;
        *fout << "zstep: " << hiter->second->get_zstep() << std::endl;
        break;
      case PHG4CylinderCellDefs::etaphibinning:
        *fout << "etabins: " << hiter->second->get_etabins() << std::endl;
        *fout << "etastep: " << hiter->second->get_etastep() << std::endl;
        *fout << "etamin: " << hiter->second->get_etamin() << std::endl;
        break;
      default:
        break;
      }
      *fout << "phibins: " << hiter->second->get_phibins() << std::endl;
      *fout << "phistep: " << hiter->second->get_phistep() << std::endl;
      *fout << "phimin: " << hiter->second->get_phimin() << std::endl;
    }
  }
  return 0;
}
