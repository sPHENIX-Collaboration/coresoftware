#include "DumpPHG4BlockCellGeomContainer.h"

#include <phool/PHIODataNode.h>

#include <g4detectors/PHG4BlockCellGeom.h>
#include <g4detectors/PHG4BlockCellGeomContainer.h>
#include <g4detectors/PHG4CylinderCellDefs.h>

#include <map>
#include <ostream>
#include <string>
#include <utility>

using MyNode_t = PHIODataNode<PHG4BlockCellGeomContainer>;

DumpPHG4BlockCellGeomContainer::DumpPHG4BlockCellGeomContainer(const std::string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpPHG4BlockCellGeomContainer::process_Node(PHNode *myNode)
{
  PHG4BlockCellGeomContainer *phg4geomcontainer = nullptr;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    phg4geomcontainer = thisNode->getData();
  }
  if (phg4geomcontainer)
  {
    PHG4BlockCellGeomContainer::ConstIterator hiter;
    PHG4BlockCellGeomContainer::ConstRange geom_begin_end = phg4geomcontainer->get_begin_end();
    *fout << "num layers: " << phg4geomcontainer->get_NLayers() << std::endl;
    for (hiter = geom_begin_end.first; hiter != geom_begin_end.second; hiter++)
    {
      *fout << "layer: " << hiter->second->get_layer() << std::endl;
      int binning = hiter->second->get_binning();
      *fout << "binning: " << binning << std::endl;
      switch (binning)
      {
      case PHG4CylinderCellDefs::sizebinning:

        *fout << "zbins: " << hiter->second->get_zbins() << std::endl;
        *fout << "zmin: " << hiter->second->get_zmin() << std::endl;
        *fout << "zstep: " << hiter->second->get_zstep() << std::endl;
        *fout << "xbins: " << hiter->second->get_xbins() << std::endl;
        *fout << "xmin: " << hiter->second->get_xmin() << std::endl;
        *fout << "xstep: " << hiter->second->get_xstep() << std::endl;
        break;
      case PHG4CylinderCellDefs::etaphibinning:
      case PHG4CylinderCellDefs::etaslatbinning:
        *fout << "etabins: " << hiter->second->get_etabins() << std::endl;
        *fout << "etastep: " << hiter->second->get_etastep() << std::endl;
        *fout << "etamin: " << hiter->second->get_etamin() << std::endl;
        *fout << "xbins: " << hiter->second->get_xbins() << std::endl;
        *fout << "xmin: " << hiter->second->get_xmin() << std::endl;
        *fout << "xstep: " << hiter->second->get_xstep() << std::endl;
        break;
      default:
        break;
      }
    }
  }
  return 0;
}
