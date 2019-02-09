#include "DumpPHG4CylinderCellGeomContainer.h"

#include <phool/PHIODataNode.h>

#include <g4detectors/PHG4CylinderCellDefs.h>
#include <g4detectors/PHG4CylinderCellGeom.h>
#include <g4detectors/PHG4CylinderCellGeomContainer.h>

#include <string>

using namespace std;

typedef PHIODataNode<PHG4CylinderCellGeomContainer> MyNode_t;

DumpPHG4CylinderCellGeomContainer::DumpPHG4CylinderCellGeomContainer(const string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpPHG4CylinderCellGeomContainer::process_Node(PHNode *myNode)
{
  PHG4CylinderCellGeomContainer *phg4geomcontainer = NULL;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    phg4geomcontainer = thisNode->getData();
  }
  if (phg4geomcontainer)
  {
    PHG4CylinderCellGeomContainer::ConstIterator hiter;
    PHG4CylinderCellGeomContainer::ConstRange geom_begin_end = phg4geomcontainer->get_begin_end();
    *fout << "num layers: " << phg4geomcontainer->get_NLayers() << endl;
    for (hiter = geom_begin_end.first; hiter != geom_begin_end.second; hiter++)
    {
      *fout << "layer: " << hiter->second->get_layer() << endl;
      *fout << "radius: " << hiter->second->get_radius() << endl;
      *fout << "thickness: " << hiter->second->get_thickness() << endl;
      int binning = hiter->second->get_binning();
      *fout << "binning: " << binning << endl;
      switch (binning)
      {
      case PHG4CylinderCellDefs::sizebinning:

        *fout << "zbins: " << hiter->second->get_zbins() << endl;
        *fout << "zmin: " << hiter->second->get_zmin() << endl;
        *fout << "zstep: " << hiter->second->get_zstep() << endl;
        break;
      case PHG4CylinderCellDefs::etaphibinning:
        *fout << "etabins: " << hiter->second->get_etabins() << endl;
        *fout << "etastep: " << hiter->second->get_etastep() << endl;
        *fout << "etamin: " << hiter->second->get_etamin() << endl;
        break;
      default:
        break;
      }
      *fout << "phibins: " << hiter->second->get_phibins() << endl;
      *fout << "phistep: " << hiter->second->get_phistep() << endl;
      *fout << "phimin: " << hiter->second->get_phimin() << endl;
    }
  }
  return 0;
}
