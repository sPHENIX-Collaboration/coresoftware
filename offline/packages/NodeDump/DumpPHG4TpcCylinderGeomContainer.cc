#include "DumpPHG4TpcCylinderGeomContainer.h"

#include <phool/PHIODataNode.h>

#include <g4detectors/PHG4TpcCylinderGeom.h>
#include <g4detectors/PHG4TpcCylinderGeomContainer.h>

#include <map>
#include <ostream>
#include <string>
#include <utility>

using MyNode_t = PHIODataNode<PHG4TpcCylinderGeomContainer>;

DumpPHG4TpcCylinderGeomContainer::DumpPHG4TpcCylinderGeomContainer(const std::string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpPHG4TpcCylinderGeomContainer::process_Node(PHNode *myNode)
{
  PHG4TpcCylinderGeomContainer *phg4geomcontainer = nullptr;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);  // NOLINT(cppcoreguidelines-pro-type-static-cast-downcast)
  if (thisNode)
  {
    phg4geomcontainer = thisNode->getData();
  }
  if (phg4geomcontainer)
  {
    PHG4TpcCylinderGeomContainer::ConstIterator hiter;
    PHG4TpcCylinderGeomContainer::ConstRange geom_begin_end = phg4geomcontainer->get_begin_end();
    *fout << "num layers: " << phg4geomcontainer->get_NLayers() << std::endl;
    for (hiter = geom_begin_end.first; hiter != geom_begin_end.second; hiter++)
    {
      *fout << *(hiter->second) << std::endl;
    }
  }
  return 0;
}
