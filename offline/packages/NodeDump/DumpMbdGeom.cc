#include "DumpMbdGeom.h"

#include <mbd/MbdGeom.h>

#include <phool/PHIODataNode.h>

#include <ostream>
#include <string>

using MyNode_t = PHIODataNode<MbdGeom>;

DumpMbdGeom::DumpMbdGeom(const std::string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpMbdGeom::process_Node(PHNode *myNode)
{
  MbdGeom *mbdgeom = nullptr;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    mbdgeom = thisNode->getData();
  }
  if (mbdgeom)
  {
    for (int i=0; i<128; i++)
    {
      *fout << "get_x(" << i << "): " << mbdgeom->get_x(i) << std::endl;
      *fout << "get_y(" << i << "): " << mbdgeom->get_y(i) << std::endl;
      *fout << "get_z(" << i << "): " << mbdgeom->get_z(i) << std::endl;
      *fout << "get_r(" << i << "): " << mbdgeom->get_r(i) << std::endl;
      *fout << "get_phi(" << i << "): " << mbdgeom->get_phi(i) << std::endl;
      *fout << "get_arm(" << i << "): " << mbdgeom->get_arm(i) << std::endl;
      *fout << "get_arm_feech(" << i << "): " << mbdgeom->get_arm_feech(i) << std::endl;
      *fout << "get_pmt(" << i << "): " << mbdgeom->get_pmt(i) << std::endl;
      *fout << "get_type(" << i << "): " << mbdgeom->get_type(i) << std::endl;
    }
  }
  return 0;
}
