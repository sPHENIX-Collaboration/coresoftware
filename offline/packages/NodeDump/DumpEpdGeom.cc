#include "DumpEpdGeom.h"

#include <phool/PHIODataNode.h>

#include <calobase/TowerInfoDefs.h>
#include <epd/EpdGeom.h>

#include <ostream>
#include <string>
#include <utility>

using MyNode_t = PHIODataNode<EpdGeom>;

DumpEpdGeom::DumpEpdGeom(const std::string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpEpdGeom::process_Node(PHNode *myNode)
{
  EpdGeom *epdgeom = nullptr;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    epdgeom = thisNode->getData();
  }
  if (epdgeom)
  {
    for (int iarm = 0; iarm<2; iarm++)
    {
      for (int irad = 0; irad < 16; irad++)
      {
	for (int iphi = 0; iphi<24; iphi++)
	{ 
	  if (irad == 0 && iphi > 11)
	  {
	    continue;
	  }
	  unsigned int key = TowerInfoDefs::encode_epd(iarm, irad, iphi);
	  *fout << "tile key: 0x" << std::hex << key << std::dec << std::endl;
	  *fout << "get_r: " << epdgeom->get_r(key) << std::endl;
	  *fout << "get_phi: " << epdgeom->get_phi(key) << std::endl;
	  *fout << "get_z: " << epdgeom->get_z(key) << std::endl;
        }
      }
    }
  }
  return 0;
}
