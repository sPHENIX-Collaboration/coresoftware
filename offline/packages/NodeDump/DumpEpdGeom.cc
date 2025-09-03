#include "DumpEpdGeom.h"

#include <phool/PHIODataNode.h>

#include <calobase/TowerInfoDefs.h>
#include <epd/EpdGeom.h>

#include <iomanip>
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
  const auto original_precision = (*fout).precision();
  EpdGeom *epdgeom = nullptr;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);  // NOLINT(cppcoreguidelines-pro-type-static-cast-downcast)
  if (thisNode)
  {
    epdgeom = thisNode->getData();
  }
  if (epdgeom)
  {
    (*fout).precision(std::numeric_limits<float>::max_digits10);
    for (int iarm = 0; iarm < 2; iarm++)
    {
      for (int irad = 0; irad < 16; irad++)
      {
        for (int iphi = 0; iphi < 24; iphi++)
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
      (*fout).precision(original_precision);
    }
  }
  return 0;
}
