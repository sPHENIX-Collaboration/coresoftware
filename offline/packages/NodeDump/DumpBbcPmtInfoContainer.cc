#include "DumpBbcPmtInfoContainer.h"

#include <bbc/BbcPmtInfoContainerV1.h>

#include <phool/PHIODataNode.h>

#include <ostream>
#include <string>

using MyNode_t = PHIODataNode<BbcPmtInfoContainerV1>;

DumpBbcPmtInfoContainer::DumpBbcPmtInfoContainer(const std::string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpBbcPmtInfoContainer::process_Node(PHNode *myNode)
{
  BbcPmtInfoContainerV1 *bbcpmtinfocontainer = nullptr;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    bbcpmtinfocontainer = thisNode->getData();
  }
  if (bbcpmtinfocontainer && bbcpmtinfocontainer->isValid())
  {
    *fout << "BbcPmtInfoV1Container->size: " << bbcpmtinfocontainer->size() << std::endl;
    for (size_t j = 0; j < bbcpmtinfocontainer->size(); j++)
    {
      BbcPmtInfoV1 *bbcpmt = bbcpmtinfocontainer->get_pmt(j);
      *fout << "BbcPmtInfoContainerV1->get_pmt(" << j << ")->get_pmt: " << bbcpmt->get_pmt() << std::endl;
      *fout << "BbcPmtInfoContainerV1->get_t(" << j << ")->get_t: " << bbcpmt->get_t() << std::endl;
      *fout << "BbcPmtInfoContainerV1->get_q(" << j << ")->get_q: " << bbcpmt->get_q() << std::endl;
    }
  }
  return 0;
}
