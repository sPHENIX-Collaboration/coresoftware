#include "DumpBbcPmtContainer.h"

#include "bbc/BbcPmtContainer.h"
#include "bbc/BbcPmtHit.h"

#include <phool/PHIODataNode.h>

#include <ostream>
#include <string>

using MyNode_t = PHIODataNode<BbcPmtContainer>;

DumpBbcPmtContainer::DumpBbcPmtContainer(const std::string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpBbcPmtContainer::process_Node(PHNode *myNode)
{
  BbcPmtContainer *bbcpmtcontainer = nullptr;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    bbcpmtcontainer = thisNode->getData();
  }
  if (bbcpmtcontainer && bbcpmtcontainer->isValid())
  {
    *fout << "BbcPmtContainer->get_npmt: " << bbcpmtcontainer->get_npmt() << std::endl;
    for (int j = 0; j < bbcpmtcontainer->get_npmt(); j++)
    {
      BbcPmtHit *bbcpmt = bbcpmtcontainer->get_pmt(j);
      *fout << "BbcPmtHit->get_pmt(" << j << "): " << bbcpmt->get_pmt() << std::endl;
      *fout << "BbcPmtHit->get_q(" << j << "): " << bbcpmt->get_q() << std::endl;
      *fout << "BbcPmtHit->get_time(" << j << "): " << bbcpmt->get_time() << std::endl;
    }
  }
  return 0;
}
