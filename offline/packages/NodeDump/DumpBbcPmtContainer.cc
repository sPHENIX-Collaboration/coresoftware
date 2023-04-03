#include "DumpBbcPmtContainer.h"


#include <bbc/BbcPmtContainer.h>

#include <phool/PHIODataNode.h>

#include <ostream>
#include <string>


typedef PHIODataNode<BbcPmtContainer> MyNode_t;

DumpBbcPmtContainer::DumpBbcPmtContainer(const std::string &NodeName): DumpObject(NodeName)
{
  return ;
}

int DumpBbcPmtContainer::process_Node(PHNode *myNode)
{
  BbcPmtContainer *bbcpmtcontainer = nullptr;
  MyNode_t *thisNode = static_cast <MyNode_t *> (myNode);
  if (thisNode)
    {
      bbcpmtcontainer = thisNode->getData();
    }
  if (bbcpmtcontainer && bbcpmtcontainer->isValid())
    {
      *fout << "BbcPmtContainer->get_npmt: " << bbcpmtcontainer->get_npmt() << std::endl;
      for (int j = 0; j < bbcpmtcontainer->get_npmt(); j++)
      {
	*fout << "BbcPmtContainer->get_pmt(" << j <<"): " << bbcpmtcontainer->get_pmt(j) << std::endl;
	*fout << "BbcPmtContainer->get_adc(" << j <<"): " << bbcpmtcontainer->get_adc(j) << std::endl;
	*fout << "BbcPmtContainer->get_tdc0(" << j << "): " << bbcpmtcontainer->get_tdc0(j) << std::endl;
	*fout << "BbcPmtContainer->get_tdc1(" << j << "): " << bbcpmtcontainer->get_tdc1(j) << std::endl;
      }
    }
  return 0;
}

