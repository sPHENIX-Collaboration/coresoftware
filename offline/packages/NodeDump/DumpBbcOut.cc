#include "DumpBbcOut.h"


#include <bbc/BbcOut.h>

#include <phool/PHIODataNode.h>

#include <ostream>
#include <string>


typedef PHIODataNode<BbcOut> MyNode_t;

DumpBbcOut::DumpBbcOut(const std::string &NodeName): DumpObject(NodeName)
{
  return ;
}

int DumpBbcOut::process_Node(PHNode *myNode)
{
  BbcOut *bbcout = nullptr;
  MyNode_t *thisNode = static_cast <MyNode_t *> (myNode);
  if (thisNode)
    {
      bbcout = thisNode->getData();
    }
  if (bbcout && bbcout->isValid())
    {
      *fout << "BbcOut->get_VertexPoint: " << bbcout->get_VertexPoint() << std::endl;
      *fout << "BbcOut->get_dVertexPoint: " << bbcout->get_dVertexPoint() << std::endl;
      *fout << "BbcOut->get_TimeZero: " << bbcout->get_TimeZero() << std::endl;
      *fout << "BbcOut->get_dTimeZero: " << bbcout->get_dTimeZero() << std::endl;

      for (int j = 0; j < 2; j++)
      {
	*fout << "BbcOut->get_nPMT(" << j <<"): " << bbcout->get_nPMT(j) << std::endl;
	*fout << "BbcOut->get_nCharge(" << j << "): " << bbcout->get_nCharge(j) << std::endl;
	*fout << "BbcOut->get_Timing(" << j << "): " << bbcout->get_Timing(j) << std::endl;
      }
    }
  return 0;
}

