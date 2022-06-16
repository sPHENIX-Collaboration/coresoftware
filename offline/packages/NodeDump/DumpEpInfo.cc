#include "DumpEpInfo.h"

#include <eventplane/EpInfo.h>

#include <phool/PHIODataNode.h>

#include <ostream>
#include <string>

using namespace std;

typedef PHIODataNode<EpInfo> MyNode_t;

DumpEpInfo::DumpEpInfo(const string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpEpInfo::process_Node(PHNode *myNode)
{
  EpInfo *epinfo = nullptr;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    epinfo = thisNode->getData();
  }
  if (epinfo)
  {
    // *fout << "EpInfo->isValid(): " << epinfo->isValid() << endl;
    // if (epinfo->isValid())
    *fout << "Max Order: " <<  epinfo->MaxOrder() << endl;
    for (unsigned int i = 0; i < epinfo->MaxOrder(); ++i)
    {
      std::pair<double, double> tmp = epinfo->ARawQ(i);
      *fout << "RawQ( " << i << "): " << tmp.first << ", " << tmp.second << endl;
      *fout << "SWRaw( " << i << "): " << epinfo->SWRaw(i) << std::endl;
      *fout << "RawPsi( " << i << "): " << epinfo->RawPsi(i) << std::endl;

    }
  }
  return 0;
}
