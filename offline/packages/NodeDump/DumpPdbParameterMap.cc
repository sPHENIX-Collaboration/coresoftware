#include "DumpPdbParameterMap.h"

#include <phool/PHIODataNode.h>

#include <pdbcalbase/PdbParameterMap.h>

#include <string>

using namespace std;

typedef PHIODataNode<PdbParameterMap> MyNode_t;

DumpPdbParameterMap::DumpPdbParameterMap(const string &NodeName): DumpObject(NodeName)
{
  return ;
}

int DumpPdbParameterMap::process_Node(PHNode *myNode)
{
  PdbParameterMap *pdbparams = NULL;
  MyNode_t *thisNode = static_cast <MyNode_t *> (myNode);
  if (thisNode)
    {
      pdbparams = thisNode->getData();
    }
  if (pdbparams)
    {
      PdbParameterMap::dIter diter;
      PdbParameterMap::dConstRange dbegin_end = pdbparams->get_dparam_iters();
      for (diter=dbegin_end.first; diter != dbegin_end.second; ++diter)
        {
          *fout << "name: " << diter->first << ": value " << diter->second << endl;
        }
    }
  return 0;
}

