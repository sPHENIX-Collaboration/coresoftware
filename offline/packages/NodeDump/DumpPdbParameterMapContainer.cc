#include "DumpPdbParameterMapContainer.h"

#include <phool/PHIODataNode.h>

#include <pdbcalbase/PdbParameterMap.h>
#include <pdbcalbase/PdbParameterMapContainer.h>

#include <map>
#include <ostream>
#include <string>
#include <utility>

using MyNode_t = PHIODataNode<PdbParameterMapContainer>;

DumpPdbParameterMapContainer::DumpPdbParameterMapContainer(const std::string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpPdbParameterMapContainer::process_Node(PHNode *myNode)
{
  PdbParameterMapContainer *pdbparams = nullptr;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    pdbparams = thisNode->getData();
  }
  if (pdbparams)
  {
    PdbParameterMapContainer::parIter piter;
    PdbParameterMapContainer::parConstRange pbegin_end = pdbparams->get_ParameterMaps();
    for (piter = pbegin_end.first; piter != pbegin_end.second; ++piter)
    {
      *fout << "PdbParameterMap # " << piter->first << std::endl;
      PdbParameterMap::dIter diter;
      PdbParameterMap::dConstRange dbegin_end = piter->second->get_dparam_iters();
      for (diter = dbegin_end.first; diter != dbegin_end.second; ++diter)
      {
        *fout << "name: " << diter->first << ": value " << diter->second << std::endl;
      }
      PdbParameterMap::iIter iiter;
      PdbParameterMap::iConstRange ibegin_end = piter->second->get_iparam_iters();
      for (iiter = ibegin_end.first; iiter != ibegin_end.second; ++iiter)
      {
        *fout << "name: " << iiter->first << ": value " << iiter->second << std::endl;
      }
      PdbParameterMap::strIter striter;
      PdbParameterMap::strConstRange strbegin_end = piter->second->get_cparam_iters();
      for (striter = strbegin_end.first; striter != strbegin_end.second; ++striter)
      {
        *fout << "name: " << striter->first << ": value " << striter->second << std::endl;
      }
    }
  }
  return 0;
}
