#include "DumpPHG4ScintillatorSlatContainer.h"

#include <phool/PHIODataNode.h>

#include <g4detectors/PHG4ScintillatorSlat.h>
#include <g4detectors/PHG4ScintillatorSlatContainer.h>

#include <map>
#include <ostream>
#include <string>
#include <utility>

using MyNode_t = PHIODataNode<PHG4ScintillatorSlatContainer>;

DumpPHG4ScintillatorSlatContainer::DumpPHG4ScintillatorSlatContainer(const std::string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpPHG4ScintillatorSlatContainer::process_Node(PHNode *myNode)
{
  PHG4ScintillatorSlatContainer *scinticontainer = nullptr;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    scinticontainer = thisNode->getData();
  }
  if (scinticontainer)
  {
    PHG4ScintillatorSlatContainer::ConstIterator hiter;
    PHG4ScintillatorSlatContainer::ConstRange begin_end = scinticontainer->getScintillatorSlats();
    *fout << "size: " << scinticontainer->size() << std::endl;
    for (hiter = begin_end.first; hiter != begin_end.second; ++hiter)
    {
      *fout << "get_key(): 0x" << std::hex << hiter->second->get_key() << std::dec << std::endl;
      *fout << "get_column(): " << hiter->second->get_column() << std::endl;
      *fout << "get_row(): " << hiter->second->get_row() << std::endl;
      *fout << "get_edep(): " << hiter->second->get_edep() << std::endl;
      *fout << "get_eion(): " << hiter->second->get_eion() << std::endl;
      *fout << "get_light_yield(): " << hiter->second->get_light_yield() << std::endl;
    }
  }
  return 0;
}
