#include "DumpPHG4ScintillatorSlatContainer.h"

#include <phool/PHIODataNode.h>

#include <g4detectors/PHG4ScintillatorSlat.h>
#include <g4detectors/PHG4ScintillatorSlatContainer.h>

#include <string>

using namespace std;

typedef PHIODataNode<PHG4ScintillatorSlatContainer> MyNode_t;

DumpPHG4ScintillatorSlatContainer::DumpPHG4ScintillatorSlatContainer(const string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpPHG4ScintillatorSlatContainer::process_Node(PHNode *myNode)
{
  PHG4ScintillatorSlatContainer *scinticontainer = NULL;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    scinticontainer = thisNode->getData();
  }
  if (scinticontainer)
  {
    PHG4ScintillatorSlatContainer::ConstIterator hiter;
    PHG4ScintillatorSlatContainer::ConstRange begin_end = scinticontainer->getScintillatorSlats();
    *fout << "size: " << scinticontainer->size() << endl;
    for (hiter = begin_end.first; hiter != begin_end.second; ++hiter)
    {
      *fout << "get_key(): 0x" << hex << hiter->second->get_key() << dec << endl;
      *fout << "get_column(): " << hiter->second->get_column() << endl;
      *fout << "get_row(): " << hiter->second->get_row() << endl;
      *fout << "get_edep(): " << hiter->second->get_edep() << endl;
      *fout << "get_eion(): " << hiter->second->get_eion() << endl;
      *fout << "get_light_yield(): " << hiter->second->get_light_yield() << endl;
    }
  }
  return 0;
}
