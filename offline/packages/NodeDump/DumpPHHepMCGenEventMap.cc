#include "DumpPHHepMCGenEventMap.h"

#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>

#include <HepMC/GenEvent.h>

#include <phool/PHIODataNode.h>

#include <map>
#include <ostream>
#include <string>
#include <utility>

using namespace std;

typedef PHIODataNode<PHHepMCGenEventMap> MyNode_t;

DumpPHHepMCGenEventMap::DumpPHHepMCGenEventMap(const string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpPHHepMCGenEventMap::process_Node(PHNode *myNode)
{
  PHHepMCGenEventMap *phhepmcgeneventmap = nullptr;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    phhepmcgeneventmap = thisNode->getData();
  }
  if (phhepmcgeneventmap)
  {
    PHHepMCGenEventMap::ConstIter iter_beg = phhepmcgeneventmap->begin();
    PHHepMCGenEventMap::ConstIter iter_end = phhepmcgeneventmap->end();
    *fout << "size: " << phhepmcgeneventmap->size() << endl;
    for (PHHepMCGenEventMap::ConstIter iter = iter_beg; iter != iter_end; ++iter)
    {
      *fout << "map entry: " << iter->first << endl; 
      PHHepMCGenEvent *genevt = iter->second;
      HepMC::GenEvent *evt = genevt->getEvent();
      *fout << "Embedding id " << genevt->get_embedding_id() << endl;
      *fout << "is simulated " << genevt->is_simulated() << endl;
      *fout << "Collision vertex x: " <<  genevt->get_collision_vertex().x() << endl;
      *fout << "Collision vertex y: " <<  genevt->get_collision_vertex().y() << endl;
      *fout << "Collision vertex z: " <<  genevt->get_collision_vertex().z() << endl;
      *fout << "Collision vertex t: " <<  genevt->get_collision_vertex().t() << endl;
      evt->print(*fout);
    }
  }
  return 0;
}
