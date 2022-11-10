#include "DumpPHHepMCGenEventMap.h"

#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <HepMC/GenEvent.h>
#pragma GCC diagnostic pop

#include <phool/PHIODataNode.h>

#include <map>
#include <ostream>
#include <string>
#include <utility>

using MyNode_t = PHIODataNode<PHHepMCGenEventMap>;

DumpPHHepMCGenEventMap::DumpPHHepMCGenEventMap(const std::string &NodeName)
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
    *fout << "size: " << phhepmcgeneventmap->size() << std::endl;
    for (PHHepMCGenEventMap::ConstIter iter = iter_beg; iter != iter_end; ++iter)
    {
      *fout << "map entry: " << iter->first << std::endl;
      PHHepMCGenEvent *genevt = iter->second;
      HepMC::GenEvent *evt = genevt->getEvent();
      *fout << "Embedding id " << genevt->get_embedding_id() << std::endl;
      *fout << "is simulated " << genevt->is_simulated() << std::endl;
      *fout << "Collision vertex x: " << genevt->get_collision_vertex().x() << std::endl;
      *fout << "Collision vertex y: " << genevt->get_collision_vertex().y() << std::endl;
      *fout << "Collision vertex z: " << genevt->get_collision_vertex().z() << std::endl;
      *fout << "Collision vertex t: " << genevt->get_collision_vertex().t() << std::endl;
      evt->print(*fout);
    }
  }
  return 0;
}
