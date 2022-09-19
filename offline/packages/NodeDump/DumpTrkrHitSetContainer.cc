#include "DumpTrkrHitSetContainer.h"

#include <phool/PHIODataNode.h>

#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>

#include <map>
#include <ostream>
#include <string>
#include <utility>

using MyNode_t = PHIODataNode<TrkrHitSetContainer>;

DumpTrkrHitSetContainer::DumpTrkrHitSetContainer(const std::string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpTrkrHitSetContainer::process_Node(PHNode *myNode)
{
  TrkrHitSetContainer *trkrhitsetcontainer = nullptr;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    trkrhitsetcontainer = thisNode->getData();
  }
  if (trkrhitsetcontainer)
  {
    TrkrHitSetContainer::ConstIterator hiter;
    TrkrHitSetContainer::ConstRange begin_end = trkrhitsetcontainer->getHitSets();
    *fout << "size: " << trkrhitsetcontainer->size() << std::endl;
    for (hiter = begin_end.first; hiter != begin_end.second; ++hiter)
    {
      TrkrHitSet *trkrhitset = hiter->second;
      TrkrHitSet::ConstIterator tsetiter;
      TrkrHitSet::ConstRange trset_begin_end = trkrhitset->getHits();
      for (tsetiter = trset_begin_end.first; tsetiter != trset_begin_end.second; ++tsetiter)
      {
        TrkrHit *hit = tsetiter->second;
        *fout << "id: " << tsetiter->first << std::endl;
        *fout << "adc: " << hit->getAdc() << std::endl;
        *fout << "energy: " << hit->getEnergy() << std::endl;
      }
    }
  }
  return 0;
}
