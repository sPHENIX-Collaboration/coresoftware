#include "DumpParticleFlowElementContainer.h"

#include <phool/PHIODataNode.h>

#include <particleflowreco/ParticleFlowElement.h>
#include <particleflowreco/ParticleFlowElementContainer.h>

#include <map>
#include <ostream>
#include <string>
#include <utility>

using MyNode_t = PHIODataNode<ParticleFlowElementContainer>;

DumpParticleFlowElementContainer::DumpParticleFlowElementContainer(const std::string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpParticleFlowElementContainer::process_Node(PHNode *myNode)
{
  ParticleFlowElementContainer *particleflowelementcontainer = nullptr;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    particleflowelementcontainer = thisNode->getData();
  }
  if (particleflowelementcontainer)
  {
    ParticleFlowElementContainer::ConstIterator hiter;
    ParticleFlowElementContainer::ConstRange begin_end = particleflowelementcontainer->getParticleFlowElements();
    *fout << "size: " << particleflowelementcontainer->size() << std::endl;
    for (hiter = begin_end.first; hiter != begin_end.second; ++hiter)
    {
      ParticleFlowElement *pfe = hiter->second;
      *fout << "get_id(): " << pfe->get_id() << std::endl;
      *fout << "get_type(): " << pfe->get_type() << std::endl;
      *fout << "get_px(): " << pfe->get_px() << std::endl;
      *fout << "get_py(): " << pfe->get_py() << std::endl;
      *fout << "get_pz(): " << pfe->get_pz() << std::endl;
      *fout << "get_e(): " << pfe->get_e() << std::endl;
      *fout << "get_p(): " << pfe->get_p() << std::endl;
      *fout << "get_pt(): " << pfe->get_pt() << std::endl;
      *fout << "get_et(): " << pfe->get_et() << std::endl;
      *fout << "get_eta(): " << pfe->get_eta() << std::endl;
      *fout << "get_phi(): " << pfe->get_phi() << std::endl;
      *fout << "get_mass(): " << pfe->get_mass() << std::endl;
    }
  }
  return 0;
}
