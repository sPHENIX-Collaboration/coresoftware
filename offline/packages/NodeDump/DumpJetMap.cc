#include "DumpJetMap.h"

#include <g4jets/Jet.h>
#include <g4jets/JetMap.h>

#include <phool/PHIODataNode.h>

#include <climits>
#include <map>
#include <ostream>
#include <string>
#include <utility>

using MyNode_t = PHIODataNode<JetMap>;

DumpJetMap::DumpJetMap(const std::string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpJetMap::process_Node(PHNode *myNode)
{
  JetMap *jetmap = nullptr;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    jetmap = thisNode->getData();
  }
  if (jetmap)
  {
    JetMap::ConstIter jiter_beg = jetmap->begin();
    JetMap::ConstIter jiter_end = jetmap->end();
    *fout << "size: " << jetmap->size() << std::endl;
    *fout << "par: " << jetmap->get_par() << std::endl;
    *fout << "algo: " << jetmap->get_algo() << std::endl;
    for (JetMap::ConstIter jiter = jiter_beg; jiter != jiter_end; ++jiter)
    {
      *fout << "id: " << jiter->second->get_id() << std::endl;
      *fout << "px: " << jiter->second->get_px() << std::endl;
      *fout << "py: " << jiter->second->get_py() << std::endl;
      *fout << "pz: " << jiter->second->get_pz() << std::endl;
      *fout << "e: " << jiter->second->get_e() << std::endl;
      *fout << "p: " << jiter->second->get_p() << std::endl;
      *fout << "pt: " << jiter->second->get_pt() << std::endl;
      *fout << "et: " << jiter->second->get_et() << std::endl;
      *fout << "eta: " << jiter->second->get_eta() << std::endl;
      *fout << "phi: " << jiter->second->get_phi() << std::endl;
      *fout << "mass: " << jiter->second->get_mass() << std::endl;
      for (auto ic = 0; ic < UCHAR_MAX; ic++)
      {
        Jet::PROPERTY prop_id = static_cast<Jet::PROPERTY>(ic);
        if (jiter->second->has_property(prop_id))
        {
          *fout << "prop id: " << static_cast<unsigned int>(ic)
                << " value: " << jiter->second->get_property(prop_id)
                << std::endl;
        }
      }
      Jet::ConstIter jetbegin = jiter->second->begin_comp();
      Jet::ConstIter jetend = jiter->second->end_comp();
      for (Jet::ConstIter jetiter = jetbegin; jetiter != jetend; ++jetiter)
      {
        *fout << "src: " << jetiter->first << " value: " << jetiter->second << std::endl;
      }
    }
  }
  return 0;
}
