#include "DumpJetMap.h"

#include <g4jets/Jet.h>
#include <g4jets/JetMap.h>

#include <phool/PHIODataNode.h>

#include <climits>
#include <map>
#include <ostream>
#include <string>
#include <utility>

using namespace std;

typedef PHIODataNode<JetMap> MyNode_t;

DumpJetMap::DumpJetMap(const string &NodeName)
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
    *fout << "size: " << jetmap->size() << endl;
    *fout << "par: " << jetmap->get_par() << endl;
    *fout << "algo: " << jetmap->get_algo() << endl;
    for (JetMap::ConstIter jiter = jiter_beg; jiter != jiter_end; ++jiter)
    {
      *fout << "id: " << jiter->second->get_id() << endl;
      *fout << "px: " << jiter->second->get_px() << endl;
      *fout << "py: " << jiter->second->get_py() << endl;
      *fout << "pz: " << jiter->second->get_pz() << endl;
      *fout << "e: " << jiter->second->get_e() << endl;
      *fout << "p: " << jiter->second->get_p() << endl;
      *fout << "pt: " << jiter->second->get_pt() << endl;
      *fout << "et: " << jiter->second->get_et() << endl;
      *fout << "eta: " << jiter->second->get_eta() << endl;
      *fout << "phi: " << jiter->second->get_phi() << endl;
      *fout << "mass: " << jiter->second->get_mass() << endl;
      for (unsigned char ic = 0; ic < UCHAR_MAX; ic++)
      {
        Jet::PROPERTY prop_id = static_cast<Jet::PROPERTY>(ic);
        if (jiter->second->has_property(prop_id))
        {
          *fout << "prop id: " << static_cast<unsigned int>(ic)
                << " value: " << jiter->second->get_property(prop_id)
                << endl;
        }
      }
      Jet::ConstIter jetbegin = jiter->second->begin_comp();
      Jet::ConstIter jetend = jiter->second->end_comp();
      for (Jet::ConstIter jetiter = jetbegin; jetiter != jetend; ++jetiter)
      {
        *fout << "src: " << jetiter->first << " value: " << jetiter->second << endl;
      }
    }
  }
  return 0;
}
