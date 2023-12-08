#include "DumpJetContainer.h"

#include <jetbase/Jet.h>
#include <jetbase/JetContainer.h>

#include <phool/PHIODataNode.h>

#include <climits>
#include <map>
#include <ostream>
#include <string>
#include <utility>

using MyNode_t = PHIODataNode<JetContainer>;

DumpJetContainer::DumpJetContainer(const std::string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpJetContainer::process_Node(PHNode *myNode)
{
  JetContainer *jets = nullptr;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    jets = thisNode->getData();
  }
  if (jets)
  {
    *fout << "size: " << jets->size() << std::endl;
    *fout << "par: " << jets->get_par() << std::endl;
    *fout << "algo: " << jets->get_algo() << std::endl;
    /* auto jet_properties = jets->get_property_vec(); */
    for (auto jet : *jets)
    {
      *fout << "id: " << jet->get_id() << std::endl;
      *fout << "px: " << jet->get_px() << std::endl;
      *fout << "py: " << jet->get_py() << std::endl;
      *fout << "pz: " << jet->get_pz() << std::endl;
      *fout << "e: " << jet->get_e() << std::endl;
      *fout << "p: " << jet->get_p() << std::endl;
      *fout << "pt: " << jet->get_pt() << std::endl;
      *fout << "et: " << jet->get_et() << std::endl;
      *fout << "eta: " << jet->get_eta() << std::endl;
      *fout << "phi: " << jet->get_phi() << std::endl;
      *fout << "mass: " << jet->get_mass() << std::endl;
      // print out the jet properties

      for (auto prop : jets->property_indices()) {
          *fout << "prop id: " << static_cast<unsigned int>(prop.first)
                << " value: " << jet->get_property(prop.second)
                << std::endl;
      }
      /* Jet::ConstIter jetbegin = jiter->second->begin_comp(); */
      /* Jet::ConstIter jetend = jiter->second->end_comp(); */
      /* for (Jet::ConstIter jetiter = jetbegin; jetiter != jetend; ++jetiter) */
      /* { */
      for (const auto& comp : jet->get_comp_vec()) {
        *fout << "src: " << comp.first << " value: " << comp.second << std::endl;
      }
    }
  }
  return 0;
}
