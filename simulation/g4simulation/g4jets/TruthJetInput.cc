
#include "TruthJetInput.h"

#include "Jet.h"
#include "Jetv1.h"

#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>

#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

// standard includes
#include <algorithm>  // std::find
#include <cmath>      // for asinh, sqrt
#include <cstdlib>
#include <iostream>
#include <map>      // for _Rb_tree_const_iterator
#include <utility>  // for pair
#include <vector>

TruthJetInput::TruthJetInput(Jet::SRC input)
  : m_Input(input)
{
}

void TruthJetInput::identify(std::ostream &os)
{
  os << "   TruthJetInput: G4TruthInfo to Jet::PARTICLE";
  if (use_embed_stream())
  {
    os << ". Processing embedded streams: ";
    for (int it : m_EmbedID)
    {
      os << it << ", ";
    }
  }
  os << std::endl;
}

std::vector<Jet *> TruthJetInput::get_input(PHCompositeNode *topNode)
{
  if (Verbosity() > 0) std::cout << "TruthJetInput::process_event -- entered" << std::endl;

  // Pull the reconstructed track information off the node tree...
  PHG4TruthInfoContainer *truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!truthinfo)
  {
    std::cout << PHWHERE << " ERROR: Can't find G4TruthInfo" << std::endl;
    return std::vector<Jet *>();
  }

  std::vector<Jet *> pseudojets;
  PHG4TruthInfoContainer::ConstRange range = truthinfo->GetPrimaryParticleRange();
  for (PHG4TruthInfoContainer::ConstIterator iter = range.first;
       iter != range.second;
       ++iter)
  {
    PHG4Particle *part = iter->second;

    if (use_embed_stream())
    {
      const int this_embed_id = truthinfo->isEmbeded(part->get_track_id());

      if (std::find(m_EmbedID.begin(), m_EmbedID.end(), this_embed_id) == m_EmbedID.end())
      {
        continue;  // reject particle as it is not in the interested embedding stream.
      }
    }

    // remove some particles (muons, taus, neutrinos)...
    // 12 == nu_e
    // 13 == muons
    // 14 == nu_mu
    // 15 == taus
    // 16 == nu_tau
    if ((abs(part->get_pid()) >= 12) && (abs(part->get_pid()) <= 16)) continue;

    // remove acceptance... _etamin,_etamax
    if ((part->get_px() == 0.0) && (part->get_py() == 0.0)) continue;  // avoid pt=0
    float eta = asinh(part->get_pz() / sqrt(pow(part->get_px(), 2) + pow(part->get_py(), 2)));
    if (eta < m_EtaMin) continue;
    if (eta > m_EtaMax) continue;

    Jet *jet = new Jetv1();
    jet->set_px(part->get_px());
    jet->set_py(part->get_py());
    jet->set_pz(part->get_pz());
    jet->set_e(part->get_e());
    jet->insert_comp(Jet::PARTICLE, part->get_track_id());
    pseudojets.push_back(jet);
  }

  if (Verbosity() > 0) std::cout << "TruthJetInput::process_event -- exited" << std::endl;

  return pseudojets;
}
