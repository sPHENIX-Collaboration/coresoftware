#include "TrackJetInput.h"

#include "Jet.h"
#include "Jetv1.h"

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>

#include <phool/getClass.h>

// standard includes
#include <iostream>
#include <map>      // for _Rb_tree_const_iterator
#include <utility>  // for pair
#include <vector>

TrackJetInput::TrackJetInput(Jet::SRC input, const std::string &nodename)
  : m_NodeName(nodename)
  , _input(input)
{
}

void TrackJetInput::identify(std::ostream &os)
{
  os << "   TrackJetInput: SvtxTrackMap to Jet::TRACK" << std::endl;
}

std::vector<Jet *> TrackJetInput::get_input(PHCompositeNode *topNode)
{
  if (Verbosity() > 0) std::cout << "TrackJetInput::process_event -- entered" << std::endl;

  // Pull the reconstructed track information off the node tree...
  SvtxTrackMap *trackmap = findNode::getClass<SvtxTrackMap>(topNode, m_NodeName);
  if (!trackmap)
  {
    return std::vector<Jet *>();
  }

  std::vector<Jet *> pseudojets;
  for (SvtxTrackMap::ConstIter iter = trackmap->begin();
       iter != trackmap->end();
       ++iter)
  {
    const SvtxTrack *track = iter->second;

    Jet *jet = new Jetv1();
    jet->set_px(track->get_px());
    jet->set_py(track->get_py());
    jet->set_pz(track->get_pz());
    jet->set_e(track->get_p());
    jet->insert_comp(Jet::TRACK, track->get_id());
    pseudojets.push_back(jet);
  }

  if (Verbosity() > 0) std::cout << "TrackJetInput::process_event -- exited" << std::endl;

  return pseudojets;
}
