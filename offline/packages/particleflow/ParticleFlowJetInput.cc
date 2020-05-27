#include "ParticleFlowJetInput.h"

#include <phool/getClass.h>

#include <g4jets/Jet.h>
#include <g4jets/Jetv1.h>

// standard includes
#include <cassert>
#include <iostream>
#include <map>                               // for _Rb_tree_const_iterator
#include <utility>                           // for pair
#include <vector>

#include "ParticleFlowElement.h"
#include "ParticleFlowElementContainer.h"

using namespace std;

ParticleFlowJetInput::ParticleFlowJetInput( )
  : _verbosity(0)
{
}

void ParticleFlowJetInput::identify(std::ostream &os)
{
  os << "   ParticleFlowJetInput: ";
  os << std::endl;
}

std::vector<Jet *> ParticleFlowJetInput::get_input(PHCompositeNode *topNode)
{
  if (_verbosity > 0) cout << "ParticleFlowJetInput::process_event -- entered" << endl;

  ParticleFlowElementContainer *pflowContainer = findNode::getClass<ParticleFlowElementContainer>(topNode, "ParticleFlowElements");
  if (!pflowContainer)
    {
      return std::vector<Jet *>();
    }

  std::vector<Jet *> pseudojets;
  ParticleFlowElementContainer::ConstRange begin_end = pflowContainer->getParticleFlowElements();
  ParticleFlowElementContainer::ConstIterator rtiter;
  for (rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter)
  {
    ParticleFlowElement *pflow = rtiter->second;

    Jet *jet = new Jetv1();
    jet->set_px( pflow->get_px() );
    jet->set_py( pflow->get_py() );
    jet->set_pz( pflow->get_pz() );
    jet->set_e( pflow->get_e() );
    jet->insert_comp( Jet::SRC::PARTICLE , pflow->get_id() );
    pseudojets.push_back( jet );
  }

  if (_verbosity > 0) cout << "ParticleFlowJetInput::process_event -- exited" << endl;

  return pseudojets;
}
