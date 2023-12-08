#include "JetContainer.h"
#include <ostream>  // for operator<<, endl, ostream, basic_ostream

class Jet;

// Dummy members for virtual classes to return in base implementation
Jet* DummyJetV2;
std::map<Jet::PROPERTY, Jet::PROPERTY> DummyMap;
std::vector<float> DummyFloatVec;
TClonesArray DummyTClonesArray("Jet", 1);
std::set<Jet::SRC> DummyJetContSourceSet;
// N.B. this will seg fault if used... so don't use it!

std::map<Jet::PROPERTY, Jet::PROPERTY /*really this is an index*/>& property_indices() { return DummyMap; }

JetContainer::ConstSrcIter JetContainer::begin_src() const
{
  return DummyJetContSourceSet.end();
}

JetContainer::ConstSrcIter JetContainer::find_src(Jet::SRC /*src*/) const
{
  return DummyJetContSourceSet.end();
}

JetContainer::ConstSrcIter JetContainer::end_src() const
{
  return DummyJetContSourceSet.end();
}

JetContainer::SrcIter JetContainer::begin_src()
{
  return DummyJetContSourceSet.end();
}

JetContainer::SrcIter JetContainer::find_src(Jet::SRC /*src*/)
{
  return DummyJetContSourceSet.end();
}

JetContainer::SrcIter JetContainer::end_src()
{
  return DummyJetContSourceSet.end();
}

void JetContainer::identify(std::ostream& os) const
{
  os << "JetContainer" << std::endl;
}

TClonesArray* JetContainer::clone_data() const
{
  return (TClonesArray*) DummyTClonesArray.Clone();
};

Jet::IterJetTCA DummyIterJetTCA{&DummyTClonesArray};

/* Jet* JetContainer::current_jet() { return DummyJetV2; }; */

/* std::vector<float>& JetContainer::jet_properties() { return DummyFloatVec; }; */

Jet::IterJetTCA JetContainer::begin()
{
  return DummyIterJetTCA;
}

Jet::IterJetTCA JetContainer::end()
{
  return DummyIterJetTCA;
}
