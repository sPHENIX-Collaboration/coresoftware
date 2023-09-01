#include "Jet.h"

#include <TClonesArray.h>
#include <iostream>
#include <map>
#include <vector>

void Jet::depmsg(const std::string& method_name, const std::string& version, std::ostream& os)
{
  os << "Method Jet::" << method_name << " implemented only in " << version << std::endl;
};

Jet::typ_comp_ids DummyJetMap;
std::vector<float> DummyJetPropVec;
Jet::TYPE_comp_vec DummyJetCompVec;

std::vector<float>& Jet::get_vec_properties()
{
  depmsg("get_vec_properties()",v2andup);
  return DummyJetPropVec;
}

// functions interfaces only for Jetv2 and above
Jet::ITER_comp_vec Jet::comp_begin(Jet::SRC /**/)
{
  depmsg("comp_begin()",v2andup);
  return DummyJetCompVec.end();
}
Jet::ITER_comp_vec Jet::comp_end(Jet::SRC /**/)
{
  depmsg("comp_end()",v2andup);
  return DummyJetCompVec.end();
}
Jet::ITER_comp_vec Jet::comp_begin()
{
  depmsg("comp_begin()",v2andup);
  return DummyJetCompVec.end();
}
Jet::ITER_comp_vec Jet::comp_end()
{
  depmsg("comp_end()",v2andup);
  return DummyJetCompVec.end();
}

Jet::TYPE_comp_vec& Jet::get_comp_vec()
{
  depmsg("get_comp_vec()",v2andup);
  return DummyJetCompVec;
}

void Jet::identify(std::ostream& os) const
{
  os << "---Jet-----------------------" << std::endl;
  return;
}

// functions interfaces only in Jetv1
Jet::ConstIter Jet::begin_comp() const
{
  depmsg("begin_comp()", v1only);
  return DummyJetMap.end();
}
Jet::ConstIter Jet::end_comp() const
{
  depmsg("end_comp()", v1only);
  return DummyJetMap.end();
}
Jet::ConstIter Jet::find(Jet::SRC /**/) const
{
  depmsg("find()", v1only);
  return DummyJetMap.end();
}

Jet::ConstIter Jet::lower_bound_comp(Jet::SRC /**/) const
{
  depmsg("lower_bound_comp()", v1only);
  return DummyJetMap.end();
}
Jet::ConstIter Jet::upper_bound_comp(Jet::SRC /**/) const
{
  depmsg("upper_bound_comp()", v1only);
  return DummyJetMap.end();
}

Jet::Iter Jet::begin_comp()
{
  depmsg("begin_comp()", v1only);
  return DummyJetMap.end();
}
Jet::Iter Jet::end_comp()
{
  depmsg("end_comp()", v1only);
  return DummyJetMap.end();
}
Jet::Iter Jet::find(Jet::SRC /**/)
{
  depmsg("find()", v1only);
  return DummyJetMap.end();
}

Jet::Iter Jet::lower_bound_comp(Jet::SRC /**/)
{
  depmsg("lower_bound_comp()", v1only);
  return DummyJetMap.end();
}
Jet::Iter Jet::upper_bound_comp(Jet::SRC /**/)
{
  depmsg("upper_bound_comp()", v1only);
  return DummyJetMap.end();
}

// structure to iterate over ther jets in a TClonesArray in the JetContainer
Jet::IterJetTCA::IterJetTCA(TClonesArray* _tca, Jet*& _in_jet)
  : tca{_tca}
  , current_jet{_in_jet}
  , size{tca->GetEntriesFast()}
{
  current_jet = (Jet*) tca->UncheckedAt(0);
}

void Jet::IterJetTCA::operator++()
{
  ++index;
  current_jet = (Jet*) tca->UncheckedAt(index);
}

bool Jet::IterJetTCA::operator!=(const IterJetTCA& rhs)
{
  if (index == rhs.size)
  {
    current_jet = (Jet*) tca->UncheckedAt(0);
    return false;
  }
  else
  {
    return true;
  }
};
