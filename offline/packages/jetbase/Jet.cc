#include "Jet.h"

#include <iostream>
#include <vector>
#include <map>
#include <TClonesArray.h>

void Jet::mv1(const std::string& method_name, std::ostream& os) {
    os << "Method Jet::" << method_name << " implemented only in Jetv1." << std::endl;
};
void Jet::mv2(const std::string& method_name, std::ostream& os) {
    os << "Method Jet::" << method_name << " implemented only in Jetv2 and above." << std::endl;
};


Jet::typ_comp_ids  DummyJetMap;
std::vector<float> DummyJetPropVec;
Jet::TYPE_comp_vec DummyJetCompVec;

std::vector<float>& Jet::get_vec_properties()
{ mv2("get_vec_properties()"); return DummyJetPropVec; }

// functions interfaces only for Jetv2 and above
Jet::ITER_comp_vec Jet::comp_begin(Jet::SRC /**/) { mv2("comp_begin()"); return DummyJetCompVec.end(); }
Jet::ITER_comp_vec Jet::comp_end(Jet::SRC /**/)   { mv2("comp_end()"  ); return DummyJetCompVec.end(); }
Jet::ITER_comp_vec Jet::comp_begin() { mv2("comp_begin()"); return DummyJetCompVec.end(); }
Jet::ITER_comp_vec Jet::comp_end()   { mv2("comp_end()"  ); return DummyJetCompVec.end(); }

Jet::TYPE_comp_vec& Jet::get_comp_vec() { mv2("get_comp_vec()"); return DummyJetCompVec; }

void Jet::identify(std::ostream& os) const
{
  os << "---Jet-----------------------" << std::endl;
  return;
}

// functions interfaces only in Jetv1
Jet::ConstIter Jet::begin_comp()       const { mv1("begin_comp()"); return DummyJetMap.end(); }
Jet::ConstIter Jet::end_comp()         const { mv1("end_comp()");   return DummyJetMap.end(); }
Jet::ConstIter Jet::find(Jet::SRC/**/) const { mv1("find()");       return DummyJetMap.end(); }

Jet::ConstIter Jet::lower_bound_comp(Jet::SRC/**/) const { mv1("lower_bound_comp()"); return DummyJetMap.end(); }
Jet::ConstIter Jet::upper_bound_comp(Jet::SRC/**/) const { mv1("upper_bound_comp()"); return DummyJetMap.end(); }


Jet::Iter Jet::begin_comp()       { mv1("begin_comp()"); return DummyJetMap.end(); }
Jet::Iter Jet::end_comp()         { mv1("end_comp()");   return DummyJetMap.end(); }
Jet::Iter Jet::find(Jet::SRC/**/) { mv1("find()");       return DummyJetMap.end(); }

Jet::Iter Jet::lower_bound_comp(Jet::SRC /**/) { mv1("lower_bound_comp()"); return DummyJetMap.end(); }
Jet::Iter Jet::upper_bound_comp(Jet::SRC /**/) { mv1("upper_bound_comp()"); return DummyJetMap.end(); }


// structure to iterate over ther jets in a TClonesArray in the JetContainer
Jet::IterJetTCA::IterJetTCA(TClonesArray* _tca, Jet*& _in_jet)
        : tca{_tca}, current_jet{_in_jet}, size { tca->GetEntriesFast() }
{ current_jet = (Jet*) tca->UncheckedAt(0); }

void Jet::IterJetTCA::operator++() { 
  ++index;
  current_jet = (Jet*) tca->UncheckedAt(index); 
}

bool Jet::IterJetTCA::operator!=(const IterJetTCA& rhs) {
  if (index == rhs.size) {
    current_jet = (Jet*) tca->UncheckedAt(0);
    return false;
  } else {
    return true;
  }
};
