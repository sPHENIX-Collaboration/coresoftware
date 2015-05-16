#include "Jet.h"

#include <cmath>

using namespace std;

ClassImp(Jet);

Jet::Jet()
  : _id(0xFFFFFFFF),
    _mom(),
    _e(NAN),
    _comp_ids()
{
  for (int i = 0; i < 3; ++i) _mom[i] = NAN;
}

void Jet::identify(ostream& os) const {
  os << "---Jet-----------------------" << endl;
  os << "jetid: " << get_id() << endl;
  os << " (px,py,pz,e) =  (" << get_px() << ", " << get_py() << ", ";
  os << get_pz() << ", " << get_e() << ") GeV" << endl;
  for (ConstIter citer = begin_comp(); citer != end_comp(); ++citer) {
    cout << citer->first << " -> " << citer->second << endl;
  }
  os << "-----------------------------------------------" << endl;
  
  return;  
}

void Jet::Reset() {
  _id = 0xFFFFFFFF;
  for (int i = 0; i < 3; ++i) _mom[i] = NAN;
  _e = NAN;
  _comp_ids.clear();  
}

int Jet::IsValid() const {
  if (_id == 0xFFFFFFFF) return 0;
  for (int i = 0; i < 3; ++i) {
    if (isnan(_mom[i])) return 0;
  }
  if (isnan(_e)) return 0;
  if (_comp_ids.empty()) return 0;
  return 1;
}
