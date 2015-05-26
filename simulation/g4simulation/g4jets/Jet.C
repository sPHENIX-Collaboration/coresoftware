#include "Jet.h"

using namespace std;

ClassImp(Jet);

Jet::Jet() {}

void Jet::identify(ostream& os) const {
  os << "---Jet-----------------------" << endl;
  return;
}
