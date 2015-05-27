#include "JetMap.h"

#include "Jet.h"

#include <cmath>

using namespace std;

ClassImp(JetMap)

JetMap::JetMap() {}

void JetMap::identify(ostream& os) const {
  os << "JetMap" << endl;
  return;
}
