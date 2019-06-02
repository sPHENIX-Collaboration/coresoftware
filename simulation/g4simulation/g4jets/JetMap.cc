#include "JetMap.h"

#include <ostream>  // for operator<<, endl, ostream, basic_ostream

using namespace std;

JetMap::JetMap() {}

void JetMap::identify(ostream& os) const
{
  os << "JetMap" << endl;
  return;
}
