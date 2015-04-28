#include "RawTower.h"

using namespace std;

std::pair< std::map<unsigned int,float>::const_iterator, std::map<unsigned int,float>::const_iterator >
RawTower::get_g4cells()
{
  map <unsigned int, float> dummy;
  return make_pair(dummy.begin(), dummy.end());
}
