#include "TowerBackground_v1.h"

using namespace std;

TowerBackground_v1::TowerBackground_v1()
{

  _v2[0] = 0.0;
  _v2[1] = 0.0;
  _v2[2] = 0.0;

  _Psi2[0] = 0.0;
  _Psi2[1] = 0.0;
  _Psi2[2] = 0.0;
  
  _UE.resize(3, std::vector<float>(1, 0) );

}

TowerBackground_v1::~TowerBackground_v1()
{
}

void TowerBackground_v1::identify(ostream& os) const
{

  os << "TowerBackground: " << std::endl;
  for (int n = 0; n < 3; n++) {
    
    os << " layer " << n << " : UE in " << _UE[ n ].size() << " eta bins: ";
    for (unsigned int eta = 0; eta < _UE[ n ].size(); eta++) 
      os << _UE[ n ].at( eta ) << " ";
    os << std::endl;
    os << " layer " << n << " : v2 = " << _v2[ n ] << ", Psi2 = " << _Psi2[ n ] << std::endl;
  }

  return;
}
