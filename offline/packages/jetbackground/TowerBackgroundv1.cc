#include "TowerBackgroundv1.h"

#include <memory>
#include <ostream>

TowerBackgroundv1::TowerBackgroundv1()
{
  _UE.resize(3, std::vector<float>(1, 0));
}

void TowerBackgroundv1::identify(std::ostream& os) const
{
  os << "TowerBackground: " << std::endl;
  for (int n = 0; n < 3; n++)
  {
    os << " layer " << n << " : UE in " << _UE[n].size() << " eta bins: ";
    for (float eta : _UE[n])
    {
      os << eta << " ";
    }
    os << std::endl;
  }

  os << " v2 = " << _v2 << ", Psi2 = " << _Psi2 << ", # towers used for bkg = " << _nTowers << " , # strips used for flow = " << _nStrips << std::endl;

  return;
}
