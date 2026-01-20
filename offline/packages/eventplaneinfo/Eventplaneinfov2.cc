#include "Eventplaneinfov2.h"

#include <cmath>

void Eventplaneinfov2::identify(std::ostream& os) const
{
  os << "---------Eventplaneinfov2------------------" << std::endl;
  return;
}

double Eventplaneinfov2::GetPsi(const double Qx, const double Qy, const unsigned int order) const
{
  double temp;
  if ((Qx == 0.0) && (Qy == 0.0))
  {
    temp = NAN;
  }
  else
  {
    temp = atan2(Qy, Qx) / ((double) order);
  }
  return temp;
}
