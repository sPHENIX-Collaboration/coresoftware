#include "Eventplaneinfov2.h"

#include <cmath>

void Eventplaneinfov2::identify(std::ostream& os) const
{
  os << "---------Eventplaneinfov2------------------" << std::endl;
  return;
}

double Eventplaneinfov2::GetPsi(const double Qx, const double Qy, const unsigned int order) const
{
  if (order == 0)
  {
    return NAN;
  }
  if ((Qx == 0.0) && (Qy == 0.0))
  {
    return NAN;
  }
  return atan2(Qy, Qx) / static_cast<double>(order);
}
