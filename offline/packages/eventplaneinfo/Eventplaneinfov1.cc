#include "Eventplaneinfov1.h"

#include <cmath>
#include <limits>

void Eventplaneinfov1::identify(std::ostream& os) const
{
  os << "---------Eventplaneinfov1------------------" << std::endl;
  return;
}

double Eventplaneinfov1::GetPsi(const double Qx, const double Qy, const unsigned int order) const
{
  if (order == 0)
  {
    return std::numeric_limits<double>::quiet_NaN();
  }
  if ((Qx == 0.0) && (Qy == 0.0))
  {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return std::atan2(Qy, Qx) / order;
}
