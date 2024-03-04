#include "Eventplaneinfov1.h"

#include <cmath>

void Eventplaneinfov1::identify(std::ostream& os) const
{
  os << "---------Eventplaneinfov1------------------" << std::endl;
  os << "\t second order event plane angle is " << get_psi(2) << std::endl;
  return;
}

double Eventplaneinfov1::GetPsi(double Qx, double Qy, unsigned int order) const
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
