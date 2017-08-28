
#include "PHFieldConst.h"

//root framework
#include <TFile.h>
#include <TNtuple.h>

#include <iostream>
#include <set>

using namespace std;
using namespace CLHEP;  // units

PHFieldConst::PHFieldConst(
    double field_mag_x,
    double field_mag_y,
    double field_mag_z):
    field_mag_x_(field_mag_x),
    field_mag_y_(field_mag_y),
    field_mag_z_(field_mag_z)
{
}

void PHFieldConst::GetFieldValue(const double point[4], double *Bfield) const
{
  Bfield[0] = field_mag_x_;
  Bfield[1] = field_mag_y_;
  Bfield[2] = field_mag_z_;

  return;
}
