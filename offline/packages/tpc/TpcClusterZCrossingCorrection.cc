/*!
 * \file TpcClusterZCrossingCorrection.cc
 * \brief applies correction to TPC cluster Z for bunch crossing time offset
 * \author Tony Frawley, March 2022
 */

#include "TpcClusterZCrossingCorrection.h"

#include <phool/sphenix_constants.h>

#include <cmath>
#include <iostream>
#include <limits>

// default value, override from macro (cm/ns)
float TpcClusterZCrossingCorrection::_vdrift = 8.0e-03;

// ns, same value as in pileup generator
float TpcClusterZCrossingCorrection::_time_between_crossings = sphenix_constants::time_between_crossings;

//______________________________________________________________________________________________
float TpcClusterZCrossingCorrection::correctZ(float zinit, unsigned int side, short int crossing)
{
  if (crossing == std::numeric_limits<short>::max())
  {
    return std::numeric_limits<float>::quiet_NaN();
  }

  float z_bunch_separation = _time_between_crossings * _vdrift;

  //    +ve crossing occurs in the future relative to time zero
  //        -ve z side (south, side 0), cluster arrives late, so z seems more positive
  //        +ve z side (north, side 1), cluster arrives late, so z seems more negative
  float corrected_z;
  if (side == 0)
  {
    corrected_z = zinit - (float) crossing * z_bunch_separation;
  }
  else
  {
    corrected_z = zinit + (float) crossing * z_bunch_separation;
  }

  //  std::cout << "         crossing " << crossing << " _vdrift " << _vdrift << " zinit " << zinit << " side " << side << " corrected_z " << corrected_z << std::endl;

  return corrected_z;
}
