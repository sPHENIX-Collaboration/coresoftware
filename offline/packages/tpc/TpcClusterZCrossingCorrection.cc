/*!
 * \file TpcClusterZCrossingCorrection.cc
 * \brief applies correction to TPC cluster Z for bunch crossing time offset
 * \author Tony Frawley, March 2022 
 */

#include "TpcClusterZCrossingCorrection.h"

#include <cmath>
#include <iostream>
#include <climits>

float TpcClusterZCrossingCorrection::_vdrift = 8.0e-03;  // default value, override from macro

TpcClusterZCrossingCorrection::TpcClusterZCrossingCorrection()
{ }

float TpcClusterZCrossingCorrection::correctZ(float zinit, unsigned int side, short int crossing) const
{
  if(crossing == SHRT_MAX) return NAN;

  float z_bunch_separation = _time_between_crossings * _vdrift;  

  //    +ve crossing occurs in the future relative to time zero
  //        -ve z side (south, side 0), cluster arrives late, so z seems more positive
  //        +ve z side (north, side 1), cluster arrives late, so z seems more negative
  float corrected_z;
  if(side == 0) 
    corrected_z = zinit - (float) crossing * z_bunch_separation;  
  else
    corrected_z = zinit + (float) crossing * z_bunch_separation;  
    
  //  std::cout << "         crossing " << crossing << " _vdrift " << _vdrift << " zinit " << zinit << " side " << side << " corrected_z " << corrected_z << std::endl;

  return corrected_z;
}


