/*!
 * \file TpcClusterZCrossingCorrection.cc
 * \brief applies correction to TPC cluster Z for bunch crossing time offset
 * \author Tony Frawley, March 2022 
 */

#include "TpcClusterZCrossingCorrection.h"

#include <cmath>
#include <iostream>
#include <climits>

float TpcClusterZCrossingCorrection::correctZ(float zinit, short int crossing)
{
  if(crossing == SHRT_MAX) return NAN;

  double vdrift = 8.00;  // cm /microsecond
  double z_bunch_separation = 0.1064 * vdrift;  // 106.4 ns bunch crossing interval

  // assume measured z is in the correct TPC side
  // negative crossings/times are in the past
  //    +ve crossing occurs in the future relative to time zero
  //        +ve z, cluster arrives late, so z seems smaller (more negative)
  //        -ve z, cluster arrives late, so z seems  more positive
  float corrected_z;
  if(zinit > 0)
    corrected_z = zinit - (float) crossing * z_bunch_separation;  
  else
    corrected_z = zinit + (float) crossing * z_bunch_separation;  
    
  std::cout << "         crossing " << crossing << " zinit " << zinit << " corrected_z " << corrected_z << std::endl;

  return corrected_z;
}


