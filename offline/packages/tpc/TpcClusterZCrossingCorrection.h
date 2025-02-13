#ifndef TPC_TPCCLUSTERZCROSSINGCORRECTION_H
#define TPC_TPCCLUSTERZCROSSINGCORRECTION_H

/*!
 * \file TpcClusterZCrossingCorrection.h
 * \brief applies correction to TPC cluster Z for bunch crossing time offset
 * \author Tony Frawley, March 2022
 */

class TpcClusterZCrossingCorrection
{
  public:

  TpcClusterZCrossingCorrection() = default;

  static float correctZ(float zinit, unsigned int side, short int crossing);

  static float _vdrift;

  static float _time_between_crossings;  // ns, same value as in pileup generator

};

#endif
