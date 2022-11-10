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
  
TpcClusterZCrossingCorrection();

  float correctZ(float zinit, unsigned int side, short int crossing) const;

  static float _vdrift;

private:

  float _time_between_crossings = 106; // ns, same value as in pileup generator  
};

#endif
