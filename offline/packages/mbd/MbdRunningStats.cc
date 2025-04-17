#include "MbdRunningStats.h"
#include <limits>
#include <cmath>
#include <iostream>
 
MbdRunningStats::MbdRunningStats(const unsigned int imaxnum) :
  maxnum{imaxnum}
{
  Clear();
}

void MbdRunningStats::Clear()
{
  while ( !values.empty() )
  {
    values.pop();
  }
  S1 = S2 = 0.0;
}

void MbdRunningStats::Push(double x)
{
  if ( Size() == maxnum )
  {
    double lastval = values.front();
    values.pop();
    S1 -= lastval;
    S2 -= (lastval*lastval);
  }

  values.push(x);
  S1 += x;
  S2 += (x*x);
}

double MbdRunningStats::Mean() const
{
  if ( values.empty() ) 
  {
    //return std::numeric_limits<double>::infinity();
    //return std::numeric_limits<float>::quiet_NaN();
    return 0.;
  }

  return S1/values.size();
}

double MbdRunningStats::Variance() const
{
  if ( values.empty() ) 
  {
    //return std::numeric_limits<double>::infinity();
    //return std::numeric_limits<float>::quiet_NaN();
    return 0.;
  }

  double var = (S2/values.size()) - (Mean()*Mean());
  /*
  std::cout << "RMS " << S2 << "\t" << values.size() << "\t" << Mean() << "\t" << Mean()*Mean() << "\t" << var << std::endl;
  */
  /*
  if ( var>4 )
  {
    var = 4.;
  }
  */

  return var;
}

double MbdRunningStats::StandardDeviation() const
{
  return sqrt( Variance() );
}

double MbdRunningStats::RMS() const
{
  return sqrt( Variance() );
}

