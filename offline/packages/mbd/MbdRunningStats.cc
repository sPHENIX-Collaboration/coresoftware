#include "MbdRunningStats.h"
#include <limits>
#include <cmath>
 
MbdRunningStats::MbdRunningStats(const std::size_t imaxnum) :
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
  if ( values.size() == maxnum )
  {
    double lastval = values.front();
    values.pop();
    S1 -= lastval;
    S2 -= lastval*lastval;
  }

  values.push(x);
  S1 += x;
  S2 += x*x;
}

long long MbdRunningStats::NumDataValues() const
{
  return values.size();
}

double MbdRunningStats::Mean() const
{
  if ( values.size()==0 ) 
  {
    //return std::numeric_limits<double>::infinity();
    return std::numeric_limits<float>::quiet_NaN();
  }

  return S1/values.size();
}

double MbdRunningStats::Variance() const
{
  if ( values.size()==0 ) 
  {
    //return std::numeric_limits<double>::infinity();
    return std::numeric_limits<float>::quiet_NaN();
  }
  else if ( values.size()==1 )
  {
    return 0.;
  }

  double var = S2/(values.size()-1.0) - Mean()*Mean();
  if ( var>4 )
  {
    var = 4.;
  }

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

