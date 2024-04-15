#ifndef __MBDRUNNINGSTATS_H__
#define __MBDRUNNINGSTATS_H__
 
#include <queue>

/**
 * Class to calculate running average and RMS
 */
class MbdRunningStats
{
public:
  MbdRunningStats(const std::size_t imaxnum = -1);
  void Clear();
  void Push(double x);

  long long NumDataValues() const;
  std::size_t MaxNum() const;

  double Mean() const;
  double Variance() const;
  double StandardDeviation() const;
  double RMS() const;


private:
  std::queue<int> values;
  std::size_t maxnum; // max values in sum
  double S1{0.};  // sum of values
  double S2{0.};  // sum of squares
};

#endif // __MBDRUNNINGSTATS_H__
