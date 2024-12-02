#ifndef __MBDRUNNINGSTATS_H__
#define __MBDRUNNINGSTATS_H__
 
#include <queue>

/**
 * Class to calculate running average and RMS
 */
class MbdRunningStats
{
public:
  MbdRunningStats(const unsigned int imaxnum = 100);
  void Clear();
  void Push(double x);

  unsigned int Size() const { return static_cast<unsigned int>(values.size()); }
  unsigned int MaxNum() const { return maxnum; }

  double Mean() const;
  double Variance() const;
  double StandardDeviation() const;
  double RMS() const;

private:
  std::queue<int> values;
  unsigned int maxnum; // max values in sum
  double S1{0.};  // sum of values
  double S2{0.};  // sum of squares
};

#endif // __MBDRUNNINGSTATS_H__
