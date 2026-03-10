#ifndef TPC_LASEREVENTINFOV1_H
#define TPC_LASEREVENTINFOV1_H

#include "LaserEventInfo.h"

#include <iostream>

class LaserEventInfov1 : public LaserEventInfo
{
 public:
  LaserEventInfov1() = default;
  ~LaserEventInfov1();

  void identify(std::ostream &os = std::cout) const override;
  void Reset() override;

  PHObject* CloneMe() const override { return new LaserEventInfov1(*this); }
  void CopyTo(LaserEventInfo *info) override;

  bool isLaserEvent() const override { return m_isLaserEvent; }
  void setIsLaserEvent(const bool isLaserEvent) override { m_isLaserEvent = isLaserEvent; }

  int getPeakSample(const bool side) const override { return m_peakSample[side]; }
  void setPeakSample(const bool side, const int sample) override { m_peakSample[side] = sample; }

  double getPeakWidth(const bool side) const override { return m_peakWidth[side]; }
  void setPeakWidth(const bool side, const double width) override { m_peakWidth[side] = width; }

 protected:

  bool m_isLaserEvent{false};
  int m_peakSample[2] = {std::numeric_limits<int>::max(), std::numeric_limits<int>::max()};
  double m_peakWidth[2] = {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};

  ClassDefOverride(LaserEventInfov1, 1);
};

#endif
