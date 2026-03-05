#include "LaserEventInfov1.h"

LaserEventInfov1::~LaserEventInfov1() = default;

void LaserEventInfov1::identify(std::ostream& os) const
{
  os << "LaserEventInfov1: " << std::endl;
  os << "   isLaserEvent? " << isLaserEvent() << std::endl;
  if (isLaserEvent())
  {
    os << "      South peak sample: " << getPeakSample(false) << "   width: " << getPeakWidth(false) << std::endl;
    os << "      North peak sample: " << getPeakSample(true) << "   width: " << getPeakWidth(true) << std::endl;
  }
  return;
}

void LaserEventInfov1::Reset()
{
  m_isLaserEvent = false;
  for (int i = 0; i < 2; i++)
  {
    m_peakSample[i] = std::numeric_limits<int>::max();
    m_peakWidth[i] = std::numeric_limits<double>::quiet_NaN();
  }

  return;
}

void LaserEventInfov1::CopyTo(LaserEventInfo* info)
{
  info->setIsLaserEvent(isLaserEvent());

  for (int side = 0; side <= 1; side++)
  {
    info->setPeakSample(side, getPeakSample(side));
    info->setPeakWidth(side, getPeakWidth(side));
  }

  return;
}
