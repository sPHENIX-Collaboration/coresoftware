#include "LaserEventInfov2.h"

LaserEventInfov2::~LaserEventInfov2() = default;

void LaserEventInfov2::identify(std::ostream& os) const
{
  os << "LaserEventInfov2: " << std::endl;
  os << "   isGl1LaserEvent? " << isGl1LaserEvent() << std::endl;
  os << "   isGl1LaserPileupEvent? " << isGl1LaserPileupEvent() << std::endl;
  os << "   isLaserEvent? " << isLaserEvent() << std::endl;
  if (isLaserEvent())
  {
    os << "      South peak sample: " << getPeakSample(false) << "   width: " << getPeakWidth(false) << std::endl;
    os << "      North peak sample: " << getPeakSample(true) << "   width: " << getPeakWidth(true) << std::endl;
  }

  return;
}

void LaserEventInfov2::Reset()
{
  m_isLaserEvent = false;
  m_isGl1LaserEvent = false;
  for (int i = 0; i < 2; i++)
  {
    m_peakSample[i] = std::numeric_limits<int>::max();
    m_peakWidth[i] = std::numeric_limits<float>::quiet_NaN();
  }

  return;
}

void LaserEventInfov2::CopyTo(LaserEventInfo* info)
{
  info->setIsLaserEvent(isLaserEvent());
  info->setIsGl1LaserEvent(isGl1LaserEvent());
  info->setIsGl1LaserPileupEvent(isGl1LaserPileupEvent());

  for (int side = 0; side <= 1; side++)
  {
    info->setPeakSample(side, getPeakSample(side));
    info->setPeakWidth(side, getPeakWidth(side));
  }

  return;
}
