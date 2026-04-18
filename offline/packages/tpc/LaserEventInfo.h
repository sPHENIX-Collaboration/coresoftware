#ifndef TPC_LASEREVENTINFO_H
#define TPC_LASEREVENTINFO_H

//==================================
/// \file LaserEventInfo.h
/// \brief Laser event information storage node
/// \author Ben Kimelman
//==================================

#include <phool/PHObject.h>

#include <iostream>
#include <limits>

class LaserEventInfo : public PHObject
{
 public:
  ~LaserEventInfo() override{};

  void identify(std::ostream &os = std::cout) const override {os << "LaserEventInfo base class" << std::endl; };
  virtual void CopyTo(LaserEventInfo *) {return;}

  virtual bool isLaserEvent() const { return false; }
  virtual void setIsLaserEvent(const bool /*isLaserEvent*/) {}

  virtual bool isGl1LaserEvent() const { return false; }
  virtual void setIsGl1LaserEvent(const bool /*isLaserEvent*/) {}

  virtual bool isGl1LaserPileupEvent() const { return false; }
  virtual void setIsGl1LaserPileupEvent(const bool /*isLaserEvent*/) {}

  virtual int getPeakSample(const bool /*side*/) const { return std::numeric_limits<int>::max(); }
  virtual void setPeakSample(const bool /*side*/, const int /*sample*/) {}

  virtual float getPeakWidth(const bool /*side*/) const { return std::numeric_limits<float>::quiet_NaN(); }
  virtual void setPeakWidth(const bool /*side*/, const float /*width*/) {}

 protected:
  LaserEventInfo() = default;
  ClassDefOverride(LaserEventInfo, 1);
};

#endif
