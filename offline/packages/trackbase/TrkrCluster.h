/**
 * @file trackbase/TrkrCluster.h
 * @author D. McGlinchey
 * @date June 2018
 * @brief Base class for cluster object
 */
#ifndef TRACKBASE_TRKRCLUSTER_H
#define TRACKBASE_TRKRCLUSTER_H

#include "TrkrDefs.h"

#include <phool/PHObject.h>

#include <climits>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <memory>

/**
 * @brief Base class for cluster object
 *
 * Virtual base class for cluster object used for all tracking subsystems
 */
class TrkrCluster : public PHObject
{
 public:
  //! dtor
  ~TrkrCluster() override = default;

  // PHObject virtual overloads
  void identify(std::ostream& os = std::cout) const override
  {
    os << "TrkrCluster base class" << std::endl;
  }
  void Reset() override {}
  int isValid() const override { return 0; }

  //! import PHObject CopyFrom, in order to avoid clang warning
  using PHObject::CopyFrom;

  //! copy content from base class
  virtual void CopyFrom(const TrkrCluster&)
  {
  }

  //! copy content from base class
  virtual void CopyFrom(TrkrCluster*)
  {
  }

  //
  // cluster position
  //
  virtual float getLocalX() const { return std::numeric_limits<float>::quiet_NaN(); }
  virtual void setLocalX(float) {}
  virtual float getLocalY() const { return std::numeric_limits<float>::quiet_NaN(); }
  virtual void setLocalY(float) {}

  //
  // cluster info
  //
  virtual void setAdc(uint16_t) {}
  virtual uint16_t getAdc() const { return UINT16_MAX; }
  virtual void setMaxAdc(uint16_t) {}
  virtual uint16_t getMaxAdc() const { return UINT16_MAX; }
  virtual uint8_t getOverlap() const { return std::numeric_limits<uint8_t>::max(); }
  virtual void setOverlap(uint8_t) {}
  virtual uint8_t getEdge() const { return std::numeric_limits<uint8_t>::max(); }
  virtual void setEdge(uint8_t) {}
  virtual void setTime(const float) {}
  virtual float getTime() const { return std::numeric_limits<float>::quiet_NaN(); }
  virtual uint8_t getSize() const { return std::numeric_limits<uint8_t>::max(); }

  //
  // convenience interface
  //
  virtual float getPhiSize() const { return std::numeric_limits<float>::quiet_NaN(); }
  virtual float getZSize() const { return std::numeric_limits<float>::quiet_NaN(); }
  virtual float getPhiError() const { return std::numeric_limits<float>::quiet_NaN(); }
  virtual float getRPhiError() const { return std::numeric_limits<float>::quiet_NaN(); }
  virtual float getZError() const { return std::numeric_limits<float>::quiet_NaN(); }
  virtual uint16_t getCenAdc() const { return UINT16_MAX; }
  virtual float getPadCen() const { return std::numeric_limits<float>::quiet_NaN(); }
  virtual float getTBinCen() const { return std::numeric_limits<float>::quiet_NaN(); }
  virtual float getPadMax() const { return std::numeric_limits<float>::quiet_NaN(); }
  virtual float getTBinMax() const { return std::numeric_limits<float>::quiet_NaN(); }
  virtual uint8_t getSLEdge() const { return std::numeric_limits<uint8_t>::max(); }
  virtual uint8_t getSREdge() const { return std::numeric_limits<uint8_t>::max(); }
  virtual uint8_t getTLEdge() const { return std::numeric_limits<uint8_t>::max(); }
  virtual uint8_t getTREdge() const { return std::numeric_limits<uint8_t>::max(); }
  virtual uint8_t getDLEdge() const { return std::numeric_limits<uint8_t>::max(); }
  virtual uint8_t getDREdge() const { return std::numeric_limits<uint8_t>::max(); }
  virtual uint8_t getHLEdge() const { return std::numeric_limits<uint8_t>::max(); }
  virtual uint8_t getHREdge() const { return std::numeric_limits<uint8_t>::max(); }
  virtual int getSLMix() const { return std::numeric_limits<int>::max(); }
  virtual int getSRMix() const { return std::numeric_limits<int>::max(); }
  virtual int getTLMix() const { return std::numeric_limits<int>::max(); }
  virtual int getTRMix() const { return std::numeric_limits<int>::max(); }
  virtual int getPhiBinLo() const { return std::numeric_limits<int>::max(); }
  virtual int getPhiBinHi() const { return std::numeric_limits<int>::max(); }
  virtual int getTBinLo() const { return std::numeric_limits<int>::max(); }
  virtual int getTBinHi() const { return std::numeric_limits<int>::max(); }
  virtual float getPadPhase() const { return std::numeric_limits<float>::quiet_NaN(); }
  virtual float getTBinPhase() const { return std::numeric_limits<float>::quiet_NaN(); }
  virtual float getRSize() const { return std::numeric_limits<float>::quiet_NaN(); }

  /// Acts functions, for Acts modules use only
  virtual void setActsLocalError(unsigned int /*i*/, unsigned int /*j*/, float /*value*/) {}
  virtual float getActsLocalError(unsigned int /*i*/, unsigned int /*j*/) const { return std::numeric_limits<float>::quiet_NaN(); }
  virtual TrkrDefs::subsurfkey getSubSurfKey() const { return TrkrDefs::SUBSURFKEYMAX; }
  virtual void setSubSurfKey(TrkrDefs::subsurfkey /*id*/) {}

  // Global coordinate functions are deprecated, use local
  // coordinate functions only
  virtual float getX() const { return std::numeric_limits<float>::quiet_NaN(); }
  virtual void setX(float) {}
  virtual float getY() const { return std::numeric_limits<float>::quiet_NaN(); }
  virtual void setY(float) {}
  virtual float getZ() const { return std::numeric_limits<float>::quiet_NaN(); }
  virtual void setZ(float) {}
  virtual float getPosition(int /*coor*/) const { return std::numeric_limits<float>::quiet_NaN(); }
  virtual void setPosition(int /*coor*/, float /*xi*/) {}
  virtual void setGlobal() {}
  virtual void setLocal() {}
  virtual bool isGlobal() const { return true; }
  virtual float getError(unsigned int /*i*/, unsigned int /*j*/) const { return std::numeric_limits<float>::quiet_NaN(); }
  virtual void setError(unsigned int /*i*/, unsigned int /*j*/, float /*value*/) {}
  virtual float getSize(unsigned int /*i*/, unsigned int /*j*/) const { return std::numeric_limits<float>::quiet_NaN(); }
  virtual void setSize(unsigned int /*i*/, unsigned int /*j*/, float /*value*/) {}

 protected:
  TrkrCluster() = default;
  ClassDefOverride(TrkrCluster, 1)
};

#endif  // TRACKBASE_TRKRCLUSTER_H
