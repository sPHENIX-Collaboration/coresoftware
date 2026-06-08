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
  virtual float getLocalX() const { return NAN; }
  virtual void setLocalX(float) {}
  virtual float getLocalY() const { return NAN; }
  virtual void setLocalY(float) {}

  //
  // cluster info
  //
  virtual void setAdc(unsigned int) {}
  virtual unsigned int getAdc() const { return UINT_MAX; }
  virtual void setMaxAdc(uint16_t) {}
  virtual unsigned int getMaxAdc() const { return UINT_MAX; }
  virtual char getOverlap() const { return std::numeric_limits<char>::max(); }
  virtual void setOverlap(char) {}
  virtual char getEdge() const { return std::numeric_limits<char>::max(); }
  virtual void setEdge(char) {}
  virtual void setTime(const float) {}
  virtual float getTime() const { return NAN; }
  virtual char getSize() const { return std::numeric_limits<char>::max(); }

  //
  // convenience interface
  //
  virtual float getPhiSize() const { return NAN; }
  virtual float getZSize() const { return NAN; }
  virtual float getPhiError() const { return NAN; }
  virtual float getRPhiError() const { return NAN; }
  virtual float getZError() const { return NAN; }
  virtual unsigned int getCenAdc() const { return UINT_MAX; }
  virtual float getPadCen() const { return NAN; }
  virtual float getTBinCen() const { return NAN; }
  virtual float getPadMax() const { return NAN; }
  virtual float getTBinMax() const { return NAN; }
  virtual char getSLEdge() const { return std::numeric_limits<char>::max(); }
  virtual char getSREdge() const { return std::numeric_limits<char>::max(); }
  virtual char getTLEdge() const { return std::numeric_limits<char>::max(); }
  virtual char getTREdge() const { return std::numeric_limits<char>::max(); }
  virtual char getDLEdge() const { return std::numeric_limits<char>::max(); }
  virtual char getDREdge() const { return std::numeric_limits<char>::max(); }
  virtual char getHLEdge() const { return std::numeric_limits<char>::max(); }
  virtual char getHREdge() const { return std::numeric_limits<char>::max(); }
  virtual int getSLMix() const { return std::numeric_limits<int>::max(); }
  virtual int getSRMix() const { return std::numeric_limits<int>::max(); }
  virtual int getTLMix() const { return std::numeric_limits<int>::max(); }
  virtual int getTRMix() const { return std::numeric_limits<int>::max(); }
  virtual float getPhiBinLo() const { return NAN; }
  virtual float getPhiBinHi() const { return NAN; }
  virtual float getTBinLo() const { return NAN; }
  virtual float getTBinHi() const { return NAN; }
  virtual float getPadPhase() const { return NAN; }
  virtual float getTBinPhase() const { return NAN; }
  virtual float getRSize() const { return NAN; }

  /// Acts functions, for Acts modules use only
  virtual void setActsLocalError(unsigned int /*i*/, unsigned int /*j*/, float /*value*/) {}
  virtual float getActsLocalError(unsigned int /*i*/, unsigned int /*j*/) const { return NAN; }
  virtual TrkrDefs::subsurfkey getSubSurfKey() const { return TrkrDefs::SUBSURFKEYMAX; }
  virtual void setSubSurfKey(TrkrDefs::subsurfkey /*id*/) {}

  // Global coordinate functions are deprecated, use local
  // coordinate functions only
  virtual float getX() const { return NAN; }
  virtual void setX(float) {}
  virtual float getY() const { return NAN; }
  virtual void setY(float) {}
  virtual float getZ() const { return NAN; }
  virtual void setZ(float) {}
  virtual float getPosition(int /*coor*/) const { return NAN; }
  virtual void setPosition(int /*coor*/, float /*xi*/) {}
  virtual void setGlobal() {}
  virtual void setLocal() {}
  virtual bool isGlobal() const { return true; }
  virtual float getError(unsigned int /*i*/, unsigned int /*j*/) const { return NAN; }
  virtual void setError(unsigned int /*i*/, unsigned int /*j*/, float /*value*/) {}
  virtual float getSize(unsigned int /*i*/, unsigned int /*j*/) const { return NAN; }
  virtual void setSize(unsigned int /*i*/, unsigned int /*j*/, float /*value*/) {}

 protected:
  TrkrCluster() = default;
  ClassDefOverride(TrkrCluster, 1)
};

#endif  // TRACKBASE_TRKRCLUSTER_H
