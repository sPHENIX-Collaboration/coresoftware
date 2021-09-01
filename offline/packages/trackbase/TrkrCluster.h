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
  ~TrkrCluster() override {}
  // PHObject virtual overloads
  void identify(std::ostream& os = std::cout) const override
  {
    os << "TrkrCluster base class" << std::endl;
  }
  void Reset() override {}
  int isValid() const override { return 0; }
  //
  // cluster id
  //
  virtual TrkrDefs::cluskey getClusKey() const { return TrkrDefs::CLUSKEYMAX; }
  virtual void setClusKey(TrkrDefs::cluskey) {}
 
  //
  // cluster position
  //
  virtual float getX() const { return NAN; }
  virtual void setX(float) {}
  virtual float getY() const { return NAN; }
  virtual void setY(float) {}
  virtual float getZ() const { return NAN; }
  virtual void setZ(float) {}
  virtual float getPosition(int) const { return NAN; }
  virtual void setPosition(int, float) {}
  virtual void setGlobal() {}
  virtual void setLocal() {}
  virtual bool isGlobal() { return true; }
  virtual float getLocalX() const { return NAN; }
  virtual void setLocalX(float) {}
  virtual float getLocalY() const { return NAN; }
  virtual void setLocalY(float) {}
  //
  // cluster info
  //
  virtual void setAdc(unsigned int) {}
  virtual unsigned int getAdc() const { return UINT_MAX; }
  virtual float getSize(unsigned int, unsigned int) const { return NAN; }
  virtual void setSize(unsigned int, unsigned int, float) {}
  virtual float getError(unsigned int, unsigned int) const { return NAN; }
  virtual void setError(unsigned int, unsigned int, float) {}
  //
  // convenience interface
  //
  virtual float getPhiSize() const { return NAN; }
  virtual float getZSize() const { return NAN; }
  virtual float getPhiError() const { return NAN; }
  virtual float getRPhiError() const { return NAN; }
  virtual float getZError() const { return NAN; }

  /// Acts functions, for Acts modules use only
  virtual void setActsLocalError(unsigned int, unsigned int, float){}
  virtual float getActsLocalError(unsigned int, unsigned int) const { return NAN; }
  virtual TrkrDefs::subsurfkey getSubSurfKey() const { return TrkrDefs::SUBSURFKEYMAX; }
  virtual void setSubSurfKey(TrkrDefs::subsurfkey) {}

 protected:
  TrkrCluster() = default;
  ClassDefOverride(TrkrCluster, 1)
};

#endif //TRACKBASE_TRKRCLUSTER_H
