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
  virtual void setClusKey(TrkrDefs::cluskey id) {}
 
  //
  // cluster position
  //
  virtual float getX() const { return NAN; }
  virtual void setX(float x) {}
  virtual float getY() const { return NAN; }
  virtual void setY(float y) {}
  virtual float getZ() const { return NAN; }
  virtual void setZ(float z) {}
  virtual float getPosition(int coor) const { return NAN; }
  virtual void setPosition(int coor, float xi) {}
  virtual void setGlobal() {}
  virtual void setLocal() {}
  virtual bool isGlobal() { return true; }
  virtual float getLocalX() const { return NAN; }
  virtual void setLocalX(float x) {}
  virtual float getLocalY() const { return NAN; }
  virtual void setLocalY(float y) {}
  //
  // cluster info
  //
  virtual void setAdc(unsigned int adc) {}
  virtual unsigned int getAdc() const { return UINT_MAX; }
  virtual float getSize(unsigned int i, unsigned int j) const { return NAN; }
  virtual void setSize(unsigned int i, unsigned int j, float value) {}
  virtual float getError(unsigned int i, unsigned int j) const { return NAN; }
  virtual void setError(unsigned int i, unsigned int j, float value) {}
  //
  // convenience interface
  //
  virtual float getPhiSize() const { return NAN; }
  virtual float getZSize() const { return NAN; }
  virtual float getPhiError() const { return NAN; }
  virtual float getRPhiError() const { return NAN; }
  virtual float getZError() const { return NAN; }

  /// Acts functions, for Acts modules use only
  virtual void setActsLocalError(unsigned int i, unsigned int j, float value){}
  virtual float getActsLocalError(unsigned int i, unsigned int j) const { return NAN; }
  virtual TrkrDefs::subsurfkey getSubSurfKey() const { return TrkrDefs::SUBSURFKEYMAX; }
  virtual void setSubSurfKey(TrkrDefs::subsurfkey id) {}

 protected:
  TrkrCluster() = default;
  ClassDefOverride(TrkrCluster, 1)
};

#endif //TRACKBASE_TRKRCLUSTER_H
