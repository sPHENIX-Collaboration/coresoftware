#ifndef __TrkrCluster_H__
#define __TrkrCluster_H__

#include <phool/PHObject.h>
#include "TrkrDefUtil.h"

#include <limits.h>
#include <cmath>
#include <iostream>
#include <set>

class TrkrCluster : public PHObject
{
 public:
  //! dtor
  virtual ~TrkrCluster() {}
  // PHObject virtual overloads
  virtual void identify(std::ostream& os = std::cout) const
  {
    os << "TrkrCluster base class" << std::endl;
  }
  virtual void Reset() {}
  virtual int isValid() const { return 0; }
  //
  // cluster id
  //
  virtual TrkrDefs::cluskey GetClusKey() const { return TrkrDefs::CLUSKEYMAX; }
  virtual void SetClusKey(TrkrDefs::cluskey id) {}
  //
  // cluster position
  //
  virtual float GetX() const { return NAN; }
  virtual void SetX(float x) {}
  virtual float GetY() const { return NAN; }
  virtual void SetY(float y) {}
  virtual float GetZ() const { return NAN; }
  virtual void SetZ(float z) {}
  virtual float GetPosition(int coor) const { return NAN; }
  virtual void SetPosition(int coor, float xi) {}
  virtual void SetGlobal() {}
  virtual void SetLocal() {}
  virtual bool IsGlobal() { return true; }
  //
  // cluster info
  //
  virtual void SetAdc(unsigned int adc) {}
  virtual unsigned int GetAdc() const { return UINT_MAX; }
  virtual float GetSize(unsigned int i, unsigned int j) const { return NAN; }
  virtual void SetSize(unsigned int i, unsigned int j, float value) {}
  virtual float GetError(unsigned int i, unsigned int j) const { return NAN; }
  virtual void SetError(unsigned int i, unsigned int j, float value) {}
  //
  // convenience interface
  //
  virtual float GetPhiSize() const { return NAN; }
  virtual float GetZSize() const { return NAN; }
  virtual float GetPhiError() const { return NAN; }
  virtual float GetRPhiError() const { return NAN; }
  virtual float GetZError() const { return NAN; }
 protected:
  TrkrCluster() {}
  ClassDef(TrkrCluster, 1);
};

#endif
