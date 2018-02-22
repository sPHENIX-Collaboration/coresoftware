#ifndef __TrkrClusterv1_H__
#define __TrkrClusterv1_H__

#include "TrkrCluster.h"
#include "TrkrDefUtil.h"

#include <limits.h>
#include <cmath>
#include <iostream>
#include <set>

class TrkrClusterv1 : public TrkrCluster
{
 public:
  //! ctor
  TrkrClusterv1();

  //!dtor
  virtual ~TrkrClusterv1() {}
  // PHObject virtual overloads
  void identify(std::ostream& os = std::cout) const;
  void Reset() {}
  int isValid() const;
  TrkrCluster* clone() const { return new TrkrClusterv1(*this); }
  void SetClusKey(TrkrDefs::cluskey id) { cluskey_ = id; }
  TrkrDefs::cluskey GetClusKey() const { return cluskey_; }
  //
  // cluster position
  //
  float GetX() const { return pos_[0]; }
  void SetX(float x) { pos_[0] = x; }
  float GetY() const { return pos_[1]; }
  void SetY(float y) { pos_[1] = y; }
  float GetZ() const { return pos_[2]; }
  void SetZ(float z) { pos_[2] = z; }
  float GetPosition(int coor) const { return pos_[coor]; }
  void SetPosition(int coor, float xi) { pos_[coor] = xi; }
  void SetGlobal() { is_global_ = true; }
  void SetLocal() { is_global_ = false; }
  bool IsGlobal() { return is_global_; }
  //
  // cluster info
  //
  unsigned int GetAdc() const { return adc_; }
  void SetAdc(unsigned int adc) { adc_ = adc; }
  float GetSize(unsigned int i, unsigned int j) const;        //< get cluster dimension covar
  void SetSize(unsigned int i, unsigned int j, float value);  //< set cluster dimension covar

  float GetError(unsigned int i, unsigned int j) const;        //< get cluster error covar
  void SetError(unsigned int i, unsigned int j, float value);  //< set cluster error covar

  //
  // convenience interface
  //
  float GetPhiSize() const;
  float GetZSize() const;

  float GetRPhiError() const;
  float GetPhiError() const;
  float GetZError() const;

 private:
  unsigned int CovarIndex(unsigned int i, unsigned int j) const;

  TrkrDefs::cluskey cluskey_;  //< unique identifier within container
  float pos_[3];               //< mean position x,y,z
  bool is_global_;             //< flag for coord sys (true = global)
  unsigned int adc_;           //< cluster sum adc (D. McGlinchey - Do we need this?)
  float size_[6];              //< size covariance matrix (packed storage) (+/- cm^2)
  float err_[6];               //< covariance matrix: rad, arc and z

  ClassDef(TrkrClusterv1, 1);
};

#endif
