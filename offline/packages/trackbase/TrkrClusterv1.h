#ifndef __TrkrClusterv1_H__
#define __TrkrClusterv1_H__

#include "TrkrCluster.h"
#include "TrkrDefs.h"

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
  void setClusKey(TrkrDefs::cluskey id) { m_cluskey = id; }
  TrkrDefs::cluskey getClusKey() const { return m_cluskey; }
  //
  // cluster position
  //
  float getX() const { return m_pos[0]; }
  void setX(float x) { m_pos[0] = x; }
  float getY() const { return m_pos[1]; }
  void setY(float y) { m_pos[1] = y; }
  float getZ() const { return m_pos[2]; }
  void setZ(float z) { m_pos[2] = z; }
  float getPosition(int coor) const { return m_pos[coor]; }
  void setPosition(int coor, float xi) { m_pos[coor] = xi; }
  void setGlobal() { m_isGlobal = true; }
  void setLocal() { m_isGlobal = false; }
  bool isGlobal() { return m_isGlobal; }
  //
  // cluster info
  //
  unsigned int getAdc() const { return m_adc; }
  void setAdc(unsigned int adc) { m_adc = adc; }
  float getSize(unsigned int i, unsigned int j) const;        //< get cluster dimension covar
  void setSize(unsigned int i, unsigned int j, float value);  //< set cluster dimension covar

  float getError(unsigned int i, unsigned int j) const;        //< get cluster error covar
  void setError(unsigned int i, unsigned int j, float value);  //< set cluster error covar

  //
  // convenience interface
  //
  float getPhiSize() const;
  float getZSize() const;

  float getRPhiError() const;
  float getPhiError() const;
  float getZError() const;

 private:
  unsigned int covarIndex(unsigned int i, unsigned int j) const;

  TrkrDefs::cluskey m_cluskey;  //< unique identifier within container
  float m_pos[3];               //< mean position x,y,z
  bool m_isGlobal;             //< flag for coord sys (true = global)
  unsigned int m_adc;           //< cluster sum adc (D. McGlinchey - Do we need this?)
  float m_size[6];              //< size covariance matrix (packed storage) (+/- cm^2)
  float m_err[6];               //< covariance matrix: rad, arc and z

  ClassDef(TrkrClusterv1, 1);
};

#endif
