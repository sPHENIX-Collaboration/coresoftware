/**
 * @file trackbase/TrkrClusterv1.h
 * @author D. McGlinchey
 * @date June 2018
 * @brief Version 1 of TrkrCluster
 */
#ifndef TRACKBASE_TRKRCLUSTERV1_H
#define TRACKBASE_TRKRCLUSTERV1_H

#include "TrkrCluster.h"
#include "TrkrDefs.h"

#include <iostream>

class PHObject;

/**
 * @brief Version 1 of TrkrCluster
 *
 * Note - D. McGlinchey June 2018:
 *   CINT does not like "override", so ignore where CINT
 *   complains. Should be checked with ROOT 6 once
 *   migration occurs.
 */
class TrkrClusterv1 : public TrkrCluster
{
 public:
  //! ctor
  TrkrClusterv1();

  //!dtor
  virtual ~TrkrClusterv1() {}
  // PHObject virtual overloads
  virtual void identify(std::ostream& os = std::cout) const override;
  virtual void Reset() override {}
  virtual int isValid() const override;
  virtual PHObject* CloneMe() const override { return new TrkrClusterv1(*this); }

  virtual void setClusKey(TrkrDefs::cluskey id) override { m_cluskey = id; }
  virtual TrkrDefs::cluskey getClusKey() const override { return m_cluskey; }
  //
  // cluster position
  //
  virtual float getX() const override { return m_pos[0]; }
  virtual void setX(float x) override { m_pos[0] = x; }
  virtual float getY() const override { return m_pos[1]; }
  virtual void setY(float y) override { m_pos[1] = y; }
  virtual float getZ() const override { return m_pos[2]; }
  virtual void setZ(float z) override { m_pos[2] = z; }
  virtual float getPosition(int coor) const override { return m_pos[coor]; }
  virtual void setPosition(int coor, float xi) override { m_pos[coor] = xi; }
  virtual void setGlobal() override { m_isGlobal = true; }
  virtual void setLocal() override { m_isGlobal = false; }
  virtual bool isGlobal() override { return m_isGlobal; }
  //
  // cluster info
  //
  virtual unsigned int getAdc() const override { return m_adc; }
  virtual void setAdc(unsigned int adc) override { m_adc = adc; }
  virtual float getSize(unsigned int i, unsigned int j) const override;        //< get cluster dimension covar
  virtual void setSize(unsigned int i, unsigned int j, float value) override;  //< set cluster dimension covar

  virtual float getError(unsigned int i, unsigned int j) const override;        //< get cluster error covar
  virtual void setError(unsigned int i, unsigned int j, float value) override;  //< set cluster error covar

  //
  // convenience interface
  //
  virtual float getPhiSize() const override;
  virtual float getZSize() const override;

  virtual float getRPhiError() const override;
  virtual float getPhiError() const override;
  virtual float getZError() const override;

 protected:

  TrkrDefs::cluskey m_cluskey;  //< unique identifier within container
  float m_pos[3];               //< mean position x,y,z
  bool m_isGlobal;             //< flag for coord sys (true = global)
  unsigned int m_adc;           //< cluster sum adc (D. McGlinchey - Do we need this?)
  float m_size[6];              //< size covariance matrix (packed storage) (+/- cm^2)
  float m_err[6];               //< covariance matrix: rad, arc and z

  ClassDefOverride(TrkrClusterv1, 1)
};

#endif //TRACKBASE_TRKRCLUSTERV1_H
