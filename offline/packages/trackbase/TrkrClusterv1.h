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
  ~TrkrClusterv1() override = default;
  // PHObject virtual overloads
  void identify(std::ostream& os = std::cout) const override;
  void Reset() override {}
  int isValid() const override;
  PHObject* CloneMe() const override { return new TrkrClusterv1(*this); }
 
  //! copy content from base class
  void CopyFrom( const TrkrCluster& ) override;

  //! copy content from base class
  void CopyFrom( TrkrCluster* source ) override
  { CopyFrom( *source ); }

  //
  // cluster position
  //
  float getX() const override { return m_pos[0]; }
  void setX(float x) override { m_pos[0] = x; }
  float getY() const override { return m_pos[1]; }
  void setY(float y) override { m_pos[1] = y; }
  float getZ() const override { return m_pos[2]; }
  void setZ(float z) override { m_pos[2] = z; }
  float getPosition(int coor) const override { return m_pos[coor]; }
  void setPosition(int coor, float xi) override { m_pos[coor] = xi; }
  void setGlobal() override { m_isGlobal = true; }
  void setLocal() override { m_isGlobal = false; }
  bool isGlobal() const override { return m_isGlobal; }
  //
  // cluster info
  //
  unsigned int getAdc() const override { return m_adc; }
  void setAdc(unsigned int adc) override { m_adc = adc; }
  float getSize(unsigned int i, unsigned int j) const override;        //< get cluster dimension covar
  void setSize(unsigned int i, unsigned int j, float value) override;  //< set cluster dimension covar

  float getError(unsigned int i, unsigned int j) const override;        //< get cluster error covar
  void setError(unsigned int i, unsigned int j, float value) override;  //< set cluster error covar

  //
  // convenience interface
  //
  float getPhiSize() const override;
  float getZSize() const override;

  float getRPhiError() const override;
  float getPhiError() const override;
  float getZError() const override;

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
