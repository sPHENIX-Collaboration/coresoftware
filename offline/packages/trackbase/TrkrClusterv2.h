/**
 * @file trackbase/TrkrClusterv2.h
 * @author J. Osborn
 * @date March 2021
 * @brief Version 2 of TrkrCluster
 */
#ifndef TRACKBASE_TRKRCLUSTERV2_H
#define TRACKBASE_TRKRCLUSTERV2_H

#include "TrkrCluster.h"
#include "TrkrDefs.h"
#include <iostream>

class PHObject;

/**
 * @brief Version 2 of TrkrCluster
 *
 * This version of TrkrCluster contains Acts source link objects
 * as member variables, to join the Svtx and Acts worlds
 */
class TrkrClusterv2 : public TrkrCluster
{
 public:
 
  //! ctor
  TrkrClusterv2();

  //!dtor
  ~TrkrClusterv2() override = default;
  
  // PHObject virtual overloads
  void identify(std::ostream& os = std::cout) const override;
  void Reset() override {}
  int isValid() const override;
  PHObject* CloneMe() const override { return new TrkrClusterv2(*this); }
  
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

  float getLocalX() const override { return m_local[0]; }
  void setLocalX(float loc0) override { m_local[0] = loc0; }
  float getLocalY() const override { return m_local[1]; }
  void setLocalY(float loc1) override { m_local[1] = loc1; }

  /// Acts functions, for Acts module use only
  void setActsLocalError(unsigned int i, unsigned int j, float value) override;
  float getActsLocalError(unsigned int i, unsigned int j) const override { return m_actsLocalErr[i][j]; }
  TrkrDefs::subsurfkey getSubSurfKey() const override { return m_subsurfkey; }
  void setSubSurfKey(TrkrDefs::subsurfkey id) override { m_subsurfkey = id; }

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
  TrkrDefs::subsurfkey m_subsurfkey; //< unique identifier for hitsetkey-surface maps
  float m_pos[3];               //< mean position x,y,z
  bool m_isGlobal;             //< flag for coord sys (true = global)
  unsigned int m_adc;           //< cluster sum adc (D. McGlinchey - Do we need this?)
  float m_size[6];              //< size covariance matrix (packed storage) (+/- cm^2)
  float m_err[6];               //< covariance matrix: rad, arc and z

  float m_local[2];             //< 2D local position [cm]
  float m_actsLocalErr[2][2];   //< 2D local error for Acts [cm]

  ClassDefOverride(TrkrClusterv2, 2)
};

#endif //TRACKBASE_TRKRCLUSTERV2_H
