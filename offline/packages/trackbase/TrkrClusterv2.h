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
#include <Acts/Surfaces/Surface.hpp>
#include <ActsExamples/EventData/TrkrClusterSourceLink.hpp>
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
  virtual ~TrkrClusterv2() {}
  // PHObject virtual overloads
  virtual void identify(std::ostream& os = std::cout) const;
  virtual void Reset() {}
  virtual int isValid() const;
  virtual PHObject* CloneMe() const { return new TrkrClusterv2(*this); }
  virtual void setClusKey(TrkrDefs::cluskey id) { m_cluskey = id; }
  virtual TrkrDefs::cluskey getClusKey() const { return m_cluskey; }
  //
  // cluster position
  //
  virtual float getX() const { return m_pos[0]; }
  virtual void setX(float x) { m_pos[0] = x; }
  virtual float getY() const { return m_pos[1]; }
  virtual void setY(float y) { m_pos[1] = y; }
  virtual float getZ() const { return m_pos[2]; }
  virtual void setZ(float z) { m_pos[2] = z; }
  virtual float getPosition(int coor) const { return m_pos[coor]; }
  virtual void setPosition(int coor, float xi) { m_pos[coor] = xi; }
  virtual void setGlobal() { m_isGlobal = true; }
  virtual void setLocal() { m_isGlobal = false; }
  virtual bool isGlobal() { return m_isGlobal; }

  virtual float getLocalX() const { return m_local[0]; }
  virtual void setLocalX(float loc0) { m_local[0] = loc0; }
  virtual float getLocalY() const { return m_local[1]; }
  virtual void setLocalY(float loc1) { m_local[1] = loc1; }

  /// Acts functions, for Acts module use only
  void setActsLocalError(unsigned int i, unsigned int j, float value);
  Surface getActsSurface() const { return m_surface; }
  void setActsSurface(Surface surface) { m_surface = surface; }
  SourceLinkPtr getActsSourceLink() const;


  //
  // cluster info
  //
  virtual unsigned int getAdc() const { return m_adc; }
  virtual void setAdc(unsigned int adc) { m_adc = adc; }
  virtual float getSize(unsigned int i, unsigned int j) const;        //< get cluster dimension covar
  virtual void setSize(unsigned int i, unsigned int j, float value);  //< set cluster dimension covar

  virtual float getError(unsigned int i, unsigned int j) const;        //< get cluster error covar
  virtual void setError(unsigned int i, unsigned int j, float value);  //< set cluster error covar
  
  //
  // convenience interface
  //
  virtual float getPhiSize() const;
  virtual float getZSize() const;

  virtual float getRPhiError() const;
  virtual float getPhiError() const;
  virtual float getZError() const;

 protected:

  TrkrDefs::cluskey m_cluskey;  //< unique identifier within container
  float m_pos[3];               //< mean position x,y,z
  bool m_isGlobal;             //< flag for coord sys (true = global)
  unsigned int m_adc;           //< cluster sum adc (D. McGlinchey - Do we need this?)
  float m_size[6];              //< size covariance matrix (packed storage) (+/- cm^2)
  float m_err[6];               //< covariance matrix: rad, arc and z

  float m_local[2];             //< 2D local position [cm]
  Surface m_surface;            //< acts surface that cluster lies on
  float m_actsLocalErr[2][2];   //< 2D local error for Acts [cm]

  ClassDef(TrkrClusterv2, 1)
};

#endif //TRACKBASE_TRKRCLUSTERV2_H
