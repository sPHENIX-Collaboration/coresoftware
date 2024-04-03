/**
 * @file trackbase/LaserClusterv1.h
 * @author Ben Kimelman
 * @date February 2024
 * @brief Version 1 of CMFLashCluster
 */
#ifndef TRACKBASE_LASERCLUSTERV1_H
#define TRACKBASE_LASERCLUSTERV1_H

#include "LaserCluster.h"

#include <iostream>

class PHObject;

/**
 * @brief Version 1 of LaserCluster
 *
 * Note - D. McGlinchey June 2018:
 *   CINT does not like "override", so ignore where CINT
 *   complains. Should be checked with ROOT 6 once
 *   migration occurs.
 */
class LaserClusterv1 : public LaserCluster
{
 public:
  //! ctor
  LaserClusterv1() = default;

  // PHObject virtual overloads
  void identify(std::ostream& os = std::cout) const override;
  void Reset() override {}
  int isValid() const override;
  PHObject* CloneMe() const override { return new LaserClusterv1(*this); }
 
  //! copy content from base class
  void CopyFrom( const LaserCluster& ) override;

  //! copy content from base class
  void CopyFrom( LaserCluster* source ) override
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

  float getLayer() const override { return m_posHardware[0]; }
  void setLayer(float layer) override { m_posHardware[0] = layer; }
  float getIPhi() const override { return m_posHardware[1]; }
  void setIPhi(float iphi) override { m_posHardware[1] = iphi; }
  float getIT() const override { return m_posHardware[2]; }
  void setIT(float it) override { m_posHardware[2] = it; }

  unsigned int getNhits() const override {return m_hitVec.size();}

  //
  // cluster info
  //
  unsigned int getAdc() const override { return m_adc; }
  void setAdc(unsigned int adc) override { m_adc = adc; }

  void addHit() override {m_hitVec.push_back({0.0,0.0,0.0,0.0}); m_hitVecHardware.push_back({0.0,0.0,0.0}); }

  void setHitLayer(int i, float layer) override {m_hitVecHardware[i][0] = layer; }
  float getHitLayer(int i) const override { return m_hitVecHardware[i][0]; }

  void setHitIPhi(int i, float iphi) override {m_hitVecHardware[i][1] = iphi; }
  float getHitIPhi(int i) const override { return m_hitVecHardware[i][1]; }

  void setHitIT(int i, float it) override {m_hitVecHardware[i][2] = it; }
  float getHitIT(int i) const override { return m_hitVecHardware[i][2]; }

  void setHitX(int i, float x) override {m_hitVec[i][0] = x; }
  float getHitX(int i) const override { return m_hitVec[i][0]; }

  void setHitY(int i, float y) override {m_hitVec[i][1] = y; }
  float getHitY(int i) const override { return m_hitVec[i][1]; }

  void setHitZ(int i, float z) override {m_hitVec[i][2] = z; }
  float getHitZ(int i) const override { return m_hitVec[i][2]; }

  void setHitAdc(int i, float adc) override {m_hitVec[i][3] = adc; }
  float getHitAdc(int i) const override { return m_hitVec[i][3]; }


 protected:

  /// mean cluster position
  float m_pos[3] = {NAN, NAN, NAN};          
  float m_posHardware[3] = {NAN, NAN, NAN};

  std::vector<std::vector<float>> m_hitVec;
  std::vector<std::vector<float>> m_hitVecHardware;

  /// cluster sum adc
  unsigned int m_adc = 0xFFFFFFFF;

  /// number of TPC clusters used to create this central mebrane cluster
  unsigned int m_nhits = UINT_MAX;

  ClassDefOverride(LaserClusterv1, 1)
};

#endif //TRACKBASE_CMFLASHCLUSTERV1_H
