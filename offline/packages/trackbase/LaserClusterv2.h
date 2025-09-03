/**
 * @file trackbase/LaserClusterv2.h
 * @author Ben Kimelman
 * @date July 2025
 * @brief Version 2 of CMFLashCluster
 */
#ifndef TRACKBASE_LASERCLUSTERV2_H
#define TRACKBASE_LASERCLUSTERV2_H

#include "LaserCluster.h"

#include <iostream>
#include <vector>

class PHObject;

/**
 * @brief Version 2 of LaserCluster
 *
 * Note - D. McGlinchey June 2018:
 *   CINT does not like "override", so ignore where CINT
 *   complains. Should be checked with ROOT 6 once
 *   migration occurs.
 */
class LaserClusterv2 : public LaserCluster
{
 public:
  //! ctor
  LaserClusterv2() = default;

  // PHObject virtual overloads
  void identify(std::ostream& os = std::cout) const override;
  void Reset() override {}
  int isValid() const override;
  PHObject* CloneMe() const override { return new LaserClusterv2(*this); }
 
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

  bool getFitMode() const override { return m_fitMode; }
  void setFitMode(bool fitMode) override { m_fitMode = fitMode; }

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

  void setNLayers(unsigned int nLayers) override { m_nLayers = nLayers; }
  unsigned int getNLayers() const override { return m_nLayers; }

  void setNIPhi(unsigned int nIPhi) override { m_nIPhi = nIPhi; }
  unsigned int getNIPhi() const override { return m_nIPhi; }

  void setNIT(unsigned int nIT) override { m_nIT = nIT; }
  unsigned int getNIT() const override { return m_nIT; }

  void setSDLayer(float SDLayer) override { m_SDLayer = SDLayer; }
  float getSDLayer() const override { return m_SDLayer; }

  void setSDIPhi(float SDIPhi) override { m_SDIPhi = SDIPhi; }
  float getSDIPhi() const override { return m_SDIPhi; }

  void setSDIT(float SDIT) override { m_SDIT = SDIT; }
  float getSDIT() const override { return m_SDIT; }

  void setSDWeightedLayer(float SDLayer) override { m_SDWeightedLayer = SDLayer; }
  float getSDWeightedLayer() const override { return m_SDWeightedLayer; }

  void setSDWeightedIPhi(float SDIPhi) override { m_SDWeightedIPhi = SDIPhi; }
  float getSDWeightedIPhi() const override { return m_SDWeightedIPhi; }

  void setSDWeightedIT(float SDIT) override { m_SDWeightedIT = SDIT; }
  float getSDWeightedIT() const override { return m_SDWeightedIT; }

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
  float m_pos[3] = {std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN()};
  float m_posHardware[3] = {std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN()};

  bool m_fitMode{false};

  std::vector<std::vector<float>> m_hitVec;
  std::vector<std::vector<float>> m_hitVecHardware;

  /// cluster sum adc
  unsigned int m_adc = std::numeric_limits<unsigned int>::max();
  
  /// number of TPC clusters used to create this central mebrane cluster
  unsigned int m_nhits = std::numeric_limits<unsigned int>::max();
  unsigned int m_nLayers = std::numeric_limits<unsigned int>::max();
  unsigned int m_nIPhi = std::numeric_limits<unsigned int>::max();
  unsigned int m_nIT = std::numeric_limits<unsigned int>::max();
  float m_SDLayer = std::numeric_limits<float>::quiet_NaN();
  float m_SDIPhi = std::numeric_limits<float>::quiet_NaN();
  float m_SDIT = std::numeric_limits<float>::quiet_NaN();
  float m_SDWeightedLayer = std::numeric_limits<float>::quiet_NaN();
  float m_SDWeightedIPhi = std::numeric_limits<float>::quiet_NaN();
  float m_SDWeightedIT = std::numeric_limits<float>::quiet_NaN();

  ClassDefOverride(LaserClusterv2, 1)
};

#endif //TRACKBASE_LASERCLUSTERV2_H
