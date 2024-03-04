/**
 * @file trackbase/CMFlashClusterv3.h
 * @author Ben Kimelman
 * @date March 2023
 * @brief Version 3 of CMFlashCluster
 */
#ifndef TRACKBASE_CMFLASHCLUSTERV3_H
#define TRACKBASE_CMFLASHCLUSTERV3_H

#include "CMFlashCluster.h"

#include <iostream>

class PHObject;

/**
 * @brief Version 3 of CMFlashCluster
 *
 *Adding variable to keep track of clusters
 *put into metaclusters
 *
 */
class CMFlashClusterv3 : public CMFlashCluster
{
 public:
  //! ctor
  CMFlashClusterv3() = default;

  // PHObject virtual overloads
  void identify(std::ostream& os = std::cout) const override;
  void Reset() override {}
  int isValid() const override;
  PHObject* CloneMe() const override { return new CMFlashClusterv3(*this); }

  //! copy content from base class
  void CopyFrom(const CMFlashCluster&) override;

  //! copy content from base class
  void CopyFrom(CMFlashCluster* source) override
  {
    CopyFrom(*source);
  }

  //
  // cluster position
  //
  float getX() const override { return m_pos[0]; }
  void setX(float x) override { m_pos[0] = x; }
  float getY() const override { return m_pos[1]; }
  void setY(float y) override { m_pos[1] = y; }
  float getZ() const override { return m_pos[2]; }
  void setZ(float z) override { m_pos[2] = z; }

  float getX1() const override { return m_pos1[0]; }
  void setX1(float x) override { m_pos1[0] = x; }
  float getY1() const override { return m_pos1[1]; }
  void setY1(float y) override { m_pos1[1] = y; }
  float getZ1() const override { return m_pos1[2]; }
  void setZ1(float z) override { m_pos1[2] = z; }

  float getX2() const override { return m_pos2[0]; }
  void setX2(float x) override { m_pos2[0] = x; }
  float getY2() const override { return m_pos2[1]; }
  void setY2(float y) override { m_pos2[1] = y; }
  float getZ2() const override { return m_pos2[2]; }
  void setZ2(float z) override { m_pos2[2] = z; }

  unsigned int getLayer1() const override { return m_layer1; }
  void setLayer1(unsigned int layer) override { m_layer1 = layer; }
  unsigned int getLayer2() const override { return m_layer2; }
  void setLayer2(unsigned int layer) override { m_layer2 = layer; }

  unsigned int getNclusters() const override { return m_nclusters; }
  void setNclusters(unsigned int n) override { m_nclusters = n; }
  bool getIsRGap() const override { return m_isRGap; }
  void setIsRGap(bool isRGap) override { m_isRGap = isRGap; }
  bool getIsPhiGap() const override { return m_isPhiGap; }
  void setIsPhiGap(bool isPhiGap) override { m_isPhiGap = isPhiGap; }

  //
  // cluster info
  //
  unsigned int getAdc() const override { return m_adc; }
  void setAdc(unsigned int adc) override { m_adc = adc; }

  unsigned int getAdc1() const override { return m_adc1; }
  void setAdc1(unsigned int adc) override { m_adc1 = adc; }

  unsigned int getAdc2() const override { return m_adc2; }
  void setAdc2(unsigned int adc) override { m_adc2 = adc; }

 protected:
  /// mean cluster position
  float m_pos[3] = {NAN, NAN, NAN};

  float m_pos1[3] = {NAN, NAN, NAN};
  float m_pos2[3] = {NAN, NAN, NAN};

  /// cluster sum adc
  unsigned int m_adc = 0xFFFFFFFF;
  unsigned int m_adc1 = 0xFFFFFFFF;
  unsigned int m_adc2 = 0xFFFFFFFF;

  unsigned int m_layer1 = UINT_MAX;
  unsigned int m_layer2 = UINT_MAX;

  /// number of TPC clusters used to create this central mebrane cluster
  unsigned int m_nclusters = UINT_MAX;

  /// bools to identify if meta-cluster is across sector/module gaps
  bool m_isRGap = false;
  bool m_isPhiGap = false;

  ClassDefOverride(CMFlashClusterv3, 1)
};

#endif  // TRACKBASE_CMFLASHCLUSTERV3_H
