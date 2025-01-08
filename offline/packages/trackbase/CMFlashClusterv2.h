/**
 * @file trackbase/CMFlashClusterv2.h
 * @author Ben Kimelman
 * @date March 2023
 * @brief Version 2 of CMFlashCluster
 */
#ifndef TRACKBASE_CMFLASHCLUSTERV2_H
#define TRACKBASE_CMFLASHCLUSTERV2_H

#include "CMFlashCluster.h"

#include <iostream>

class PHObject;

/**
 * @brief Version 2 of CMFlashCluster
 *
 * Adding bool variables to identify if
 * meta-cluster contains sub-clusters from
 * both sides of sector (phi) or module (R)
 * gaps
 *
 */
class CMFlashClusterv2 : public CMFlashCluster
{
 public:
  //! ctor
  CMFlashClusterv2() = default;

  // PHObject virtual overloads
  void identify(std::ostream& os = std::cout) const override;
  void Reset() override {}
  int isValid() const override;
  PHObject* CloneMe() const override { return new CMFlashClusterv2(*this); }

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

 protected:
  /// mean cluster position
  float m_pos[3] = {NAN, NAN, NAN};

  /// cluster sum adc
  unsigned int m_adc = 0xFFFFFFFF;

  /// number of TPC clusters used to create this central mebrane cluster
  unsigned int m_nclusters = UINT_MAX;

  /// bools to identify if meta-cluster is across sector/module gaps
  bool m_isRGap = false;
  bool m_isPhiGap = false;

  ClassDefOverride(CMFlashClusterv2, 1)
};

#endif  // TRACKBASE_CMFLASHCLUSTERV2_H
