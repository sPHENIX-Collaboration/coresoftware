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
  unsigned int getNhits() const override {return m_nhits;}
  void setNhits(unsigned int n) override { m_nhits = n;}

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
  unsigned int m_nhits = UINT_MAX;

  ClassDefOverride(LaserClusterv1, 1)
};

#endif //TRACKBASE_CMFLASHCLUSTERV1_H
