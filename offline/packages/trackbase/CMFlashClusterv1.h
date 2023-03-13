/**
 * @file trackbase/CMFlashClusterv1.h
 * @author Tony Frawley
 * @date January 2022
 * @brief Version 1 of CMFLashCluster
 */
#ifndef TRACKBASE_CMFLASHCLUSTERV1_H
#define TRACKBASE_CMFLASHCLUSTERV1_H

#include "CMFlashCluster.h"

#include <iostream>

class PHObject;

/**
 * @brief Version 1 of CMFlashCluster
 *
 * Note - D. McGlinchey June 2018:
 *   CINT does not like "override", so ignore where CINT
 *   complains. Should be checked with ROOT 6 once
 *   migration occurs.
 */
class CMFlashClusterv1 : public CMFlashCluster
{
 public:
  //! ctor
  CMFlashClusterv1() = default;

  // PHObject virtual overloads
  void identify(std::ostream& os = std::cout) const override;
  void Reset() override {}
  int isValid() const override;
  PHObject* CloneMe() const override { return new CMFlashClusterv1(*this); }
 
  //! copy content from base class
  void CopyFrom( const CMFlashCluster& ) override;

  //! copy content from base class
  void CopyFrom( CMFlashCluster* source ) override
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
  unsigned int getNclusters() const override {return m_nclusters;}
  void setNclusters(unsigned int n) override { m_nclusters = n;}

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

  ClassDefOverride(CMFlashClusterv1, 1)
};

#endif //TRACKBASE_CMFLASHCLUSTERV1_H
