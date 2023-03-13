/**
 * @file trackbase/CMFlashDifferencev1.h
 * @author Tony Frawley
 * @date January 2022
 * @brief Version 1 of CMFLashDifference
 */
#ifndef TRACKBASE_CMFLASHDIFFERENCEV1_H
#define TRACKBASE_CMFLASHDIFFERENCEV1_H

#include "CMFlashDifference.h"

#include <iostream>

class PHObject;

/**
 * @brief Version 1 of CMFlashDifference
 *
 */
class CMFlashDifferencev1 : public CMFlashDifference
{
 public:
  //! ctor
  CMFlashDifferencev1() = default;

  // PHObject virtual overloads
  void identify(std::ostream& os = std::cout) const override;
  void Reset() override {}
  int isValid() const override;
  PHObject* CloneMe() const override { return new CMFlashDifferencev1(*this); }
 
  //! copy content from base class
  void CopyFrom( const CMFlashDifference& ) override;

  //! copy content from base class
  void CopyFrom( CMFlashDifference* source ) override
  { CopyFrom( *source ); }

  void setNclusters(unsigned int n) override { m_nclusters = n; }
  unsigned int getNclusters() const override { return m_nclusters; }
  //
  // difference position
  //
  float getTruthPhi() const override { return m_Phi[0]; }
  void setTruthPhi(float phi) override { m_Phi[0] = phi; }

  float getRecoPhi() const override { return m_Phi[1]; }
  void setRecoPhi(float phi) override { m_Phi[1] = phi; }

  float getTruthR() const override { return m_R[0]; }
  void setTruthR(float r) override { m_R[0] = r; }

  float getRecoR() const override { return m_R[1]; }
  void setRecoR(float r) override { m_R[1] = r; }

  float getTruthZ() const override { return m_Z[0]; }
  void setTruthZ(float z) override { m_Z[0] = z; }

  float getRecoZ() const override { return m_Z[1]; }
  void setRecoZ(float z) override { m_Z[1] = z; }


  /*
  float getTruthX() const override { return m_pos_truth[0]; }
  void setTruthX(float x) override { m_pos_truth[0] = x; }
  float getTruthY() const override { return m_pos_truth[1]; }
  void setTruthY(float y) override { m_pos_truth[1] = y; }
  float getRecoX() const override { return m_pos_reco[0]; }
  void setRecoX(float x) override { m_pos_reco[0] = x; }
  float getRecoY() const override { return m_pos_reco[1]; }
  void setRecoY(float y) override { m_pos_reco[1] = y; }
  */

 protected:
  unsigned int m_nclusters = UINT_MAX;

  float m_Phi[2] = {NAN, NAN};  // (truth, reco)
  float m_R[2] = {NAN, NAN};
  float m_Z[2] = {NAN, NAN};

  ClassDefOverride(CMFlashDifferencev1, 1)
};

#endif //TRACKBASE_CMFLASHDIFFERENCEV1_H
