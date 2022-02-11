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
  CMFlashDifferencev1();

  //!dtor
  ~CMFlashDifferencev1() override {}
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

  void setKey(unsigned int id) override { m_key = id; }
  unsigned int getKey() const override { return m_key; }
  //
  // difference position
  //
  float getTruthX() const override { return m_pos_truth[0]; }
  void setTruthX(float x) override { m_pos_truth[0] = x; }
  float getTruthY() const override { return m_pos_truth[1]; }
  void setTruthY(float y) override { m_pos_truth[1] = y; }

  float getRecoX() const override { return m_pos_reco[0]; }
  void setRecoX(float x) override { m_pos_reco[0] = x; }
  float getRecoY() const override { return m_pos_reco[1]; }
  void setRecoY(float y) override { m_pos_reco[1] = y; }


 protected:
  unsigned int m_key;
  float m_pos_truth[2];              //< truth x,y
  float m_pos_reco[2];               //< reco x,y


  ClassDefOverride(CMFlashDifferencev1, 1)
};

#endif //TRACKBASE_CMFLASHDIFFERENCEV1_H
