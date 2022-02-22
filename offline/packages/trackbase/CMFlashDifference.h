/**
 * @file trackbase/CMFlashDifference.h
 * @author Tony Frawley
 * @date January 2022
 * @brief Base class for central membrane flash cluster difference due to distortions 
 */
#ifndef TRACKBASE_CMFLASHDIFFERENCE_H
#define TRACKBASE_CMFLASHDIFFERENCE_H

#include <phool/PHObject.h>

#include <climits>
#include <cmath>
#include <iostream>
#include <memory>


/**
 * @brief Base class for Central Membrane flassh combined difference object
 *
 * Virtual base class for TPC CM flash difference object 
 */
class CMFlashDifference : public PHObject
{
 public:
  //! dtor
  ~CMFlashDifference() override {}
  // PHObject virtual overloads
  void identify(std::ostream& os = std::cout) const override
  {
    os << "CMFlashDifference base class" << std::endl;
  }
  void Reset() override {}
  int isValid() const override { return 0; }

  
  //! import PHObject CopyFrom, in order to avoid clang warning
  using PHObject::CopyFrom;
  
  //! copy content from base class
  virtual void CopyFrom( const CMFlashDifference& )  {}

  //! copy content from base class
  virtual void CopyFrom( CMFlashDifference* )  {}

  //
  // difference id
  //
  virtual unsigned int getKey() const { return UINT_MAX; }
  virtual void setKey(unsigned int) {}
 
  //
  // difference position
  //
  virtual float getTruthX() const { return NAN; }
  virtual void setTruthX(float) {}
  virtual float getTruthY() const { return NAN; }
  virtual void setTruthY(float) {}

  virtual float getRecoX() const { return NAN; }
  virtual void setRecoX(float) {}
  virtual float getRecoY() const { return NAN; }
  virtual void setRecoY(float) {}

 protected:
  CMFlashDifference() = default;
  ClassDefOverride(CMFlashDifference, 1)
};

#endif //TRACKBASE_CMFLASHDIFFERENCE_H
