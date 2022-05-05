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
  // difference position
  //
 virtual float getTruthPhi() const { return NAN; }
  virtual  void setTruthPhi(float) {}

  virtual float getRecoPhi() const { return NAN; }
  virtual void setRecoPhi(float) {}

  virtual float getTruthR() const { return NAN; }
  virtual void setTruthR(float) {}

  virtual float getRecoR() const { return NAN; }
  virtual void setRecoR(float) {}

  virtual float getTruthZ() const { return NAN; }
  virtual void setTruthZ(float) {}

  virtual float getRecoZ() const { return NAN; }
  virtual void setRecoZ(float) {}

  virtual unsigned int getNclusters() const { return UINT_MAX; }
  virtual void setNclusters(unsigned int) {}

 protected:
  CMFlashDifference() = default;
  ClassDefOverride(CMFlashDifference, 1)
};

#endif //TRACKBASE_CMFLASHDIFFERENCE_H
