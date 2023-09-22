// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef BBC_BBCNORTHSOUTHV1_H
#define BBC_BBCNORTHSOUTHV1_H

#include "BbcNorthSouth.h"

#include <iostream>
#include <limits>

class BbcNorthSouthV1 : public BbcNorthSouth
{
 public:
  BbcNorthSouthV1() = default;
  BbcNorthSouthV1(const short npmt, const float chargesum, const float timing);
  virtual ~BbcNorthSouthV1() {}
  void identify(std::ostream& os = std::cout) const override;

  short get_nPMT() const override { return bn; }
  float get_nCharge() const override { return bq; }
  float get_MeanTime() const override { return bt; }

 protected:
  virtual void Clear(Option_t* /*option*/ = "") override {}

  short bn = 0;
  float bq = std::numeric_limits<float>::quiet_NaN();
  float bt = std::numeric_limits<float>::quiet_NaN();

 private:
  ClassDefOverride(BbcNorthSouth, 1)
};

#endif
