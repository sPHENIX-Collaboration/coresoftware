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

  short get_nPMT() const override { return nPmt; }
  float get_nCharge() const override { return nCharge; }
  float get_MeanTime() const override { return MeanTime; }

 protected:
  virtual void Clear(Option_t* /*option*/ = "") override {}

  short nPmt = 0;
  float nCharge = std::numeric_limits<float>::quiet_NaN();
  float MeanTime = std::numeric_limits<float>::quiet_NaN();

 private:
  ClassDefOverride(BbcNorthSouth, 1)
};

#endif
