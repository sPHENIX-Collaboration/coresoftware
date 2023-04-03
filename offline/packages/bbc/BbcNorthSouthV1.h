// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef BBC_BBCNORTHSOUTHV1_H
#define BBC_BBCNORTHSOUTHV1_H

#include "BbcNorthSouth.h"

#include <iostream>

class BbcNorthSouthV1 : public BbcNorthSouth
{
 public:
  BbcNorthSouthV1() = delete;
  BbcNorthSouthV1(const short npmt, const float chargesum, const float timing);
  virtual ~BbcNorthSouthV1() {}
  void identify(std::ostream& os = std::cout) const override;

  short get_nPMT() const override { return nPmt; }
  float get_nCharge() const override { return nCharge; }
  float get_MeanTime() const override { return MeanTime; }

 protected:
  virtual void Clear(Option_t* /*option*/ = "") override {}

  short nPmt;
  float nCharge;
  float MeanTime;

 private:
  ClassDefOverride(BbcNorthSouth, 1)
};

#endif
