// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef BBC_BBCNORTHSOUTHV1_H
#define BBC_BBCNORTHSOUTHV1_H

#include "BbcNorthSouth.h"


class BbcNorthSouthV1 : public BbcNorthSouth
{
public:
  BbcNorthSouthV1() { }
  BbcNorthSouthV1(const Short_t npmt, const float chargesum, const Float_t timing);
  virtual ~BbcNorthSouthV1() { }
  void identify(std::ostream& os = std::cout) const override;

  Short_t get_nPMT() const override { return nPmt; }
  float get_nCharge() const override { return nCharge; }
  float get_MeanTime() const override { return MeanTime; }

protected:
  virtual void Clear(Option_t * /*option*/ = "") override { }

  Short_t nPmt;
  float nCharge;
  float MeanTime;

private:
  ClassDefOverride(BbcNorthSouth,1)
};

#endif
