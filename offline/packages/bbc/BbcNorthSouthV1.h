#ifndef __BBCNORTHSOUTHV1_H__
#define __BBCNORTHSOUTHV1_H__

#include "BbcNorthSouth.h"


class BbcNorthSouthV1 : public BbcNorthSouth
{
public:
  BbcNorthSouthV1() { }
  BbcNorthSouthV1(const Short_t npmt, const Float_t chargesum, const Float_t timing);
  virtual ~BbcNorthSouthV1() { }
  void identify(std::ostream& os = std::cout) const override;

  Short_t get_nPMT() const override { return nPmt; }
  Float_t get_nCharge() const override { return nCharge; }
  Float_t get_MeanTime() const override { return MeanTime; }

protected:
  virtual void Clear(Option_t * /*option*/ = "") override { }

  Short_t nPmt;
  Float_t nCharge;
  Float_t MeanTime;

private:
  ClassDefOverride(BbcNorthSouth,1)
};

#endif
