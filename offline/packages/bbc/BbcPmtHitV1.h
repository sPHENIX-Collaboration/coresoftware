// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef BBC_BBCPMTHITV1_H
#define BBC_BBCPMTHITV1_H

#include "BbcPmtHit.h"

#include <iostream>

class BbcPmtHitV1 : public BbcPmtHit
{
 public:
  BbcPmtHitV1() = delete;
  BbcPmtHitV1(const short pmt, const float adc, const float tdc0, const float tdc1);
  ~BbcPmtHitV1() override = default;

  short get_pmt() const override { return pmt; }
  float get_adc() const override { return adc; }
  float get_tdc0() const override { return tdc0; }
  float get_tdc1() const override { return tdc1; }

  void identify(std::ostream& os = std::cout) const override;

 private:
  short pmt;
  float adc;
  float tdc0;
  float tdc1;

  ClassDefOverride(BbcPmtHitV1, 1)
};

#endif
