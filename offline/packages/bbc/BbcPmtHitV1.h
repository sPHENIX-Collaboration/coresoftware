// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef BBC_BBCPMTHITV1_H
#define BBC_BBCPMTHITV1_H

#include "BbcPmtHit.h"

#include <iostream>

class BbcPmtHitV1 : public BbcPmtHit
{

public:
  BbcPmtHitV1() { }
  BbcPmtHitV1(const Short_t pmt, const Float_t adc, const Float_t tdc0, const Float_t tdc1);
  ~BbcPmtHitV1() override = default;

  Short_t get_pmt()  const override {return pmt;}
  Float_t get_adc()  const override {return adc;}
  Float_t get_tdc0() const override {return tdc0;}
  Float_t get_tdc1() const override {return tdc1;}

  void identify(std::ostream& os = std::cout) const override;

private:
  Short_t pmt;
  Float_t adc;
  Float_t tdc0;
  Float_t tdc1;

  ClassDefOverride(BbcPmtHitV1,1)
};

#endif
