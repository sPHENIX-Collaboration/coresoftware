// TPC DIGIT class
// Stores vector of digits associated to one readout pad
// Author: Carlos Perez
#ifndef __TPCDIGIT_H__
#define __TPCDIGIT_H__

#include <map>
#include "TPCDataTypes.h"
using namespace TPCDataTypes;

class TPCDigit {
 public:
  TPCDigit();
  virtual ~TPCDigit();
  void AddContent(const Time_t bin, Adc_t v);
  Adc_t GetContent(const Time_t bin);
  Int_t GetEntries() {return static_cast<Int_t>(fTrainOfDigits.size());};
  void Reset() {fTrainOfDigits.clear();}
 protected:
  std::map<Time_t,Adc_t> fTrainOfDigits; // train of adc

  ClassDef(TPCDigit,1);
};

#endif
