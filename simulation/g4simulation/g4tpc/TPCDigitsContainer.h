// TPC DIGITSCONTAINER class
// Stores TPC DIGITS per module
// note to self: it stores a digits in vector and not in map
// no need of index access since reconstruction will expand per module
// i.e. stores compressed info
// features avoidance of ctor dtor calls
// Author: Carlos Perez

#ifndef __TPCDIGITSCONTAINER_H__
#define __TPCDIGITSCONTAINER_H__

#include <map>
#include <phool/PHObject.h>
#include "TPCDigit.h"
#include "TPCDataTypes.h"

using namespace TPCDataTypes;

class TPCDigitsContainer : public PHObject {
 public:
  TPCDigitsContainer();
  virtual ~TPCDigitsContainer();
  void Add(Module_t mod, Pad_t pad, BinTime_t tim, Adc_t w, int tid=-999);
  TPCDigit* GetDigit(Module_t mod, Pad_t pad);
  void AddOrSetDigit(Module_t mod, Pad_t pad, TPCDigit* dig) {(fDigits[mod])[pad] = dig;}
  int  GetNDigits(Module_t mod) {return fNDigits[mod];}
  void Reset();

 protected:
  int fNDigits[72];
  std::map<Pad_t,TPCDigit*> fDigits[72];

  ClassDef(TPCDigitsContainer,1);
};

#endif
