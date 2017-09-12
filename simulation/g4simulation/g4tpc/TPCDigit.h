// TPC DIGIT class
// Stores vector of digits associated to one readout pad
// Author: Carlos Perez
#ifndef __TPCDIGIT_H__
#define __TPCDIGIT_H__

#include <map>
#include <g4detectors/PHG4Cell.h>
#include "TPCDataTypes.h"

using namespace TPCDataTypes;

class TPCDigit {
 public:
  TPCDigit();
  virtual ~TPCDigit();
  void AddContent(const BinTime_t bin, Adc_t v, int Tid=-999);
  Adc_t GetContent(const BinTime_t bin);
  int   GetEntries();
  void SetCell(PHG4Cell *cel) {fCell=cel;}
  void Reset();
 protected:
  PHG4Cell *fCell;

  ClassDef(TPCDigit,1);
};

#endif
