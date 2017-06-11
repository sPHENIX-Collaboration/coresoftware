// TPC DIGIT class
// Stores vector of digits associated to one readout pad
// Author: Carlos Perez
#include <map>
#include "TPCDigit.h"
#include "TPCDataTypes.h"

using namespace TPCDataTypes;

//=====
TPCDigit::TPCDigit() {}
//=====
TPCDigit::~TPCDigit() {
  fTrainOfDigits.clear();
}
//=====
Adc_t TPCDigit::GetContent(const Time_t bin) {
  std::map<Time_t,Adc_t>::const_iterator iter = fTrainOfDigits.find(bin);
  if(iter==fTrainOfDigits.end()) return 0;
  return (*iter).second;
}
//=====
void TPCDigit::AddContent(const Time_t bin, Adc_t val) {
  std::map<Time_t,Adc_t>::iterator iter = fTrainOfDigits.find(bin);
  if(iter==fTrainOfDigits.end()) fTrainOfDigits[bin] = val;
  else (*iter).second += val;
}
